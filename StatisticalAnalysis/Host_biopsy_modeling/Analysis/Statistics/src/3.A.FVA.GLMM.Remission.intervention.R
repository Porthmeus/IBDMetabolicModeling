# Porthmeus
# 17.03.22


# required packages
require(data.table)
require(lme4)
require(lmerTest)
require(foreach)
require(parallel)
require(doMC)
require(DHARMa)



# load FVA results
FVA.min <- read.csv("../../../results/FVA/minFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
FVA.max <- read.csv("../../../results/FVA/maxFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
# convert the values to ranges and centers
FVA.range <- FVA.max-FVA.min
FVA.center <- (FVA.max+FVA.min)/2

# load metadata
clinic <- fread("../../../resources/ClinicalData.csv")
meta <- fread("../../../resources/META_emed_future.csv")
blcklst <- fread("../../../resources/sampleBlacklist.csv")

# merge meta and clincal data
meta[,PatientID_Time := paste(PatientID, Time_seq, sep ="_")]
meta <- meta[Time_seq != -1,]
meta <- merge(meta,clinic, by="PatientID_Time", all.x=TRUE, suffix = c(".meta",""))

# remove controls
meta <- meta[Remission != "C",]

# remove degraded samples
blckSmpls <- blcklst[relevance < 2, Sample]
meta <- meta[!(SeqID %in% blckSmpls),]
meta <- meta[SeqID %in% colnames(FVA.min),]
sel<-!(colnames(FVA.min) %in% blckSmpls)
FVA.range <- FVA.range[,sel]
FVA.center <- FVA.center[,sel]

# load the clustering 
cluster <- fread("../results/FVA_DBSCAN_Cluster.csv")
cluster[,rxn := rep.rxn]


# remove rxns which were clustered
rxns <- unique(cluster[,rxn])
FVA.range <- FVA.range[rxns,]
FVA.center <- FVA.center[rxns,]

# scale FVA
FVA.range <- t(scale(t(FVA.range)))
FVA.center <- t(scale(t(FVA.center)))


# melt the FVA data to one data.table to make the associations
FVA.range.melt <- melt(
                       data.table(
                                  FVA.range,
                                  keep.rownames = TRUE
                                  ),
                       id.vars = "rn",
                       variable.name = "SeqID",
                       value.name = "range"
                       )
FVA.center.melt <- melt(
                       data.table(
                                  FVA.center,
                                  keep.rownames = TRUE
                                  ),
                       id.vars = "rn",
                       variable.name = "SeqID",
                       value.name = "center"
                       )
FVA <- merge(FVA.center.melt, FVA.range.melt)

# merge the meta data
FVA <- merge(FVA, meta, by="SeqID")

# load the other coef table to construct the formula
sigs.FVA <- fread("../results/FVA.LMM.HBMayo.etiology_coefs.csv")

# fit the model
threads <- detectCores()-1
registerDoMC(threads)

mods <- foreach(reac = rxns) %dopar% {
#for(reac in rxns[1:9]){
    # remove the last time point, because this encodes the remission state in the HB/Mayo setup
    dat <- FVA[rn == reac & Time_seq != max(Time_seq),]
    form <- as.formula(paste0("factor(Remission) ~ HB_Mayo_impu +",paste(sigs.FVA[rxn == reac, coef], collapse = "+"), "+ (1|PatientID)"))
    mod <- tryCatch({glmer(data = dat,
          formula = form,
          family = "binomial")},
                    error = function(e){
                        print(paste0("Error to build the model for ",reac))
                        print(e)
                        return(NA)})
}
names(mods) <- rxns
saveRDS(mods, file = "../results/FVA.GLMM.Remission.intervention_models.RDS")

# get the stats
stats <- list() 
for(rxn in names(mods)){
    mod <- mods[[rxn]]
    if(!is.na(mod)){
        # get all the interesting coefficients
        form_split <- colnames(coef(mod)[[1]])[-(1:2)]
        if(length(form_split) >1){
            dat <- as.data.frame(coef(summary(mod))[form_split,])
        } else {
            dat <- as.data.frame(t(coef(summary(mod))[form_split,]))
        }
        stats[[rxn]] <- cbind(
                              rxn = rxn, # add the reaction name to the table
                              coef = form_split,
                              dat
                              )
    }
}
stats <- do.call(rbind, stats)
stats <- data.table(stats)
stats[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
stats[,OR := round(exp(Estimate),3)]



# save the models and the table
write.csv(stats, row.names=FALSE, file = "../results/FVA.GLMM.Remission.intervention_coefs.csv")


# create diagnostic plots for the significant reactions
rxns.sig <- unique(stats[padj < 0.05, rxn])
pdf("../results/FVA.GLMM.Remission.intervention_sigDiagPlots.pdf")
for(rxn in rxns.sig){
    sim <- simulateResiduals(mods[[rxn]])
    plot(sim)
    mtext(rxn)
}
dev.off()
