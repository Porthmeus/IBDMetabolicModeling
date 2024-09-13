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



rxnCount <- read.csv("../../../results/ModelAnalysis/rxnIdMat.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)

# make rxnCounts a integer matrix
rxnCount2 <- rxnCount
rxnCount2[,] <- 0
rxnCount2[rxnCount == "True"] <- 1
rxnCount <- rxnCount2

clinic <- fread("../../../resources/ClinicalData.csv")
meta <- fread("../../../resources/META_emed_future.csv")
blcklst <- fread("../../../resources/sampleBlacklist.csv")

# merge meta and clincal data
meta[,PatientID_Time := paste(PatientID, Time_seq, sep ="_")]
meta <- meta[Time_seq != -1,]
meta <- merge(meta,clinic, by="PatientID_Time", all.x=TRUE, suffix = c(".meta",""))

# remove degraded samples
blckSmpls <- blcklst[relevance < 2, Sample]
meta <- meta[!(SeqID %in% blckSmpls),]
rxnCount <- rxnCount[,!(colnames(rxnCount) %in% blckSmpls)]
meta <- meta[SeqID %in% colnames(rxnCount),]

# remove control samples
meta <- meta[Response != "C",]
rxnCount <- rxnCount[, meta[,SeqID]]

# load the clustering 
cluster <- fread("../results/PA_DBSCAN_Cluster.csv")
cluster[,rxn := rep.rxn]


# remove rxns which were clustered
rxns <- unique(cluster[,rxn])
rxnCount <- rxnCount[rxns,]


# melt the FVA data to one data.table to make the associations
rxnCount.melt <- melt(
                       data.table(
                                  rxnCount,
                                  keep.rownames = TRUE
                                  ),
                       id.vars = "rn",
                       variable.name = "SeqID",
                       value.name = "value"
                       )

# merge the meta data
rxnCount.meta <- merge(rxnCount.melt, meta, by="SeqID")


# fit the model
threads <- detectCores()-1
registerDoMC(threads)

mods <- foreach(reac = rxns) %dopar% {
#for(reac in rxns){
    dat <- rxnCount.meta[rn == reac & Time_seq <= 14,]
    mod <- tryCatch({glmer(data = dat,
          formula = factor(Response) ~value + (1|PatientID),
          family = "binomial")},
                    error = function(e){
                        print(paste0("Error to build the model for ",reac))
                        print(e)
                        return(NA)})
}
names(mods) <- rxns
saveRDS(mods, file = "../results/PA.GLMM.Response.biomarker_models.RDS")

# get the stats
stats <- list() 
for(rxn in names(mods)){
    mod <- mods[[rxn]]
    if(!is.na(mod)){
        # get all the interesting coefficients
        dat <- as.data.frame(t(coef(summary(mod))["value",]))
        stats[[rxn]] <- cbind(
                              rxn = rxn, # add the reaction name to the table
                              coef = "rxnCount",
                              dat
                              )
    }
}
stats <- do.call(rbind, stats)
stats <- data.table(stats)
stats[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
stats[,OR := round(exp(Estimate),3)]



# save the models and the table
write.csv(stats, row.names=FALSE, file = "../results/PA.GLMM.Response.biomarker_coefs.csv")

# create diagnostic plots for the significant reactions
rxns.sig <- stats[padj < 0.05, rxn]
pb <- txtProgressBar(min = 0, max = length(rxns), style = 3)
i<-0
pdf("../results/PA.GLMM.Response.biomarker_sigDiagPlots.pdf")
for(rxn in rxns.sig){
    sim <- simulateResiduals(mods[[rxn]])
    plot(sim)
    mtext(rxn)
    i<-1+i
    setTxtProgressBar(pb,i)
}
dev.off()
