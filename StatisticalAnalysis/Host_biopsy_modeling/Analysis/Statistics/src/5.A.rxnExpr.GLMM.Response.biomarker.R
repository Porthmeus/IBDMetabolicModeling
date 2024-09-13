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



rxnExpr <- read.csv("../../../results/RxnExpression/rxnExpr.GL10-L50-GU90.colormore3D.csv", row.names = 1)

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
rxnExpr <- rxnExpr[,!(colnames(rxnExpr) %in% blckSmpls)]
meta <- meta[SeqID %in% colnames(rxnExpr),]

# remove control samples
meta <- meta[Response != "C",]
rxnExpr <- rxnExpr[, meta[,SeqID]]

# load the clustering 
cluster <- fread("../results/rxnExpr_DBSCAN_Cluster.csv")
cluster[,rxn := rep.rxn]


# remove rxns which were clustered
rxns <- unique(cluster[,rxn])
rxnExpr <- rxnExpr[rxns,]

# scale rxnExpr
rxnExpr <- t(scale(t(rxnExpr)))


# melt the FVA data to one data.table to make the associations
rxnExpr.melt <- melt(
                       data.table(
                                  rxnExpr,
                                  keep.rownames = TRUE
                                  ),
                       id.vars = "rn",
                       variable.name = "SeqID",
                       value.name = "value"
                       )

# merge the meta data
rxnExpr.meta <- merge(rxnExpr.melt, meta, by="SeqID")


# fit the model
threads <- detectCores()-1
registerDoMC(threads)

mods <- foreach(reac = rxns) %dopar% {
#for(reac in rxns){
    dat <- rxnExpr.meta[rn == reac & Time_seq <= 14,]
    mod <- tryCatch({glmer(data = dat,
          formula = factor(Response) ~ value + (1|PatientID),
          family = "binomial")},
                    error = function(e){
                        print(paste0("Error to build the model for ",reac))
                        print(e)
                        return(NA)})
}
names(mods) <- rxns

# get the stats
stats <- list() 
for(rxn in names(mods)){
    mod <- mods[[rxn]]
    if(!is.na(mod)){
        # get all the interesting coefficients
        dat <- as.data.frame(t(coef(summary(mod))["value",]))
        stats[[rxn]] <- cbind(
                              rxn = rxn, # add the reaction name to the table
                              coef = "rxnExpr",
                              dat
                              )
    }
}
stats <- do.call(rbind, stats)
stats <- data.table(stats)
stats[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
stats[,OR := round(exp(Estimate),3)]



# save the models and the table
write.csv(stats, row.names=FALSE, file = "../results/rxnExpr.GLMM.Response.biomarker_coefs.csv")
saveRDS(mods, file = "../results/rxnExpr.GLMM.Response.biomarker_models.RDS")

# create diagnostic plots for the significant reactions
rxns.sig <- stats[padj < 0.05, rxn]
pb <- txtProgressBar(min = 0, max = length(rxns), style = 3)
i<-0
pdf("../results/rxnExpr.GLMM.Response.biomarker_sigDiagPlots.pdf")
for(rxn in rxns.sig){
    sim <- simulateResiduals(mods[[rxn]])
    plot(sim)
    mtext(rxn)
    i<-1+i
    setTxtProgressBar(pb,i)
}
dev.off()
