# Porthmeus
# 21.01.22

# take the best explaining data set and prepare result tables for the presence/absence of reactions to HB/Mayo scores


# libraries
## @knitr loadLibraries
require(data.table)
require(lme4)
require(lmerTest)
require(caret)
require(foreach)
require(parallel)
require(doMC)
require(fpc)
require(car)
require(cowplot)


# load data
## @knitr loadData_stat
rxnCount <- read.csv("../../../results/ModelAnalysis/rxnIdMat.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)

clinic <- fread("../../../resources/ClinicalData.csv")
meta <- fread("../../../resources/META_emed_future.csv")
blcklst <- fread("../../../resources/sampleBlacklist.csv")

# make rxnCounts a integer matrix
rxnCount2 <- rxnCount
rxnCount2[,] <- 0
rxnCount2[rxnCount == "True"] <- 1
rxnCount <- rxnCount2

# merge meta and clincal data
meta[,PatientID_Time := paste(PatientID, Time_seq, sep ="_")]
meta <- meta[Time_seq != -1,]
meta <- merge(meta,clinic, by="PatientID_Time", all.x=TRUE, suffix = c(".meta",""))


# remove degraded samples
blckSmpls <- blcklst[relevance < 2, Sample]
meta <- meta[!(SeqID %in% blckSmpls),]
rxnCount <- rxnCount[,!(colnames(rxnCount) %in% blckSmpls)]
meta <- meta[SeqID %in% colnames(rxnCount),]

# remove near zero variance reactions
#near0var_idx <- nearZeroVar(t(rxnCount))
#rxnCount <- rxnCount[-near0var_idx,]
#
## cluster the rxns to reduce the dimensionality and unecessary tests
#corre <- cor(t(rxnCount))
#corre <- 1-sqrt(corre^2)
#
#set.seed(23501839)
#corre.cluster <- dbscan(corre, method = "dist", eps = 0.1, MinPts = 3)
#
## create a cluster table to export and read later on
#tbl.cluster <- data.table( cluster = corre.cluster[["cluster"]],
#                          eps = corre.cluster[["eps"]],
#                          MinPts = corre.cluster[["MinPts"]],
#                          isseed = corre.cluster[["isseed"]],
#                          rxn = rownames(corre),
#                          ClusterID = paste0("cl",corre.cluster[["cluster"]]))
#
## get a representative reaction for each of the clusters
#tbl.cluster[cluster == 0 ,rep.rxn := rxn]
#tbl.cluster[cluster != 0, rep.rxn := names(which.min(rowSums(corre[rxn,rxn]))), by = cluster]
#write.csv(tbl.cluster, file = "../results/PA_DBSCAN_Cluster.csv")

tbl.cluster <- fread("../results/PA_DBSCAN_Cluster.csv")

rxnCount2 <- rxnCount[unique(tbl.cluster[, rep.rxn]),]

# sort meta data to fit the samples
meta <- data.frame(meta)
rownames(meta) <- meta[,"SeqID"]
meta <- meta[colnames(rxnCount2),]


# create a list of lmer models
# for parallel computation register threads
threads <- detectCores()-1
registerDoMC(threads)

models <- list()
for(diag in unique(meta[,"Diagnosis"])){
    mods <- foreach(rxn = rownames(rxnCount2)) %dopar% {
        dat <- data.table(HB_Mayo_impu = meta[["HB_Mayo_impu"]],
                          Diagnosis = meta[["Diagnosis"]],
                          rxnCount = as.factor(unlist(rxnCount2[rxn,])),
                          PatientID = meta[["PatientID"]])
        dat <- dat[Diagnosis == diag,]
        lmer(data = dat,
             HB_Mayo_impu~ rxnCount+(1|PatientID))
    }
    names(mods) <- rownames(rxnCount2)
    models[[diag]] <- mods
}
saveRDS(models, file = "../results/PA.LMM.HBMayo.etiology.DiagnosisIndividual_models.RDS")

#names(mods) <- rownames(rxnCount)[-linCombos[["remove"]]]

# extract the relevant information from the list
# generate a vector with the information what diagnosis set it was:
diag <- unlist(lapply(names(models), function(x) rep(x,length(models[[x]]))))
mods <- unlist(models)
sums <- lapply(mods, function(x) coef(summary(x))["rxnCount1",])
sums <- do.call(rbind, sums)
sums <- data.table(sums, keep.rownames = TRUE)

# adjust p.value
sums[,padj := p.adjust(`Pr(>|t|)`, method = "BH")]
colnames(sums)[1] <- "diag.rxn"
sums[,Diagnosis := gsub("(.*)\\..*","\\1",diag.rxn)]
sums[,rxn := gsub(".*\\.(.*)","\\1",diag.rxn)]


# save the results
gz <- gzfile(description = "../results/PA.LMM.HBMayo.etiology.DiagnosisIndividual_coefs.csv", "w")
write.csv(sums, file = gz, row.names = FALSE)
close(gz)

#save diagnostic plots for all the significant reactions
rxns <- sums[padj < 0.05, diag.rxn]
pbar <- txtProgressBar(min = 0 , max = length(rxns), style = 3)
pdf("../results/PA.LMM.HBMayo.etiology.DiagnosisIndividual_sigDiagPlots.pdf", width = 8, height= 5)
i <- 0 
for(rxn in rxns){
    print(plot_grid(~qqPlot(resid(mods[[rxn]]), main = rxn), plot(mods[[rxn]])))
    i <- i+1
    setTxtProgressBar(pbar, i)
}
dev.off()

