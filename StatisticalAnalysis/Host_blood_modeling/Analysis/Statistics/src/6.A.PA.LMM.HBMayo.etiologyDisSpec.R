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

# remove all reactions where no information of expression is present
rxnCount <- rxnCount[rowSums(rxnCount) != 0,]

# remove all nearZeroVariables
NZV <- nearZeroVar(t(rxnCount))
rxnCount <- rxnCount[-NZV,]

# scale the data
rxnCount <- t(scale(t(rxnCount)))

tbl.cluster<- fread(file = "../results/PA_DBSCAN_Cluster.csv")

rxnCount2 <- rxnCount[unique(tbl.cluster[, rep.rxn]),]

# sort meta data to fit the samples
meta <- data.frame(meta)
rownames(meta) <- meta[,"SeqID"]
meta <- meta[colnames(rxnCount2),]


# create a list of lmer models
# for parallel computation register threads
threads <- detectCores()-1
registerDoMC(threads)

mods <- foreach(rxn = rownames(rxnCount2)) %dopar% {
    dat <- data.frame(HB_Mayo_impu = meta[["HB_Mayo_impu"]],
                      rxnCount = as.factor(unlist(rxnCount2[rxn,])),
                      PatientID = meta[["PatientID"]])
    lmer(data = dat,
         HB_Mayo_impu~ rxnCount+(1|PatientID))
}
names(mods) <- rownames(rxnCount2)
saveRDS(mods, file = "../results/PA.LMM.HBMayo.etiologyDisSpec_models.RDS")
#names(mods) <- rownames(rxnCount)[-linCombos[["remove"]]]

# extract the relevant information from the list
sums <- lapply(mods, function(x) coef(summary(x))["rxnCount1",])
sums <- do.call(rbind, sums)
sums <- data.table(sums, keep.rownames = TRUE)

# adjust p.value
sums[,padj := p.adjust(`Pr(>|t|)`, method = "BH")]
colnames(sums)[1] <- "rxn"

# save the results
write.csv(sums, file = "../results/PA.LMM.HBMayo.etiologyDisSpec_coefs.csv", row.names = FALSE)

#save diagnostic plots for all the significant reactions
rxns <- sums[padj < 0.05, rxn]
pbar <- txtProgressBar(min = 0 , max = length(rxns), style = 3)
pdf("../results/PA.LMM.HBMayo.etiologyDisSpec_sigDiagPlots.pdf", width = 8, height= 5)
i <- 0 
for(rxn in rxns){
    print(plot_grid(~qqPlot(resid(mods[[rxn]]), main = rxn), plot(mods[[rxn]])))
    i <- i+1
    setTxtProgressBar(pbar, i)
}
dev.off()

