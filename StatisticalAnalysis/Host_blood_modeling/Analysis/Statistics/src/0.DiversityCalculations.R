# Porthmeus
# 10.03.22

# create some diversity measures for the different analysis values we have

require(data.table)
require(vegan)
require(lme4)
require(lmerTest)

# load data
rxnCount <- read.csv("../../../results/ModelAnalysis/rxnIdMat.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
rxnExpr <- read.csv("../../../results/RxnExpression/rxnExpr.GL10-L50-GU90.colormore3D.csv", row.names = 1)
# load FVA results
FVA.min <- read.csv("../../../results/FVA/minFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
FVA.max <- read.csv("../../../results/FVA/maxFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
# convert the values to ranges and centers
FVA.range <- FVA.max-FVA.min
FVA.range[FVA.range < 1E-6] <- 0


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
meta <- meta[SeqID %in% colnames(rxnCount),]

rxnExpr <- rxnExpr[,!(colnames(rxnExpr) %in% blckSmpls)]
rxnCount <- rxnCount[,!(colnames(rxnCount) %in% blckSmpls)]
sel<-!(colnames(FVA.min) %in% blckSmpls)
FVA.range <- FVA.range[,sel]



# simply count the reactions retrieved per sample
shannon.rxnExpr <- diversity(t(rxnExpr), index = "shannon")
richness.PA <- colSums(rxnCount)/nrow(rxnCount)
shannon.FVA <-diversity(t(FVA.range), index = "shannon")
 
meta[, richness.PA:= richness.PA[SeqID]]
meta[, shannon.rxnExpr:= shannon.rxnExpr[SeqID]]
meta[, shannon.FVA:= shannon.FVA[SeqID]]

# save the meta data with the diverstiy measures
write.csv(meta, file = "../results/diversityMeasures.csv", row.names = FALSE)


#psych::pairs.panels(meta[,.(HB_Mayo_impu, richness.PA, shannon.rxnExpr, shannon.FVA)],method = "spearman", stars = TRUE)
#
#mod <- lmer(data= meta,
#            HB_Mayo_impu ~ richness.PA + (1|PatientID))
#summary(mod)
#mod <- lmer(data= meta,
#            HB_Mayo_impu ~ shannon.rxnExpr + (1|PatientID))
#summary(mod)
#mod <- lmer(data= meta,
#            HB_Mayo_impu ~ shannon.FVA + (1|PatientID))
#summary(mod)
