# Porthmeus
# 21.02.22

# libraries
require(data.table)
require(lme4)
require(lmerTest)
require(caret)
require(foreach)
require(parallel)
require(doParallel)
require(fpc)
require(car)
require(cowplot)
require(pbapply)


# load metadata
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


# load FVA results
FVA.min <- read.csv("../../../results/FVA/minFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
FVA.max <- read.csv("../../../results/FVA/maxFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
# convert the values to ranges and centers
FVA.range <- FVA.max-FVA.min
FVA.center <- (FVA.max+FVA.min)/2

sel<-!(colnames(FVA.min) %in% blckSmpls)
FVA.range <- FVA.range[,sel]
FVA.center <- FVA.center[,sel]

# remove nearZeroVariance rxns
NZV.range <- nearZeroVar(t(FVA.range))
NZV.center <- nearZeroVar(t(FVA.center))
NZV <- intersect(NZV.range, NZV.center)
FVA.range <- FVA.range[-NZV,]
FVA.center <- FVA.center[-NZV,]

# scale the data
FVA.range <- t(scale(t(FVA.range)))
FVA.center <- t(scale(t(FVA.center)))


# correlate the values of each sample to each other for both matrices
FVA.range.corr <- cor(t(FVA.range))
FVA.center.corr <- cor(t(FVA.center))
# replace NaN's with the value of the other correlation matrix 
FVA.range.corr[is.na(FVA.range.corr)] <- FVA.center.corr[is.na(FVA.range.corr)]
FVA.center.corr[is.na(FVA.center.corr)] <- FVA.range.corr[is.na(FVA.center.corr)]
# get the sign of the correlations
FVA.center.corr.sign <- melt(data.table(sign(FVA.center.corr),
                                        keep.rownames = TRUE),
                             id.var = "rn",
                             variable.name = "rep.rxn",
                             value.name = "sign.cor")
colnames(FVA.center.corr.sign)[1] <- "rxn" 
FVA.range.corr.sign <- melt(data.table(sign(FVA.range.corr),
                                        keep.rownames = TRUE),
                             id.var = "rn",
                             variable.name = "rep.rxn",
                             value.name = "sign.cor")
colnames(FVA.range.corr.sign)[1] <- "rxn" 
# transform correlation to distance measure
FVA.range.corr <- 1-sqrt(FVA.range.corr^2)
FVA.center.corr <- 1-sqrt(FVA.center.corr^2)
# take the mean of both correlation matrices
FVA.corr <- (FVA.center.corr+FVA.range.corr)/2

# cluster the reactions with dbscan
set.seed(6218732)
cluster <- dbscan(FVA.corr,method = "dist", eps = 0.1, MinPts = 3)

# create a reasonable table to export
tbl.cluster <- data.table( cluster = cluster[["cluster"]],
                          eps = cluster[["eps"]],
                          MinPts = cluster[["MinPts"]],
                          isseed = cluster[["isseed"]],
                          rxn = rownames(FVA.corr),
                          clusterID = paste0("cl",cluster[["cluster"]]))

# find the most representative reaction in each cluster
# to do so find the reaction in the cluster with the smallest correlation distance to all other clusters
tbl.cluster[cluster == 0, rep.rxn := rxn]
tbl.cluster[cluster != 0, rep.rxn := names(which.min(rowSums(FVA.corr[rxn,rxn]))), by = cluster]

# add the signs of the correlation again to the cluster table
tbl.cluster <- merge(tbl.cluster, FVA.range.corr.sign)
tbl.cluster <- merge(tbl.cluster, FVA.center.corr.sign, suffix = c(".range",".center"))

# save cluster table
write.csv(tbl.cluster, file = "../results/FVA_DBSCAN_Cluster.csv")

### PA 
# load data
## @knitr loadData_stat
rxnCount <- read.csv("../../../results/ModelAnalysis/rxnIdMat.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)

# make rxnCounts a integer matrix
rxnCount2 <- rxnCount
rxnCount2[,] <- 0
rxnCount2[rxnCount == "True"] <- 1
rxnCount <- rxnCount2

rxnCount <- rxnCount[,!(colnames(rxnCount) %in% blckSmpls)]
near0var_idx <- nearZeroVar(t(rxnCount))
rxnCount <- rxnCount[-near0var_idx,]

# cluster the rxns to reduce the dimensionality and unecessary tests
corre <- cor(t(rxnCount))
# save the signs of the correlation
corr.sign <- melt(data.table(sign(corre),
                                        keep.rownames = TRUE),
                             id.var = "rn",
                             variable.name = "rep.rxn",
                             value.name = "sign.cor")
colnames(corr.sign)[1] <- "rxn" 
corre <- 1-sqrt(corre^2)

set.seed(23501839)
corre.cluster <- dbscan(corre, method = "dist", eps = 0.1, MinPts = 3)

# create a cluster table to export and read later on
tbl.cluster <- data.table( cluster = corre.cluster[["cluster"]],
                          eps = corre.cluster[["eps"]],
                          MinPts = corre.cluster[["MinPts"]],
                          isseed = corre.cluster[["isseed"]],
                          rxn = rownames(corre),
                          ClusterID = paste0("cl",corre.cluster[["cluster"]]))

# get a representative reaction for each of the clusters
tbl.cluster[cluster == 0 ,rep.rxn := rxn]
tbl.cluster[cluster != 0, rep.rxn := names(which.min(rowSums(corre[rxn,rxn]))), by = cluster]

# add back the sign of the correlation
tbl.cluster <- merge(tbl.cluster, corr.sign)

write.csv(tbl.cluster, file = "../results/PA_DBSCAN_Cluster.csv")


### rxnExpression
# load data
rxnExpr <- read.csv("../../../results/RxnExpression/rxnExpr.GL10-L50-GU90.colormore3D.csv", row.names = 1)
rxnExpr <- rxnExpr[,!(colnames(rxnExpr) %in% blckSmpls)]

# remove all reactions where no information of expression is present
rxnExpr <- rxnExpr[rowSums(rxnExpr) != 0,]

# remove all nearZeroVariables
NZV <- nearZeroVar(t(rxnExpr))
rxnExpr <- rxnExpr[-NZV,]

# scale the data
rxnExpr <- t(scale(t(rxnExpr)))

# cluster the reactions by correlation
mat.cor <- cor(t(rxnExpr))
# save the sign of correlation for later
corr.sign <- melt(data.table(sign(mat.cor),
                                        keep.rownames = TRUE),
                             id.var = "rn",
                             variable.name = "rep.rxn",
                             value.name = "sign.cor")
colnames(corr.sign)[1] <- "rxn" 

mat.cor <- 1-sqrt(mat.cor^2)

set.seed(6841) # change in your code
cluster <- dbscan(mat.cor, # the matrix with the correlation coefficients converted to distance
                  method = "dist", # telling dbscan not to calculate the distances by itself
                  eps = 0.1, # the threshold how distant members are allowed to be to each other to count a cluster
                  MinPts = 3) # the number of how many reactions have to make up a cluster at least - 3 is the minimum

# create a reasonable table to export
tbl.cluster <- data.table( cluster = cluster[["cluster"]],
                          eps = cluster[["eps"]],
                          MinPts = cluster[["MinPts"]],
                          isseed = cluster[["isseed"]],
                          rxn = rownames(mat.cor))

# get a representative reaction for each cluster
tbl.cluster[cluster == 0 ,rep.rxn := rxn]
tbl.cluster[cluster != 0, rep.rxn := names(which.min(rowSums(mat.cor[rxn,rxn]))), by = cluster]

# add back the sign of the correlation
tbl.cluster <- merge(tbl.cluster, corr.sign)

write.csv(tbl.cluster, file = "../results/rxnExpr_DBSCAN_Cluster.csv")
