# Porthmeus
# 08.03.22

# use the mapping of the TPM values of the expressed genes to draw conclusions on the usage of the metabolism

# load libraries
require(data.table)
require(lme4)
require(lmerTest)
require(caret)
require(foreach)
require(parallel)
require(fpc)
require(car)
require(cowplot)
require(pbapply)

# load data
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


# melt the data and merge it with the meta data
rxnExpr.melt <- melt(data.table(rxnExpr[unique(tbl.cluster[,rep.rxn]),], keep.rownames = TRUE),
                     id.vars = "rn",
                     variable.name = "SeqID",
                     value.name = "Expression")

rxnExpr.melt <- merge(rxnExpr.melt, meta, by = "SeqID")
rxnExpr.melt <- rxnExpr.melt[Remission != "C",]

# fit the models

createMods <-function(rxn){
    #foreach(rxn = rownames(rxnCount2)) %dopar% {
    dat.mod <- rxnExpr.melt[rn == rxn,]
    mod <- tryCatch(lmer(data = dat.mod, Expression ~ Remission*Time_seq+(1|PatientID)),
                    error = function(e){
                        print(e)
                        return(NA)
                    })
#    is.singular <- tryCatch(isSingular(mod),
#                            error = function(e){
#                                print(rxn)
#                                print(e)
#                                return(TRUE)
#                            })
    if(is.na(mod)){
        return(NA)
    }
#    if(!is.singular){
        stat <- data.table(as.data.frame(coef(summary(mod)))["RemissionR:Time_seq",],keep.rownames = TRUE)
        colnames(stat)[1] <- "coef"
        stat <- cbind(rxn = rxn,
                      stat)
        return(stat)
#    }else{
#        return(NA)
#    }
}

rxns <- unique(rxnExpr.melt[,rn])

threads <- detectCores() - 1
cl <- makeForkCluster(threads)
stats <- pblapply(cl = cl, FUN = createMods, rxns)
stopCluster(cl)
stats.coef <- data.table(do.call(rbind, stats[!is.na(stats)]))
stats.coef[,padj := p.adjust(`Pr(>|t|)`,method = "BH")]

write.csv(stats.coef, file = "../results/rxnExpr.LMM.Remission.TimeSeq_coefs.csv")

#save diagnostic plots for all the significant reactions
rxns <- stats.coef[padj < 0.05, rxn]
pbar <- txtProgressBar(min = 0 , max = length(rxns), style = 3)
pdf("../results/rxnExpr.LMM.Remission.TimeSeq_sigDiagPlots.pdf", width = 8, height= 5)
i <- 0 
for(rxn in rxns){
    mod <- lmer(data = rxnExpr.melt[rn == rxn,],Expression ~ Remission*Time_seq+(1|PatientID))
    print(plot_grid(~qqPlot(resid(mod), main = rxn), plot(mod)))
    i <- i+1
    setTxtProgressBar(pbar, i)
}
dev.off()

