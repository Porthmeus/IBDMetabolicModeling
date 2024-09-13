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
require(pbapply)
require(DHARMa)


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

rxnCount2 <- rxnCount[unique(tbl.cluster[, rep.rxn]),]
rxnCount2.melt <- melt(data.table(rxnCount2, keep.rownames=TRUE), id.vars = "rn")

# merge with meta
dat <- merge(rxnCount2.melt, meta, by.x = "variable",by.y= "SeqID")
dat <- dat[Response != "C",]


# create a list of lmer models
# for parallel computation register threads

createMods <-function(rxn){
    #foreach(rxn = rownames(rxnCount2)) %dopar% {
    dat.mod <- dat[rn == rxn,]
    mod <- tryCatch(glmer(data = dat.mod, value ~ Response*Time_seq+(1|PatientID), family = "binomial"),
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
        stat <- data.table(as.data.frame(coef(summary(mod)))["ResponseR:Time_seq",],keep.rownames = TRUE)
        colnames(stat)[1] <- "coef"
        stat <- cbind(rxn = rxn,
                      stat)
        return(stat)
#    }else{
#        return(NA)
#    }
}

rxns <- dat[,unique(rn)]

threads <- detectCores()-1
cl <- makeForkCluster(threads)
#clusterExport(cl, "dat")
#registerDoParallel(threads)
mods <- pblapply(rxns, FUN = createMods, cl = cl)
stopCluster(cl)
names(mods) <- rxns

# remove singular and invalid models
sel <- !is.na(mods)
mods <- mods[sel]
stats <- do.call(rbind, mods)
stats[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]

write.csv(stats, file = "../results/PA.LMM.Response.TimeSeq_coefs.csv", row.names = FALSE)


## create diagnostic plots for the significant reactions
#rxns.sig <- stats[padj < 0.05, rxn]
#pdf("../results/PA.LMM.Response_sigModsDiagnostics.pdf")
#for(rxn in rxns.sig){
#    dat.mod <- dat[rn == rxn,]
#    mod <-glmer(data = dat.mod, value ~ Response*Time_seq+(1|PatientID), family = "binomial")
#    sim <- simulateResiduals(mod)
#    plot(sim)
#    mtext(rxn)
#}
#dev.off()
