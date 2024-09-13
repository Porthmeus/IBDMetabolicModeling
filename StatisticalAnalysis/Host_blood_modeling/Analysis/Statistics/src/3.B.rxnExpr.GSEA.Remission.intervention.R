# Porthmeus
# 09.03.22

# load libraries
require(data.table)
require(clusterProfiler)
require(lme4)
require(lmerTest)
require(caret)
require(foreach)
require(parallel)
require(doMC)

# load data
rxnExpr <- read.csv("../../../results/RxnExpression/rxnExpr.GL10-L50-GU90.colormore3D.csv", row.names = 1)
subsystems <- fread("../dat/subsystems.csv")
cluster <- fread("../results/rxnExpr_DBSCAN_Cluster.csv")
stats <- fread("../results/rxnExpr.GLMM.Remission.intervention_coefs.csv")

# create a vector for the GSEA from the coefficients of the LMM
## first expand the cluster
stats.exp <- merge(stats, cluster, by.x = "rxn", by.y = "rep.rxn", suffix = c("","_cluster"))

## get the estimates
estimates <- stats.exp[,Estimate]
names(estimates) <- stats.exp[,rxn_cluster]

## add 0 for the reactions which have been removed due to low variance
## don't do this anymore, causes trouble downstream
#rmvd.rxns <- setdiff(rownames(rxnExpr), stats.exp[,rxn_cluster])
#estimates.rmvd <- sapply(rmvd.rxns, function(x) 0 )

## merge the two vectors and sort them decreasingly
#estimates.all <- sort(c(estimates, estimates.rmvd), decreasing = TRUE)
estimates.all <- sort(estimates, decreasing = TRUE)

# due to the many ties I need to make some stability analysis
# I simply repeat the GSEA 100 times and count how often a subsystem is enriched
threads <- detectCores() -1
registerDoMC(threads)
#
n <- 100
stability <- foreach(i = 1:n) %dopar% {
    tryCatch({gsea <- GSEA(estimates.all,
                 TERM2GENE = subsystems,
                 minGSSize = 3,
                 maxGSSize = nrow(subsystems))
                print(paste0("Job: ",i,"/",n))
                data.frame(gsea)[,"ID"]},
             error = function(e){ print(as.character(e))
                                  print(paste("Will ignore error in job",i))
                                  NULL
             })
}

# get those subsystems which appear at least in 80% of the analysis
stability.i <- !sapply(stability, is.null)
n <- sum(stability.i)
if(n > 0){
    stability <- stability[stability.i]
    stab <- as.data.frame(table(unlist(stability))/length(stability))
    subs <- stab[stab[,2]>0.8, 1]
} else {
#    warning("None of the GSEAs delivered a result - usually a problem in the structure of the vector for the analysis\n Will fall back to correlations")
    stop("\nNone of the GSEAs delivered a result - usually a problem in the structure of the vector for the analysis.\n")
       
## load metadata
#clinic <- fread("../../../resources/ClinicalData.csv")
#meta <- fread("../../../resources/META_emed_future.csv")
#blcklst <- fread("../../../resources/sampleBlacklist.csv")
#subsystems <- fread("../dat/subsystems.csv")
#
## merge meta and clincal data
#meta[,PatientID_Time := paste(PatientID, Time_seq, sep ="_")]
#meta <- meta[Time_seq != -1,]
#meta <- merge(meta,clinic, by="PatientID_Time", all.x=TRUE, suffix = c(".meta",""))
#
## remove degraded samples
#blckSmpls <- blcklst[relevance < 2, Sample]
#meta <- meta[!(SeqID %in% blckSmpls),]
#meta <- meta[SeqID %in% colnames(rxnExpr),]
#rxnExpr <- rxnExpr[,!(colnames(rxnExpr) %in% blckSmpls)]
#cor.x <- as.integer(as.factor(meta[Remission != "C", Remission]))
#rxnExpr.y <- rxnExpr[,meta[Remission != "C",SeqID]]
#estimate.all <- sort(apply(rxnExpr.y, 1, function(x) {
#                        tryCatch(cor(x, cor.x),
#                                 warning = function(w) 0,
#                                 error = function(e) 0)}),
#                     decreasing = TRUE)
#
#threads <- detectCores() -1
#registerDoMC(threads)
##
#n <- 10
#stability <- foreach(i = 1:n) %dopar% {
#    tryCatch({gsea <- GSEA(estimates.all,
#                 TERM2GENE = subsystems,
#                 minGSSize = 3,
#                 maxGSSize = nrow(subsystems))
#                print(paste0("Job: ",i,"/",n))
#                data.frame(gsea)[,"ID"]},
#             error = function(e){ print(as.character(e))
#                                  print(paste("Will ignore error in job",i))
#                                  NULL
#             })
#}

}
# get those subsystems which appear at least in 80% of the analysis
stability.i <- !sapply(stability, is.null)
n <- sum(stability.i)
if(n > 0){
    stability <- stability[stability.i]
    stab <- as.data.frame(table(unlist(stability))/length(stability))
    subs <- stab[stab[,2]>0.8, 1]

}

# find a representative gsea analysis set to save, while updating the stability anlalysis.
gsea <- NULL
while(is.null(gsea)){
    tryCatch({gsea <- GSEA(estimates.all,
                 TERM2GENE = subsystems,
                 minGSSize = 3,
                 maxGSSize = nrow(subsystems))
                 n <- n+1
                 stability[[n]] <- data.frame(gsea)[,"ID"]
                 stab <- as.data.frame(table(unlist(stability))/n)
                 subs <- stab[stab[,2]>0.8, 1]
                 },
             error = function(e){ print(as.character(e))
                                  NULL
             })
}

while(is.null(gsea) | !all(sort(data.frame(gsea)[,"ID"]) == sort(subs))){
    tryCatch({gsea <- GSEA(estimates.all,
                 TERM2GENE = subsystems,
                 minGSSize = 3,
                 maxGSSize = nrow(subsystems))},
             error = function(e){ print(as.character(e))
                                  print(paste("Will ignore error in job",i))
                                  NULL
             })
    if(!is.null(gsea)){
        n <- length(stability) +1
        stability[[n]] <- data.frame(gsea)[,"ID"]
        stab <- as.data.frame(table(unlist(stability))/n)
        subs <- stab[stab[,2]>0.8, 1]
        # in case there is no valid solution after 50 attempts, take the current result and remove all extra subsystems
        if(n > 150){
            if(all(subs %in% as.data.frame(gsea)[,"ID"])){
                gsea@result <- as.data.frame(gsea)[as.data.frame(gsea)[,"ID"] %in% subs,]
                break
            }
        }
    }
}

# reformat the stability analysis table
colnames(stab) <- c("subsystem","FractionObtained")
stab <- cbind(stab, reps = n)

# save the results
write.csv(stab, file = "../results/rxnExpr.GSEA.Remission.intervention_stabilityAnalysis.csv")
write.csv(as.data.frame(estimates.all), file = "../results/rxnExpr.GSEA.Remission.intervention_rankingVector.csv")
saveRDS(gsea, file = "../results/rxnExpr.GSEA.Remission.intervention_GSEAResult.RDS")

