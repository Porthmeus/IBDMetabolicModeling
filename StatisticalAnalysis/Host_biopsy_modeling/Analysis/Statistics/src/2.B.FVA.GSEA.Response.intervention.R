# Porthmeus
# 28.02.22

# make a GSEA for FVA results

# load libraries
require(data.table)
require(clusterProfiler)
require(lme4)
require(lmerTest)
require(caret)
require(foreach)
require(parallel)
require(doMC)

# load FVA results
FVA.min <- read.csv("../../../results/FVA/minFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)

subsystems <- fread("../dat/subsystems.csv")
cluster <- fread("../results/FVA_DBSCAN_Cluster.csv")
stats <- fread("../results/FVA.GLMM.Response.intervention_coefs.csv")
# sum the estimates for the reactions, thus there is only one for each of the reactions in case both center and range are calculated
stats <- stats[,.(Estimate = sum(Estimate)), by = rxn]

# create a vector for the GSEA from the coefficients of the LMM
## first expand the cluster
stats.exp <- merge(stats, cluster, by.x = "rxn", by.y = "rep.rxn")

## get the estimates
estimates <- stats.exp[,Estimate]
names(estimates) <- stats.exp[,rxn]

## add 0 for the reactions which have been removed due to low variance
## don't do it anymore, because it will raise problems downstream
#rmvd.rxns <- setdiff(rownames(FVA.min), stats.exp[,rxn])
#estimates.rmvd <- sapply(rmvd.rxns, function(x) 0 )

## merge the two vectors and sort them decreasingly
#estimates.all <- sort(c(estimates, estimates.rmvd), decreasing = TRUE)
estimates.all <- sort(estimates, decreasing =TRUE)

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
                 maxGSSize = nrow(subsystems),
                 verbose = FALSE
                 )

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
stability <- stability[stability.i]
stab <- as.data.frame(table(unlist(stability))/length(stability))
subs <- stab[stab[,2]>0.8, 1]

# find a representative gsea analysis set to save, while updating the stability anlalysis.
gsea <- NULL
while(is.null(gsea)){
    tryCatch({gsea <- GSEA(estimates.all,
                 TERM2GENE = subsystems,
                 minGSSize = 3,
                 maxGSSize = nrow(subsystems),
                 verbose = FALSE
                 )

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
                 maxGSSize = nrow(subsystems),
                 verbose = FALSE
                 )
            },
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
write.csv(stab, file = "../results/FVA.GSEA.Response.intervention_stabilityAnalysis.csv")
write.csv(as.data.frame(estimates.all), file = "../results/FVA.GSEA.Response.intervention_rankingVector.csv")
saveRDS(gsea, file = "../results/FVA.GSEA.Response.intervention_GSEAResult.RDS")
