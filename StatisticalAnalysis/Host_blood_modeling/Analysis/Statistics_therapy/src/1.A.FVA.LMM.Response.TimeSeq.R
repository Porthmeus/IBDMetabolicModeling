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

# load FVA results
FVA.min <- read.csv("../../../results/FVA/minFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
FVA.max <- read.csv("../../../results/FVA/maxFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
# convert the values to ranges and centers
FVA.range <- FVA.max-FVA.min
FVA.center <- (FVA.max+FVA.min)/2

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
meta <- meta[SeqID %in% colnames(FVA.min),]
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

# melt the FVA data to one data.table to make the associations
FVA.range.melt <- melt(
                       data.table(
                                  FVA.range[unique(tbl.cluster[,rep.rxn]),],
                                  keep.rownames = TRUE
                                  ),
                       id.vars = "rn",
                       variable.name = "SeqID",
                       value.name = "range"
                       )
FVA.center.melt <- melt(
                       data.table(
                                  FVA.center[unique(tbl.cluster[,rep.rxn]),],
                                  keep.rownames = TRUE
                                  ),
                       id.vars = "rn",
                       variable.name = "SeqID",
                       value.name = "center"
                       )

FVA <- merge(FVA.center.melt, FVA.range.melt)

# merge the meta data
FVA <- merge(FVA, meta, by="SeqID")
# remove controls
FVA <- FVA[Response != "C",]

# fit the models
threads <- detectCores()-1
cl <- makeForkCluster(threads)
registerDoParallel(cl)

# before construction of the models, we should define, whether center and ranges are highly correlated - if so there is no need to fit a model for both parameters.
rxns <- unique(FVA[,rn])

lowCor <- foreach(rxn = rxns,.combine = "c") %dopar% {
        cntr <- FVA[rn == rxn, center]
        rng <- FVA[rn == rxn, range]
        # look for high correlations in center and range
        corCenterRange <- sqrt(cor(cntr,rng)^2)
}
names(lowCor) <- unique(rxns)
stopCluster(cl)

forms <- data.table()
for(i in 1:length(lowCor)){
        corCenterRange <- lowCor[i]
        rxn <- names(lowCor[i])
        cntr <- FVA[rn == rxn, center]
        rng <- FVA[rn == rxn, range]
        # if one of the two has variance of 0, the cor.coef will be NA
        form <- "center ~ Response*Time_seq + (1|PatientID)"
        if (is.na(corCenterRange)){
            # check which var = 0 and use the other in lmm construction
            if(any(is.na(cntr)) | var(cntr)==0){
                form <- "range ~ Response*Time_seq + (1|PatientID)"
            } else if(any(is.na(rng)) | var(rng)==0){
                form <- "center ~ Response*Time_seq + (1|PatientID)"
            } else {
                stop(paste("Zero variance for center and range of FVA for",rxn))
            }
        # if the cor.coef < 0.6 use both features in formula, otherwise use just the center
        } else if(corCenterRange< 0.6){
            form <- unique(c(form,  "range ~ Response*Time_seq + (1|PatientID)"))
        }
        # cast it into a data frame
        forms <- rbind(forms, data.table(reac = names(lowCor)[i],
                              form = form))
}

# don't know why, but parallel execution does not work

createMods <- function(i){
    form <- as.formula(forms[i, form])
    rxn <- forms[i, reac]
    dat <- FVA[rn == rxn,]
    dat[,HB_Mayo_impu := (HB_Mayo_impu+1)]
    mod <- tryCatch(lmer(data = dat, form),
                    error = function(e){
                        print(e)
                        return(NA)
                    }
                    )
}

threads <- detectCores()-1
cl <- makeForkCluster(threads)
mods <- pblapply(1:nrow(forms), FUN = createMods, cl = cl)
stopCluster(cl)


# save the models object to come back later if needed
#saveRDS(mods, file = "../results/FVA.LMM.HBMayo.etiology_models.RDS")

# select the models which are not singular
mod_rxns <- forms[,reac]
#sel <- sapply(mods, isSingular)
#mod_rxns <- mod_rxns[!sel]
#mods <- mods[!sel]
stats <- lapply(mods, summary)
stats.coef <- lapply(stats, coef)
stats.coef <- lapply(stats.coef, as.data.frame)
stats.coef <- lapply(1:length(stats.coef), function(x) cbind(rxn = mod_rxns[x],coef = as.character(formula(mods[[x]])[[2]]), coef2 = rownames(stats.coef[[x]]), stats.coef[[x]]))
stats.coef <- data.table(do.call(rbind, stats.coef))
stats.coef <- stats.coef[coef2 == "ResponseR:Time_seq",]
stats.coef[,padj := p.adjust(`Pr(>|t|)`, method = "BH")] 
# save the coefficients table
write.csv(stats.coef, file = "../results/FVA.LMM.Response.TimeSeq_coefs.csv", row.names = FALSE)


##save diagnostic plots for all the significant reactions
#rxns <- stats.coef[padj <= 0.05, rxn]
#pbar <- txtProgressBar(min = 0 , max = length(rxns), style = 3)
#pdf("../results/FVA.LMM.Response.TimeSeq_sigDiagPlots.pdf", width = 8, height= 5)
#i <- 0 
#for(rxn in rxns){
#    print(plot_grid(~qqPlot(resid(mods[[rxn]]), main = rxn), plot(mods[[rxn]])))
#    i <- i+1
#    setTxtProgressBar(pbar, i)
#}
#dev.off()

# rerun the whole thing, without the time component
#lowCor.time <- lapply(lowCor, gsub, replacement = "Time_seq +(1+Time_seq", pattern = "(1", fixed = TRUE)
#
#mods.time <- foreach(rxn = unique(FVA[,rn])) %dopar% {
#    for(rxn in rxns){
#    # check if the range and center correlate in most of the cases
#    form <- as.formula(lowCor.time[rxn])
#    mod <-lmer(data = FVA[rn == rxn,],
#             form)
#}
#


