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

# load cluster table
tbl.cluster <- fread("../results/FVA_DBSCAN_Cluster.csv")

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
FVA <- FVA[Response != "C" & Time_seq %in% c(0,14),]
FVA.center.dcast <- dcast(data = FVA, PatientID+rn~paste0("d",Time_seq ), value.var = "center")
FVA.center.dcast[,center.diff := d14 - d0]
FVA.range.dcast <- dcast(data = FVA, PatientID+rn~paste0("d",Time_seq), value.var = "range")
FVA.range.dcast[,range.diff := d14 - d0]
dat <- merge(FVA.range.dcast[,.(PatientID,rn,range.diff)], FVA.center.dcast[,.(PatientID,rn,center.diff)], by = c("PatientID","rn"))
dat <- merge(dat, unique(meta[,.(PatientID, Response)]), by = "PatientID")


# fit the models
threads <- detectCores()-1
cl <- makeForkCluster(threads)
registerDoParallel(cl)
# before construction of the models, we should define, whether center and ranges are highly correlated - if so there is no need to fit a model for both parameters.
rxns <- unique(dat[,rn])
lowCor <- foreach(rxn = rxns,.combine = "c") %dopar% {
        cntr <- dat[rn == rxn, center.diff]
        rng <- dat[rn == rxn, range.diff]
        # look for high correlations in center and range
        corCenterRange <- sqrt(cor(cntr,rng, use= "na.or.complete")^2)
}
names(lowCor) <- unique(rxns)
stopCluster(cl)

forms <- data.table()
for(i in 1:length(lowCor)){
        corCenterRange <- lowCor[i]
        rxn <- names(lowCor[i])
        cntr <- dat[rn == rxn, center.diff]
        rng <- dat[rn == rxn, range.diff]
        # if one of the two has variance of 0, the cor.coef will be NA
        form <- "factor(Response) ~ center.diff"
        if (is.na(corCenterRange)){
            # check which var = 0 and use the other in lmm construction
            if(any(is.na(cntr)) | var(cntr)==0){
                form <- "factor(Response) ~ range.diff"
            } else if(any(is.na(rng)) | var(rng)==0){
                form <- "factor(Response) ~ center.diff"
            } else {
                stop(paste("Zero variance for center and range of FVA for",rxn))
            }
        # if the cor.coef < 0.6 use both features in formula, otherwise use just the center
        } else if(corCenterRange< 0.6){
            form <- unique(c(form,  "factor(Response) ~ center.diff + range.diff"))
        }
        # cast it into a data frame
        forms <- rbind(forms, data.table(reac = names(lowCor)[i],
                              form = form))
}

# don't know why, but parallel execution does not work

createMods <- function(i){
    form <- as.formula(forms[i, form])
    rxn <- forms[i, reac]
    dat <- dat[rn == rxn,]
    mod <- tryCatch(glm(data = dat, form, family = "binomial"),
                    error = function(e){
                        print(e)
                        return(NA)
                    }
                    )
    if(length(mod)>1&all(!is.na(coef(mod)))){
        form_split <- names(coef(mod))[-(1)]
        if(length(form_split) >1){
            dat <- as.data.frame(coef(summary(mod))[form_split,])
        } else {
            dat <- as.data.frame(t(coef(summary(mod))[form_split,]))
        }
        stat <- cbind(
                          rxn = rxn, # add the reaction name to the table
                          coef = form_split,
                          dat
                          )
        return(stat)
    } else { return(NA)}
}

rxns <- dat[,unique(rn)]
threads <- detectCores()-1
cl <- makeForkCluster(threads)
stats <- pblapply(1:nrow(forms), FUN = createMods, cl = cl)
stopCluster(cl)

sel <- is.na(stats)
stats.all <- stats[!sel]
stats.all <- data.table(do.call(rbind, stats.all))
stats.all[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
stats.all[padj < 0.05,]
write.csv(stats.all, file = "../results/FVA.LMM.Response.14d0dDiff_coefs.csv", row.names = FALSE)



#save diagnostic plots for all the significant reactions
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


