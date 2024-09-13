# Porthmeus
# 05.04.22

# create LMM for the association between metabolites which are predicted to be produced by the microbiome and the patients in the FUTURE/eMed cohorts (IBD).


# load libraries
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

# load data sets
outfluxes <- read.csv("../../../resources/future.eMed.bacarena_outflux.csv", row.names = 1)
# remove internal fluxes
sel <- grep("fluxes.avail",colnames(outfluxes))
outfluxes <- outfluxes[,sel]

clinic <- fread("../../../resources/ClinicalDataMergeLongImputeCor.csv")
idconv <- fread("../../../resources/IdentifierConversion.csv")
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")
rxnsAnno <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_reactions.tsv")

# create clusters for the outfluxes
## remove near zero variance
nzv <- nearZeroVar(outfluxes)
outfluxes <- outfluxes[,-nzv]

## scale the values
outfluxes <- scale(outfluxes)
## create a correlation matrix
outfluxes.cor <- cor(outfluxes)
# save the sign of correlation for later
corr.sign <- melt(data.table(sign(outfluxes.cor),
                                        keep.rownames = TRUE),
                             id.var = "rn",
                             variable.name = "rep.rxn",
                             value.name = "sign.cor")
colnames(corr.sign)[1] <- "rxn" 
outfluxes.cor <- 1-sqrt(outfluxes.cor^2)

set.seed(528) # change in your code
cluster <- dbscan(outfluxes.cor, # the matrix with the correlation coefficients converted to distance
                  method = "dist", # telling dbscan not to calculate the distances by itself
                  eps = 0.1, # the threshold how distant members are allowed to be to each other to count a cluster
                  MinPts = 3) # the number of how many reactions have to make up a cluster at least - 3 is the minimum

tbl.cluster <- data.table( cluster = cluster[["cluster"]],
                          eps = cluster[["eps"]],
                          MinPts = cluster[["MinPts"]],
                          isseed = cluster[["isseed"]],
                          rxn = rownames(outfluxes.cor))

## get a representative reaction for each cluster
tbl.cluster[cluster == 0 ,rep.rxn := rxn]
tbl.cluster[cluster != 0, rep.rxn := names(which.min(rowSums(outfluxes.cor[rxn,rxn]))), by = cluster]
# add back the sign of the correlation
tbl.cluster <- merge(tbl.cluster, corr.sign)
write.csv(file = "../results/outflux_cluster.csv", tbl.cluster)

# merge the data sets
## melt the outfluxes 
outfluxes.melt <- melt(data.table(outfluxes[,unique(tbl.cluster[,rep.rxn])], keep.rownames = TRUE),
                     id.vars = "rn",
                     variable.name = "ex.rxn",
                     value.name = "outflux")

## rename the samples to match the expected form
outfluxes.melt[, SeqID := gsub("^X","",rn)]
outfluxes.melt <- merge(outfluxes.melt, idconv, by.x ="SeqID", by.y = "sample", all.x = TRUE)
outfluxes.melt[, PatientID_Time := paste(Cohort, Person.Numbercode, Time_seq, sep = "_")]
outfluxes.merged <- merge(outfluxes.melt, clinic, by = "PatientID_Time")


# fit the models
## register cores
threads <- detectCores() - 1
registerDoMC(threads)


rxns <- unique(outfluxes.merged[,ex.rxn])

# fit the models for the fecal communities
mods <- foreach(rxn = rxns) %dopar%{
    dat <- outfluxes.merged[ex.rxn == rxn,]# & seqtype != "cDNA",]
    mod <- lmer(data = dat,
                formula = HB_Mayo_impu ~ source + seqtype + Age + outflux + (1|PatientID))
}

names(mods) <- rxns
saveRDS(mods, file = "../results/outflux.LMM.HBMayo.etiology_mods.RDS")

# create the statistics and sav them
mc <- makeCluster(mc <- threads)
stats <- parLapply(mc, mods, summary)
stopCluster(mc)

stats.coef <- lapply(1:length(stats), function(x) {
                         cbind("Coef" = rownames(coef(stats[[x]])), "rxn" = names(stats)[x], as.data.frame(coef(stats[[x]])))}
)
stats.coef <- data.table(do.call(rbind, stats.coef))[Coef == "outflux",]
stats.coef[,padj := p.adjust(`Pr(>|t|)`, method = "BH")]
stats.coef[,cpd := gsub("^.*_(cpd[0-9]*)","\\1", rxn)]
stats.coef <- merge(stats.coef, metabolites[,.(id, name)], by.x = "cpd", by.y= "id", all.x = TRUE)
#stats.coef[,rxn_id := gsub("_.0$","",rxn)]
#stats.coef <- merge(stats.coef, rxnsAnno[,.(id, name, definition)], by.x = "rxn_id", by.y ="id", all.x = TRUE, suffix = c("","_rxn"))

write.csv(file = "../results/outflux.LMM.HBMayo.etiology_coefs.csv", stats.coef, row.names = FALSE)


rxns <- stats.coef[padj <= 0.05, rxn]
if(length(rxns) > 0){
    pbar <- txtProgressBar(min = 0 , max = length(rxns), style = 3)
    pdf("../results/outflux.LMM.HBMayo.etiology_sigDiagPlots.pdf", width = 8, height= 5)
    i <- 0 
    for(rxn in rxns){
        print(plot_grid(~qqPlot(resid(mods[[rxn]]), main = rxn), plot(mods[[rxn]])))
        i <- i+1
        setTxtProgressBar(pbar, i)
    }
    dev.off()
}
