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
fluxes <- read.csv("../../../resources/future.eMed_outflux.csv", row.names = 1)

clinic <- fread("../../../resources/ClinicalDataMergeLongImputeCor.csv")
idconv <- fread("../../../resources/IdentifierConversion.csv")
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")

# create clusters for the fluxes
## remove near zero variance
nzv <- nearZeroVar(fluxes)
fluxes <- fluxes[,-nzv]

## scale the values
fluxes <- scale(fluxes)
## create a correlation matrix
fluxes.cor <- cor(fluxes)
# save the sign of correlation for later
corr.sign <- melt(data.table(sign(fluxes.cor),
                                        keep.rownames = TRUE),
                             id.var = "rn",
                             variable.name = "rep.rxn",
                             value.name = "sign.cor")
colnames(corr.sign)[1] <- "rxn" 
fluxes.cor <- 1-sqrt(fluxes.cor^2)

set.seed(395) # change in your code
cluster <- dbscan(fluxes.cor, # the matrix with the correlation coefficients converted to distance
                  method = "dist", # telling dbscan not to calculate the distances by itself
                  eps = 0.1, # the threshold how distant members are allowed to be to each other to count a cluster
                  MinPts = 3) # the number of how many reactions have to make up a cluster at least - 3 is the minimum

tbl.cluster <- data.table( cluster = cluster[["cluster"]],
                          eps = cluster[["eps"]],
                          MinPts = cluster[["MinPts"]],
                          isseed = cluster[["isseed"]],
                          rxn = rownames(fluxes.cor))

## get a representative reaction for each cluster
tbl.cluster[cluster == 0 ,rep.rxn := rxn]
tbl.cluster[cluster != 0, rep.rxn := names(which.min(rowSums(fluxes.cor[rxn,rxn]))), by = cluster]
# add back the sign of the correlation
tbl.cluster <- merge(tbl.cluster, corr.sign)
write.csv(file = "../results/outflux_cluster.csv", tbl.cluster)

# merge the data sets
## melt the fluxes 
fluxes.melt <- melt(data.table(fluxes[,unique(tbl.cluster[,rep.rxn])], keep.rownames = TRUE),
                     id.vars = "rn",
                     variable.name = "ex.rxn",
                     value.name = "flux")

## rename the samples to match the expected form
fluxes.melt[, SeqID := gsub("^X","",rn)]
fluxes.melt <- merge(fluxes.melt, idconv, by.x ="SeqID", by.y = "sample")
fluxes.melt[, PatientID_Time := paste(Cohort, Person.Numbercode, Time_seq, sep = "_")]
fluxes.merged <- merge(fluxes.melt, clinic, by = "PatientID_Time")


# fit the models
## register cores
threads <- detectCores() - 1
registerDoMC(threads)


rxns <- unique(fluxes.merged[,ex.rxn])

# fit the models for the fecal communities
mods <- foreach(rxn = rxns) %dopar%{
    dat <- fluxes.merged[ex.rxn == rxn ,]#& seqtype != "cDNA",]
    mod <- lmer(data = dat,
                formula = HB_Mayo_impu ~ source + seqtype + Age +  flux + (1|PatientID))
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
stats.coef <- data.table(do.call(rbind, stats.coef))[Coef == "flux",]
stats.coef[,padj := p.adjust(`Pr(>|t|)`, method = "BH")]
stats.coef[,cpd := gsub("^EX_(cpd[0-9]*)_e0","\\1", rxn)]
stats.coef <- merge(stats.coef, metabolites[,.(id, name)], by.x = "cpd", by.y= "id", all.x = TRUE)

write.csv(file = "../results/outflux.LMM.HBMayo.etiology_coefs.csv", stats.coef, row.names = FALSE)


# create the diagnostic plots

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


