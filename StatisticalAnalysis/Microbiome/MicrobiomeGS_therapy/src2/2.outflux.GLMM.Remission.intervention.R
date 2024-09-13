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
require(doParallel)
require(doMC)
require(fpc)
require(car)
require(cowplot)
require(DHARMa) 
require(pbapply)

# load data sets
outfluxes <- read.csv("../../../resources/future.eMed.bacarena_outflux.csv", row.names = 1)
# remove internal fluxes
sel <- grep("fluxes.avail",colnames(outfluxes))
outfluxes <- outfluxes[,sel]

clinic <- fread("../../../resources/ClinicalDataMergeLongImputeCor.csv")
idconv <- fread("../../../resources/IdentifierConversion.csv")
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")
rxnsAnno <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_reactions.tsv")

tbl.cluster <- fread("../results/outflux_cluster.csv")

# remove unnecessary outfluxes
outfluxes <- outfluxes[,tbl.cluster[,unique(rep.rxn)]]

## scale the values
outfluxes <- scale(outfluxes)

# merge the data sets
## melt the outfluxes 
outfluxes.melt <- melt(data.table(outfluxes[,unique(tbl.cluster[,rep.rxn])], keep.rownames = TRUE),
                     id.vars = "rn",
                     variable.name = "ex.rxn",
                     value.name = "outflux")

## rename the samples to match the expected form
outfluxes.melt[, SeqID := gsub("^X","",rn)]
outfluxes.melt <- merge(outfluxes.melt, idconv, by.x ="SeqID", by.y = "sample")
outfluxes.melt[, PatientID_Time := paste(Cohort, Person.Numbercode, Time_seq, sep = "_")]
outfluxes.merged <- merge(outfluxes.melt, clinic, by = "PatientID_Time")

# remove control patients and reformat Response
outfluxes.merged <- outfluxes.merged[Remision != "C" & Time_seq.x != max(Time_seq.x),]
outfluxes.merged <- outfluxes.merged[,Remission:= factor(Remision, levels =c("NR","R"))]

# fit the models
## register cores
threads <- detectCores() - 1
cl <- makeForkCluster(threads)
registerDoParallel(cl)


rxns <- unique(outfluxes.merged[,ex.rxn])

# fit the models for the fecal communities
mods <- foreach(rxn = rxns) %dopar%{
    dat <- outfluxes.merged[ex.rxn == rxn,]#& seqtype != "cDNA",]
    mod <- tryCatch({glmer(data = dat,
                formula = Remission ~  source + seqtype + Age + outflux + (1|PatientID),
                family = "binomial")},
                    error = function(e) {
                        print(paste0("Error to build the model for ",rxn))
                        print(e)
                        return(NA)
                    }
                    )
}
stopCluster(cl)

names(mods) <- rxns
sel <- sapply(mods, is.na)
mods <- mods[!sel]
#saveRDS(mods, file = "../results/outflux.GLMM.Remission.intervention_mods.RDS")

# create the statistics and save them
mc <- makeForkCluster(threads)
stats <- pblapply(cl = mc, mods, summary)
stopCluster(mc)
stats.coef <- lapply(1:length(stats), function(x) {
                         cbind("Coef" = rownames(coef(stats[[x]])), "rxn" = names(stats)[x], as.data.frame(coef(stats[[x]])))}
)
stats.coef <- data.table(do.call(rbind, stats.coef))[Coef == "outflux",]
stats.coef[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
stats.coef[,cpd := gsub(".*_(cpd[0-9]*)","\\1", rxn)]
stats.coef <- merge(stats.coef, metabolites[,.(id, name)], by.x = "cpd", by.y= "id", all.x = TRUE)
#stats.coef[,rxn_id := gsub("_.0$","",rxn)]
#stats.coef <- merge(stats.coef, rxnsAnno[,.(id, name, definition)], by.x = "rxn_id", by.y ="id", all.x = TRUE, suffix = c("","_rxn"))

write.csv(file = "../results/outflux.Remission_coefs.csv", stats.coef, row.names = FALSE)

# do the simple linear models as well
dats <- list()
dats[["d0"]] <- outfluxes.merged[Remission != "C" & Time_seq.x == 0, .(Remission, source, seqtype, Age, value = outflux, PatientID, ex.rxn)]
dats[["d14"]] <- outfluxes.merged[Remission != "C" & Time_seq.x == 14, .(Remission, source, seqtype, Age, value = outflux, PatientID, ex.rxn)]
# calculate difference between time points
dat.temp <- outfluxes.merged[Remission != "C" & Time_seq.x %in% c(0, 14), .(Remission, source, seqtype, Age, value = outflux, PatientID, ex.rxn, Time_seq = Time_seq.x)]
dat.mat <- dcast(data=dat.temp, Age+Remission+source+seqtype+PatientID+ex.rxn~paste0("d",Time_seq), value.var = "value")
dat.mat[, value := d14-d0]
dat.mat <- dat.mat[!is.na(value),]
dats[["d14d0Diff"]] <- dat.mat

getModelStats <- function(rxn){
    dat.temp <- dat[ex.rxn == rxn,]
    mod <- glm(data = dat.temp,
               formula = factor(Remission) ~ source+seqtype+Age+value,
               family = "binomial")
    stats <- cbind(rxn = rxn, as.data.frame(coef(summary(mod)))["value",])
    return(stats)
    }


for(dat.n in names(dats)){
    dat <- dats[[dat.n]]
    cl <- makeForkCluster(threads)
    stats <- pblapply(dat[,unique(ex.rxn)], getModelStats, cl=cl)
    stopCluster(cl)
    stats <- data.table(do.call(rbind, stats))
    stats[,c("outcome", "coef","padj") := list("Remission",dat.n, p.adjust(`Pr(>|z|)`, method = "BH"))]
    fwrite(stats, file = paste0("../results/outflux.Remission.",dat.n,"_coefs.csv"))
}




## create the diagnostic plots
#
#rxns <- stats.coef[padj <= 0.05, rxn]
#if(length(rxns) > 0){
#pbar <- txtProgressBar(min = 0 , max = length(rxns), style = 3)
#pdf("../results/outflux.GLMM.Response.intervention_sigDiagPlots.pdf", width = 8, height= 5)
#i <- 0 
#for(rxn in rxns){
#    sim <- simulateResiduals(mods[[rxn]])
#    plot(sim)
#    mtext(rxn)
#    i <- i+1
#    setTxtProgressBar(pbar, i)
#}
#dev.off()
#}
#
#
