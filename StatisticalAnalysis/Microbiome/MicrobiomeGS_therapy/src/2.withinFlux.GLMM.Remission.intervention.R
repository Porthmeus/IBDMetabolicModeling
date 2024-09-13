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
require(pbapply)
require(fpc)
require(car)
require(cowplot)
require(DHARMa)

# load data sets
fluxes <- read.csv("../../../resources/future.eMed_withinOflux.csv", row.names = 1)

clinic <- fread("../../../resources/ClinicalDataMergeLongImputeCor.csv")
idconv <- fread("../../../resources/IdentifierConversion.csv")
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")
pwys <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/meta_rea_pwy-gapseq.tbl")
pwys <- pwys[,.(rxn_id = unlist(strsplit(gapseq, split = " "))), by = .(pathway.name, pathway)]
pwys <-pwys[!duplicated(pwys),]
pwys <- pwys[,.(name = paste(pathway.name, collapse = ";"), pathway = paste(pathway, collapse = ";")), by = rxn_id]
rxnsAnno <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_reactions.tsv")

tbl.cluster <- fread("../results/withinFlux_cluster.csv")

# remove unnecessary fluxes
fluxes <- fluxes[,tbl.cluster[,unique(rep.rxn)]]

## scale the values
fluxes <- scale(fluxes)

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

# remove control patients and reformat Remission
fluxes.merged <- fluxes.merged[Remision != "C" & Time_seq.x != max(Time_seq.x),]
fluxes.merged <- fluxes.merged[,Remission:= factor(Remision, levels =c("NR","R"))]

# fit the models
## register cores
threads <- detectCores() - 1
cl <- makeForkCluster(threads)
registerDoParallel(cl)


rxns <- unique(fluxes.merged[,ex.rxn])

# fit the models for the fecal communities
mods <- foreach(rxn = rxns) %dopar%{
    dat <- fluxes.merged[ex.rxn == rxn & Remission != "C" ,]#& seqtype != "cDNA",]
    mod <- tryCatch({glmer(data = dat,
                formula = Remission ~ source + seqtype + Age + flux + (1|PatientID),
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
#saveRDS(mods, file = "../results/withinFlux.GLMM.Remission.intervention_mods.RDS")

# create the statistics and sav them
mc <- makeCluster(mc <- threads)
stats <- parLapply(mc, mods, summary)
stopCluster(mc)
stats.coef <- lapply(1:length(stats), function(x) {
                         cbind("coef" = rownames(coef(stats[[x]])), "rxn" = names(stats)[x], as.data.frame(coef(stats[[x]])))}
)
stats.coef <- data.table(do.call(rbind, stats.coef))[coef == "flux",]
stats.coef[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
stats.coef[,cpd := gsub(".*_(cpd[0-9]*)_\\w\\d","\\1", rxn)]
stats.coef <- merge(stats.coef, metabolites[,.(id, name)], by.x = "cpd", by.y= "id", all.x = TRUE)
#stats.coef[,rxn_id := gsub("_.0$","",rxn)]
#stats.coef <- merge(stats.coef, rxnsAnno[,.(id, name, definition)], by.x = "rxn_id", by.y ="id", all.x = TRUE, suffix = c("","_rxn"))

write.csv(file = "../results/withinFlux.Remission_coefs.csv", stats.coef, row.names = FALSE)


# do the simple linear models as well
dats <- list()
dats[["d0"]] <- fluxes.merged[Remission != "C" & Time_seq.x == 0, .(Remission, source, seqtype, Age, value = flux, PatientID, ex.rxn)]
dats[["d14"]] <- fluxes.merged[Remission != "C" & Time_seq.x == 14, .(Remission, source, seqtype, Age, value = flux, PatientID, ex.rxn)]
# calculate difference between time points
dat.temp <- fluxes.merged[Remission != "C" & Time_seq.x %in% c(0, 14), .(Remission, source, seqtype, Age, value = flux, PatientID, ex.rxn, Time_seq = Time_seq.x)]
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
    stats <- data.table(do.call(rbind, stats))
    stats[,c("outcome", "coef","padj") := list("Remission",dat.n, p.adjust(`Pr(>|z|)`, method = "BH"))]
    fwrite(stats, file = paste0("../results/withinFlux.Remission.",dat.n,"_coefs.csv"))
}
