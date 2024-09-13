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
fluxes <- read.csv("../../../resources/future.eMed_flux.csv", row.names = 1)

clinic <- fread("../../../resources/ClinicalDataMergeLongImputeCor.csv")
idconv <- fread("../../../resources/IdentifierConversion.csv")
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")

tbl.cluster <- fread("../results/flux_cluster.csv")

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
getGLMMstats <- function(rxn){
    dat <- fluxes.merged[ex.rxn == rxn,]#& seqtype != "cDNA",]
    mod <- tryCatch({glmer(data = dat,
                formula = Remission ~ source + seqtype + Age + flux + (1|PatientID),
                family = "binomial")},
                    error = function(e) {
                        print(paste0("Error to build the model for ",rxn))
                        print(e)
                        return(NA)
                    }
                    )
    if(!is.na(mod)){
        stats <- cbind(rxn = rxn, as.data.frame(coef(summary(mod)))["flux",])
        return(stats)
    } else {
        return(NA)
    }
}

rxns <- unique(fluxes.merged[,ex.rxn])
threads <- detectCores() - 1
cl<- makeForkCluster(threads)
registerDoParallel(cl)
mods <- pblapply(cl = cl, rxns, getGLMMstats)
stopCluster(cl)

stats.coef <- data.table(do.call(rbind, mods))
stats.coef[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
stats.coef[,id := gsub("_.0$","", rxn)]
write.csv(file = "../results/flux.Remission_coefs.csv", stats.coef, row.names = FALSE)


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
    fwrite(stats, file = paste0("../results/flux.Remission.",dat.n,"_coefs.csv"))
}
