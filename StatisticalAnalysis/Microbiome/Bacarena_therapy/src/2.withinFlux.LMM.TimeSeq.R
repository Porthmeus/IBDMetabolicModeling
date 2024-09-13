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
withinFluxes <- read.csv("../../../resources/future.eMed.bacarena_outflux.csv", row.names = 1)
# remove internal fluxes
sel <- grep("fluxes.inter",colnames(withinFluxes))
withinFluxes <- withinFluxes[,sel]

clinic <- fread("../../../resources/ClinicalDataMergeLongImputeCor.csv")
idconv <- fread("../../../resources/IdentifierConversion.csv")
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")
rxnsAnno <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_reactions.tsv")

tbl.cluster <- fread("../results/withinFlux_cluster.csv")

# remove unnecessary withinFluxes
withinFluxes <- withinFluxes[,tbl.cluster[,unique(rep.rxn)]]

## scale the values
withinFluxes <- scale(withinFluxes)

# merge the data sets
## melt the withinFluxes 
withinFluxes.melt <- melt(data.table(withinFluxes[,unique(tbl.cluster[,rep.rxn])], keep.rownames = TRUE),
                     id.vars = "rn",
                     variable.name = "ex.rxn",
                     value.name = "withinFlux")

## rename the samples to match the expected form
withinFluxes.melt[, SeqID := gsub("^X","",rn)]
withinFluxes.melt <- merge(withinFluxes.melt, idconv, by.x ="SeqID", by.y = "sample")
withinFluxes.melt[, PatientID_Time := paste(Cohort, Person.Numbercode, Time_seq, sep = "_")]
withinFluxes.merged <- merge(withinFluxes.melt, clinic, by = "PatientID_Time")

# remove control patients and reformat Response
withinFluxes.merged <- withinFluxes.merged[Responder != "C",]
withinFluxes.merged <- withinFluxes.merged[,Response:= factor(Responder, levels =c("NR","R"))]
withinFluxes.merged <- withinFluxes.merged[,Remission:= factor(Remision, levels =c("NR","R"))]
withinFluxes.merged <- withinFluxes.merged[,Time_seq:= Time_seq.x]
rxns <- withinFluxes.merged[,unique(ex.rxn)]

getModelStat <- function(rxn, what){
    dat <- withinFluxes.merged[ex.rxn == rxn,]#& seqtype != "cDNA",]
    form <- paste0("withinFlux ~ source + seqtype + Age + ", what, " *Time_seq + (1|PatientID)")
    mod <- tryCatch({lmer(data = dat,
                formula = as.formula(form))},
                    error = function(e) {
                        print(paste0("Error to build the model for ",rxn))
                        print(e)
                        return(NA)
                    }
                    )
    if(suppressWarnings(!is.na(mod))){
        stats <- as.data.frame(coef(summary(mod)))
        sel <- grep(":",rownames(stats))
        stats <- data.table(cbind(rxn= rxn, stats[sel,]))
        return(stats)
    } else {
        return(NA)
    }
}

for(outcome in c("Response", "Remission")){

# fit the models
## register cores
threads <- detectCores() - 1
cl <- makeForkCluster(threads)
registerDoParallel(cl)
stats <- pblapply(cl = cl, rxns, getModelStat, what = outcome)
stopCluster(cl)

stats.coef <- do.call(rbind, stats)
stats.coef[,padj := p.adjust(`Pr(>|t|)`, method = "BH")]
fwrite(file = paste0("../results/withinFlux.", outcome,".TimeSeq_coefs.csv"), stats.coef, row.names = FALSE)
}
