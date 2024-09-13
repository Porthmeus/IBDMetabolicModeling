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
require(fpc)
require(car)
require(cowplot)
require(DHARMa) 
require(pbapply)

# load data sets
outfluxes <- read.csv("../../../resources/future.eMed_outflux.csv", row.names = 1)

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
outfluxes.merged <- outfluxes.merged[Responder != "C",]
outfluxes.merged <- outfluxes.merged[,Response:= factor(Responder, levels =c("NR","R"))]
outfluxes.merged <- outfluxes.merged[,Remission:= factor(Remision, levels =c("NR","R"))]
outfluxes.merged <- outfluxes.merged[,Time_seq:= Time_seq.x]

getModelStat <- function(rxn, what){
    dat <- outfluxes.merged[ex.rxn == rxn,]#& seqtype != "cDNA",]
    form <- paste0("outflux ~ source + seqtype + Age + ", what, " *Time_seq + (1|PatientID)")
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

rxns <- outfluxes.merged[,unique(ex.rxn)]
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
fwrite(file = paste0("../results/outflux.", outcome,".TimeSeq_coefs.csv"), stats.coef, row.names = FALSE)
}
