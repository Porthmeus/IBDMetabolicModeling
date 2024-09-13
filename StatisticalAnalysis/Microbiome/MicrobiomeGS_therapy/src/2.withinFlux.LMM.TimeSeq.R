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

# remove control patients and reformat Response
fluxes.merged <- fluxes.merged[Responder != "C",]
fluxes.merged <- fluxes.merged[,Response:= factor(Responder, levels =c("NR","R"))]
fluxes.merged <- fluxes.merged[,Remission:= factor(Remision, levels =c("NR","R"))]
fluxes.merged <- fluxes.merged[,Time_seq:= Time_seq.x]
rxns <- fluxes.merged[,unique(ex.rxn)]

getModelStat <- function(rxn, what){
    dat <- fluxes.merged[ex.rxn == rxn,]#& seqtype != "cDNA",]
    form <- paste0("flux ~ source + seqtype + Age + ", what, " *Time_seq + (1|PatientID)")
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
