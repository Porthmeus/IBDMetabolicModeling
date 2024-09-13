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
require(DHARMa)

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
outfluxes.merged <- outfluxes.merged[Responder != "C",]
outfluxes.merged <- outfluxes.merged[,Response:= factor(Responder, levels =c("NR","R"))]

# fit the models
## register cores
threads <- detectCores() - 1
registerDoMC(threads)


rxns <- unique(outfluxes.merged[,ex.rxn])

# fit the models for the fecal communities
mods <- foreach(rxn = rxns) %dopar%{
    dat <- outfluxes.merged[ex.rxn == rxn & Response != "C" ,]#& seqtype != "cDNA",]
    mod <- tryCatch({glmer(data = dat,
                formula = Response ~HB_Mayo_impu + source + seqtype + Age + outflux + (1|PatientID),
                family = "binomial")},
                    error = function(e) {
                        print(paste0("Error to build the model for ",rxn))
                        print(e)
                        return(NA)
                    }
                    )
}

names(mods) <- rxns
sel <- sapply(mods, is.na)
mods <- mods[!sel]
saveRDS(mods, file = "../results/outflux.GLMM.Response.intervention_mods.RDS")

# create the statistics and sav them
mc <- makeCluster(mc <- threads)
stats <- parLapply(mc, mods, summary)
stopCluster(mc)
stats.coef <- lapply(1:length(stats), function(x) {
                         cbind("Coef" = rownames(coef(stats[[x]])), "rxn" = names(stats)[x], as.data.frame(coef(stats[[x]])))}
)
stats.coef <- data.table(do.call(rbind, stats.coef))[Coef == "outflux",]
stats.coef[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
stats.coef[,cpd := gsub(".*_(cpd[0-9]*)_.0","\\1", rxn)]
stats.coef <- merge(stats.coef, metabolites[,.(id, name)], by.x = "cpd", by.y= "id", all.x = TRUE)
#stats.coef[,rxn_id := gsub("_.0$","",rxn)]
#stats.coef <- merge(stats.coef, rxnsAnno[,.(id, name, definition)], by.x = "rxn_id", by.y ="id", all.x = TRUE, suffix = c("","_rxn"))

write.csv(file = "../results/outflux.GLMM.Response.intervention_coefs.csv", stats.coef, row.names = FALSE)


# create the diagnostic plots

rxns <- stats.coef[padj <= 0.05, rxn]
if(length(rxns) > 0){
pbar <- txtProgressBar(min = 0 , max = length(rxns), style = 3)
pdf("../results/outflux.GLMM.Response.intervention_sigDiagPlots.pdf", width = 8, height= 5)
i <- 0 
for(rxn in rxns){
    sim <- simulateResiduals(mods[[rxn]])
    plot(sim)
    mtext(rxn)
    i <- i+1
    setTxtProgressBar(pbar, i)
}
dev.off()
}


