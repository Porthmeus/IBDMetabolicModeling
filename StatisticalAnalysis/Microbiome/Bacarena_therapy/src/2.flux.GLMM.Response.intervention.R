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
fluxes <- read.csv("../../../resources/future.eMed.bacarena_flux.csv", row.names = 1)

clinic <- fread("../../../resources/ClinicalDataMergeLongImputeCor.csv")
idconv <- fread("../../../resources/IdentifierConversion.csv")

metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")
pwys <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/meta_rea_pwy-gapseq.tbl")
pwys <- pwys[,.(rxn_id = unlist(strsplit(gapseq, split = " "))), by = .(pathway.name, pathway)]
pwys <-pwys[!duplicated(pwys),]
pwys <- pwys[,.(name = paste(pathway.name, collapse = ";"), pathway = paste(pathway, collapse = ";")), by = rxn_id]
rxnsAnno <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_reactions.tsv")

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

# remove control patients and reformat Response
fluxes.merged <- fluxes.merged[Responder != "C",]
fluxes.merged <- fluxes.merged[,Response:= factor(Responder, levels =c("NR","R"))]

# fit the models
## register cores
getGLMMstats <- function(rxn){
    dat <- fluxes.merged[ex.rxn == rxn,]#& seqtype != "cDNA",]
    mod <- tryCatch({glmer(data = dat,
                formula = Response ~ source + seqtype + Age + flux + (1|PatientID),
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
write.csv(file = "../results/flux.Response_coefs.csv", stats.coef, row.names = FALSE)

# do the simple linear models
dats <- list()
dats[["d0"]] <- fluxes.merged[Response != "C" & Time_seq.x == 0, .(Response, source, seqtype, Age, value = flux, PatientID, ex.rxn)]
dats[["d14"]] <- fluxes.merged[Response != "C" & Time_seq.x == 14, .(Response, source, seqtype, Age, value = flux, PatientID, ex.rxn)]
# calculate difference between time points
dat.temp <- fluxes.merged[Response != "C" & Time_seq.x %in% c(0, 14), .(Response, source, seqtype, Age, value = flux, PatientID, ex.rxn, Time_seq = Time_seq.x)]
dat.mat <- dcast(data=dat.temp, Response+Age+source+seqtype+PatientID+ex.rxn~paste0("d",Time_seq), value.var = "value")
dat.mat[, value := d14-d0]
dat.mat <- dat.mat[!is.na(value),]
dats[["d14d0Diff"]] <- dat.mat

getModelStats <- function(rxn){
    dat.temp <- dat[ex.rxn == rxn,]
    mod <- tryCatch({glm(data = dat.temp,
               formula = factor(Response) ~ source + seqtype + Age + value,
               family = "binomial")},
                    error = function(e){
                        print(e)
                        return(NA)
                    }
                    )
    coefs <- coef(summary(mod))
    if("value" %in% rownames(coefs)){
        stats <- cbind(rxn = rxn, as.data.frame(coefs)["value",])
        return(stats)
    } else{
        return(NA)
    }
}


for(dat.n in names(dats)){
    dat <- dats[[dat.n]]
    cl <- makeForkCluster(threads)
    rxns <- dat[,unique(ex.rxn)]
    stats <- pblapply(rxns, getModelStats, cl=cl)
    stopCluster(cl)
    stats <- data.table(do.call(rbind, stats[!is.na(stats)]))
    stats[,c("outcome", "coef","padj") := list("Response",dat.n, p.adjust(`Pr(>|z|)`, method = "BH"))]
    fwrite(stats, file = paste0("../results/flux.Response.",dat.n,"_coefs.csv"))
}

# create the diagnostic plots

#rxns <- stats.coef[padj <= 0.05, rxn]
#if(length(rxns) > 0){
#pbar <- txtProgressBar(min = 0 , max = length(rxns), style = 3)
#pdf("../results/flux.GLMM.Response.intervention_sigDiagPlots.pdf", width = 8, height= 5)
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


