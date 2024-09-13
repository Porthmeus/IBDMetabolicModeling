# Porthmeus
# 18.11.22

# make the crosslinks between the flux data of the bacteria and the host metabolism

# libraries
require(data.table)
require(lme4)
require(lmerTest)
require(car)
require(ggplot2)
#require(parallel)
require(doMC)
require(foreach)
require(cowplot)

# load flux data
fluxes <- read.csv("../../../resources/future.eMed_flux.csv", row.names = 1)

clinic <- fread("../../../resources/ClinicalDataMergeLongImputeCor.csv")
idconv <- fread("../../../resources/IdentifierConversion.csv")
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")
#pwys <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/meta_rea_pwy-gapseq.tbl")
#pwys <- pwys[,.(rxn_id = unlist(strsplit(gapseq, split = " "))), by = .(pathway.name, pathway)]
#pwys <-pwys[!duplicated(pwys),]
#pwys <- pwys[,.(name = paste(pathway.name, collapse = ";"), pathway = paste(pathway, collapse = ";")), by = rxn_id]
#rxnsAnno <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_reactions.tsv")

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


# get the host data
sets <- c(FVA.center = "../../../Pipeline_MM/results/FVA/centerFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv",
          FVA.range = "../../../Pipeline_MM/results/FVA/rangeFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv",
          PA = "../../../Pipeline_MM/results/ModelAnalysis/rxnIdMat.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv",
          rxnExpr = "../../../Pipeline_MM/results/RxnExpression/rxnExpr.GL10-L50-GU90.colormore3D.csv")
sets.host <- list()
for(setName in names(sets)){
    set <- sets[[setName]]
    mat <- read.csv(set, row.names= 1)
    cluster.tbl <- fread(paste0("../../../Pipeline_MM/manualAnalysis/Statistics/results/",gsub("\\..*$","",setName), "_DBSCAN_Cluster.csv"))
    rxn <- cluster.tbl[,unique(rep.rxn)]
    mat <- mat[rxn,]
    mat <- melt(data.table(mat, keep.rownames=TRUE),
                id.vars = "rn",
                variable.name = "Sample")
    mat[,setName := setName]
        ## correct the True/False to 1/0
    if(setName == "PA"){
        mat[,value := c(0,1)[(value == "True")+1]]
    }
    sets.host[[setName]] <- mat
}
sets.host <- do.call(rbind, sets.host)

# get the sample information
metaRNA <- fread("../../../resources/MetaData_merged_basicPatientCor_Biopsy.csv")
metaRNA[,PatientID_TimeSeq := paste(PatientID, Time_seq, sep ="_")]
sets.host <- merge(sets.host, metaRNA[,.(SeqID,PatientID_TimeSeq)], by.x = "Sample", by.y = "SeqID", all.x = TRUE)

# combine both data sets
## combining the dataset completely would exceed memory possibilities, thus I unify the design of the tables and merge subsets during the iteration process.
fluxes.merged <- fluxes.merged[,.(PatientID_Time,Rxn = ex.rxn, flux, seqtype, source, Gender, Diagnosis, Age, HB_Mayo_impu)]
fluxes.merged[,PatientID := gsub("(^.*_[0-9]*)_.*","\\1",PatientID_Time)]
sets.host <- sets.host[,.(PatientID_Time = PatientID_TimeSeq, Rxn = rn, value, setName)]

# create for each bac.rxn host.rxn pair an LMM
rxns.host <- sets.host[,unique(Rxn)]
rxns.bac <- fluxes.merged[,unique(Rxn)]
rxns.both <- data.table(
                        bac.rxn = sort(rep(rxns.bac, length(rxns.host))),
                        host.rxn = rep(rxns.host, length(rxns.bac))
                        )
rxns.both[,fl_base := paste0(bac.rxn, ".", host.rxn)]

#out.dir <- "../results/flux.LMM.hostRxns.unfiltered_crosslinks/"
#if(!dir.exists(out.dir)){
#    dir.create(out.dir)
#}

# register cores
threads <- detectCores()-1
registerDoMC(threads)

# testing purposes 
#rxns.both <- rxns.both[sample(1:nrow(rxns.both), 100)]

pvalues <- foreach(i = 1:nrow(rxns.both)) %do% {
    dat <- merge(fluxes.merged[Rxn == rxns.both[i,bac.rxn]],
                 sets.host[Rxn == rxns.both[i, host.rxn]],
                 by = "PatientID_Time",
                 suffix = c(".bac",".host"),
                 allow.cartesian=TRUE)
    if(length(dat[,unique(setName)])>2){
        mod <- lmer(data = dat,
                    formula = flux ~ source + seqtype + Gender + Diagnosis + Age + HB_Mayo_impu + value:setName+ (1| PatientID))
        setnames <- names(coef(mod)[[1]])
        setnames <- gsub("value:setName","",setnames[grepl("setName",setnames)])
    } else {
        mod <- lmer(data = dat,
                    formula = flux ~ source + seqtype + Gender + Diagnosis + Age + HB_Mayo_impu + value+ (1| PatientID))
        setnames <- dat[,unique(setName)]
    }
    coefs <- coef(summary(mod))
    pvals <- coefs[grepl("value", rownames(coefs)), ncol(coefs)]
    estimates <- coefs[grepl("value", rownames(coefs)), 1]
    #
    #if(any(pvals < 0.05)){
    #    saveRDS(mod,
    #            file = file.path(out.dir, paste0(rxns.both[i, fl_base], ".RDS")))
    #}
    cat("\r",round(i/nrow(rxns.both)*100, digits = 2) ,"% done")
    pvals <- data.table(host.rxn = rxns.both[i,host.rxn],
                        bac.rxn = rxns.both[i,bac.rxn],
                        pvals = pvals,
                        setName = setnames)
}
pvalues <- do.call(rbind, pvalues)
pvalues[,padj := p.adjust(pvals, method = "BH")]
gz <- gzfile(description = "../results/flux.LMM.hostRxns.unfiltered_crosslinks.csv.gz", "w")
write.csv(pvalues,
          row.names = FALSE,
          file = gz)
close(gz)

