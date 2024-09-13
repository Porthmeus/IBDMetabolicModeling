
# Porthmeus
# 26.09.22

# creates LMMs for changes in the production of metabolites due to the virtual addition of 1 gram in the simulation of the communities

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


changes <- fread("../../../resources/absoluteChange.longTab.csv.gz")
# remove fluxes which are not available to the host
changes <- changes[grep("_o$", metabolite.name),]

clinic <- fread("../../../resources/ClinicalDataMergeLongImputeCor.csv")
idconv <- fread("../../../resources/IdentifierConversion.csv")
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")

# add the sample information
changes[, SeqID := gsub("^X","",sample)]
changes <- merge(changes, idconv, by.x ="SeqID", by.y = "sample")
changes[, PatientID_Time := paste(Cohort, Person.Numbercode, Time_seq, sep = "_")]
changes <- merge(changes, clinic, by = "PatientID_Time", all.x = TRUE)
# correct the remission to be a factor and correctly spelled
changes[,Remission := factor(Remision)]
changes[,Response := factor(Responder)]

# read cluster
cluster <- fread("../results/fluxChange_cluster.csv.gz")

## register cores for the parts which needs parallel execution
threads <- detectCores() - 1
registerDoMC(threads)


# go through the results for each intervention and cluster them
ints <- cluster[, unique(intervention)]
models.all <- list()
errors <- list()
coefs.all <- list()
i <- 0
n <- length(ints)
for(intv in ints){
        
        # getting the relevant flux changes from the clustering
        tbl.cluster <- cluster[intervention == intv,]
        mets <- tbl.cluster[,unique(rep.met)]
        
        # fit models
        mods <- foreach(met = mets) %dopar%{
#        mods <- list(); for(met in mets){
            dat <- changes[intervention == intv & metabolite.name == met & Response != "C",]
            mod <- tryCatch({m <- glmer(data = dat,
                        formula = Response ~ source + seqtype + Age + flux.change + (1|PatientID),
                        family = "binomial")
                            },
                            error = function(e) {
                               # print(paste0("Error to build the model for ", met))
                               # print(paste0(e))
                                msg <- paste("Error to build the model for", met, "\n", e)
                                return(msg)
                            }
                            )
#            mods[[met]] <- mod
        }
        names(mods) <- mets
        sel <- sapply(mods, is.character)
        errors[[intv]] <- mods[sel]
        mods <- mods[!sel]
        models.all[[intv]] <- mods
        
        ## get stats
        mc <- makeCluster(mc <- threads)
        stats <- parLapply(mc, mods, summary)
        stopCluster(mc)

        stats.coef <- lapply(1:length(stats), function(x) {
                                 cbind("Coef" = rownames(coef(stats[[x]])), "met" = names(stats)[x], as.data.frame(coef(stats[[x]])))
                        }
        )
        stats.coef <- data.table(do.call(rbind, stats.coef))[Coef == "flux.change",]
        stats.coef[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
        stats.coef[,cpd := gsub("_o","", met)]
        stats.coef[,intervention := intv]
        stats.coef <- merge(stats.coef, metabolites[,.(id, name)], by.x = "cpd", by.y= "id", all.x = TRUE)

        stats.coef <- merge(stats.coef, metabolites[,.(id, name)], by.x = "intervention", by.y= "id", all.x = TRUE, suffix = c("",".intervention"))
        coefs.all[[intv]] <- stats.coef
    # give some feedback during calculations
    i <- i+1
    cat("\r", round(i/n, digits = 2)*100, "%")
    }
cat("\n")
print(errors)


# merge the individual tables to one large table and save the objects to the results
saveRDS(models.all, file = "../results/fluxChange.GLMM.Response.intervention_mods.RDS")
# models.all <- readRDS("../results/fluxChange.GLMM.Remission.intervention_mods.RDS")

coefs <- data.table(do.call(rbind, coefs.all))
coefs[,padj.all := p.adjust(`Pr(>|z|)`, method = "BH")] # adjust the pvalue once more, in order to get the 
gz <- gzfile("../results/fluxChange.GLMM.Response.intervention_coefs.csv.gz","w")
write.csv(file = gz, coefs, row.names = FALSE)
close(gz)
# coefs <- fread("../results/fluxChange.GLMM.Remission.intervention_coefs.csv.gz")



vals <- coefs[padj <= 0.05, .(intervention, met)]
if(nrow(vals) > 0){
    pbar <- txtProgressBar(min = 0 , max = nrow(vals), style = 3)
    pdf("../results/absChange.GLMM.Response.intervention_sigDiagPlots.pdf", width = 8, height= 5)
    for(i in 1:nrow(vals)){
        intvn <- vals[i,intervention]
        met <- vals[i, met]
        mod <- models.all[[intvn]][[met]]
        sim <- simulateResiduals(mod)
        plot(sim)
        mtext(paste0("Int.:",intvn,"; Met.:", met))
        setTxtProgressBar(pbar, i)
    }
    dev.off()
}
