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

## register cores for the parts which needs parallel execution
threads <- detectCores() - 1
registerDoMC(threads)

# go through the results for each intervention and cluster them
ints <- changes[,unique(intervention)]
clusters.all <- list()
models.all <- list()
coefs.all <- list()
i <- 0
n <- length(ints)
for(int in ints){
    fluxes <- dcast(changes[intervention == int,], sample~metabolite.name, value.var = "flux.change")
    fluxes.rn <- fluxes[,sample]
    fluxes <- as.matrix(fluxes[,-1])
    rownames(fluxes) <- fluxes.rn

    # create clusters for the fluxes
    ## remove near zero variance
    nzv <- nearZeroVar(fluxes)
    fluxes <- fluxes[,-nzv]

    # check if fluxes survived the NZV
    if(ncol(fluxes)>0){
        ## scale the values
        fluxes <- scale(fluxes)
        ## create a correlation matrix
        fluxes.cor <- cor(fluxes)
        # save the sign of the correlation for later       
        corr.sign <- melt(data.table(sign(fluxes.cor),
                                                keep.rownames = TRUE),
                                     id.var = "rn",
                                     variable.name = "rep.met",
                                     value.name = "sign.cor")
        colnames(corr.sign)[1] <- "met" 
        fluxes.cor <- 1-sqrt(fluxes.cor^2)

        set.seed(2368) # change in your code
        cluster <- dbscan(fluxes.cor, # the matrix with the correlation coefficients converted to distance
                          method = "dist", # telling dbscan not to calculate the distances by itself
                          eps = 0.1, # the threshold how distant members are allowed to be to each other to count a cluster
                          MinPts = 3) # the number of how many reactions have to make up a cluster at least - 3 is the minimum

        tbl.cluster <- data.table( cluster = cluster[["cluster"]],
                                  eps = cluster[["eps"]],
                                  MinPts = cluster[["MinPts"]],
                                  isseed = cluster[["isseed"]],
                                  met = rownames(fluxes.cor))

        ## get a representative reaction for each cluster
        tbl.cluster[cluster == 0 ,rep.met := met]
        tbl.cluster[cluster != 0, rep.met := names(which.min(rowSums(fluxes.cor[met,met]))), by = cluster]
        tbl.cluster[, cor.rep := fluxes.cor[met,met], by = met]
        tbl.cluster <- cbind(tbl.cluster, intervention = int)
        tbl.cluster <- merge(tbl.cluster, corr.sign)
        clusters.all[[int]] <- tbl.cluster

        mets <- tbl.cluster[,unique(rep.met)]
        # fit models
        mods <- foreach(met = mets) %dopar%{
            dat <- changes[intervention == int & metabolite.name == met,]
            mod <- lmer(data = dat,
                        formula = HB_Mayo_impu ~ source + seqtype + Age + flux.change + (1|PatientID))
        }
        names(mods) <- mets
        models.all[[int]] <- mods
        
        ## get stats
        mc <- makeCluster(mc <- threads)
        stats <- parLapply(mc, mods, summary)
        stopCluster(mc)

        stats.coef <- lapply(1:length(stats), function(x) {
                                 cbind("Coef" = rownames(coef(stats[[x]])), "met" = names(stats)[x], as.data.frame(coef(stats[[x]])))
                        }
        )
        stats.coef <- data.table(do.call(rbind, stats.coef))[Coef == "flux.change",]
        stats.coef[,padj := p.adjust(`Pr(>|t|)`, method = "BH")]
        stats.coef[,cpd := gsub("_o","", met)]
        stats.coef[,intervention := int]
        stats.coef <- merge(stats.coef, metabolites[,.(id, name)], by.x = "cpd", by.y= "id", all.x = TRUE)

        stats.coef <- merge(stats.coef, metabolites[,.(id, name)], by.x = "intervention", by.y= "id", all.x = TRUE, suffix = c("",".intervention"))
        coefs.all[[int]] <- stats.coef
        
    }
    # give some feedback during calculations
    i <- i+1
    cat("\r", round(i/n, digits = 2)*100, "%")

}
cat("\n")


# merge the individual tables to one large table and save the objects to the results
saveRDS(models.all, file = "../results/fluxChange.LMM.HBMayo.etiology_mods.RDS")
# models.all <- readRDS("../results/fluxChange.LMM.HBMayo.etiology_mods.RDS")

arguments <- clusters.all
arguments[["fill"]] <- TRUE # add this as an argument parameter to fill with NA in order to bind the tables, where no cluster could be detected
cluster <- data.table(do.call(rbind, arguments))
gz <- gzfile("../results/fluxChange_cluster.csv.gz","w")
write.csv(file = gz, cluster, row.names = FALSE)
close(gz)

coefs <- data.table(do.call(rbind, coefs.all))
coefs[,padj.all := p.adjust(`Pr(>|t|)`, method = "BH")] # adjust the pvalue once more, in order to get the 
gz <- gzfile("../results/fluxChange.LMM.HBMayo.etiology_coefs.csv.gz","w")
write.csv(file = gz, coefs, row.names = FALSE)
close(gz)
# coefs <- fread("../results/fluxChange.LMM.HBMayo.etiology_coefs.csv.gz")



vals <- coefs[padj <= 0.05, .(intervention, met)]
if(nrow(vals) > 0){
    pbar <- txtProgressBar(min = 0 , max = nrow(vals), style = 3)
    pdf("../results/absChange.LMM.HBMayo.etiology_sigDiagPlots.pdf", width = 8, height= 5)
    for(i in 1:nrow(vals)){
        intvn <- vals[i,intervention]
        met <- vals[i, met]
        mod <- models.all[[intvn]][[met]]
        plts <- list(~qqPlot(resid(mod),
                                main = paste0("Int:",intvn,"; Met:",met)),
                        plot(mod))
        print(plot_grid(plotlist = plts))
        setTxtProgressBar(pbar, i)
    }
    dev.off()
}
