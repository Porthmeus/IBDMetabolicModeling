# Porthmeus
# 26.09.22

# test if the different interventions cause a significant change in the fluxes from 0

require(data.table)
require(caret)
require(parallel)
require(doMC)
require(foreach)


changes <- fread("../../../resources/absoluteChange.longTab.csv.gz")
# remove fluxes which are not available to the host
changes <- changes[grep("_o$", metabolite.name),]

# add the names of the metabolites
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")

# get the threads and register workers for multthreading
threads <- detectCores() - 1
registerDoMC(threads)

# go through each of the metabolites and check difference to 0
ints <- changes[,unique(intervention)]
coefs.all <- list()
i <- 0
n <- length(ints)
for(int in ints){
    fluxes <- dcast(changes[intervention == int,], sample~metabolite.name, value.var = "flux.change")
    fluxes.rn <- fluxes[,sample]
    fluxes <- as.matrix(fluxes[,-1])
    rownames(fluxes) <- fluxes.rn

    ## remove near zero variance
    nzv <- nearZeroVar(fluxes)
    mets <- colnames(fluxes)[-nzv]

    tests <- foreach(met = mets, .combine = "rbind") %dopar% {
        vec <- changes[intervention == int & metabolite.name == met, flux.change]
        test <- t.test(vec)
        data.frame(cpd = met,
                   intervention = int,
                   mean = mean(vec),
                   sd = sd(vec),
                   median = median(vec),
                   mad = mad(vec),
                   stderr = test[["stderr"]],
                   pvalue = test[["p.value"]])
    }
    
    coefs.all[[int]] <- tests
    i=i+1
    cat("\r", round(i/n, digits = 2) * 100, "%")
}
cat("\n")

coefs <- data.table(do.call(rbind,coefs.all))
coefs[,padj := p.adjust(pvalue, method = "BH")]
write.csv(file = gzfile("../results/fluxChange.testVs0.csv.gz","w"), coefs, row.names = FALSE)
