# Porthmeus
# 20.02.22


require(data.table)
require(lme4)
require(lmerTest)
require(parallel)
require(doParallel)
require(foreach)
require(ggplot2)


# read the table
tab <- fread("../../../resources/totalFlux.longTab.csv.gz")
# correct numerical instabilities
tab[abs(flux)<1E-6, flux := 0]
# get the PatientID info
id.conv <- fread("../../../resources/IdentifierConversion.csv")
tab[,ID := gsub("[0-9]*_X","", sample)]
tab <- merge(tab, id.conv, by.x = "ID", by.y = "sample")
tab[,PatientID := paste(Person.Numbercode, Cohort, sep = "_")]
# scale the flux
tab[,log10_flux := sign(flux)*log10(abs(flux))]
tab[is.na(log10_flux), log10_flux := 0]
# relevel addRem to orig
tab[, unique(addRem)]
tab[, addRem := factor(addRem, levels = c("orig","double", "add","rem"))]



combs <- unique(tab[interventionMetabolite != "orig",.(target, interventionMetabolite)])

n <- nrow(combs)


getCoef <- function(i){
    tgt <- combs[i,target]
    inter <- combs[i,interventionMetabolite]
    dat <- tab[(interventionMetabolite == inter | interventionMetabolite == "orig") & target == tgt,]
    mod <- lmer(data = dat,
                formula = log10_flux  ~ source + seqtype + addRem + (1|PatientID))
    coefs <- coef(summary(mod))
    coefs <- cbind(target = tgt, intervention = inter, value = gsub("addRem","",rownames(coefs)), as.data.table(coefs))
    coefs <- coefs[grep("add|rem|double",value),]
    cat("\r",round((i/n)*100, digits = 2), "% done                     ")
    return(coefs)
}

threads <- detectCores() - 1
cl <- makeForkCluster(threads)
registerDoParallel(cl)

coefs <- foreach(i = 1:n) %dopar%  tryCatch(getCoef(i),
                                            error = function(e) {
                                                print(e)
                                                mod = NULL
                                            }
                                            )
stopCluster(cl)

coefs.all <- do.call(rbind, coefs[!is.null(coefs)])
coefs.all[, padj := p.adjust(`Pr(>|t|)`,method = "BH")]
gz = gzfile("../results/LMM.InterventionsTest_coefs.csv.gz", "w")
write.csv(file = gz,
          coefs.all,
          row.names = FALSE)
close(gz)
