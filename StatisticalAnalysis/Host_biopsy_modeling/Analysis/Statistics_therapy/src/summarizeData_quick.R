# Porthmeus
# 25.06.24

require(data.table)
require(ggplot2)
require(ggvenn)

res.dr <- file.path("..","results")
old.res.dr <- file.path("..","..","Statistics","results")
files.new <- list.files(res.dr, pattern = "_coefs.csv")
# add the old reference files
files.old <- list.files(old.res.dr, pattern = "intervention_coefs.csv")
files.dir <- c(file.path(res.dr, files.new), file.path(old.res.dr, files.old))
dats <- lapply(files.dir, fread)


sig.counts <- lapply(dats, function(x) length(x[padj <0.05,unique(rxn)]))
names(sig.counts) <- gsub("(.*)\\.[G]?LMM\\.(.*)_coefs.csv","\\1_\\2",c(files.new, files.old))
dat <- data.table(cnd = names(sig.counts),
                  counts = unlist(sig.counts))

dat[,c("tested.condition", "data.layer") := list(gsub("intervention","HBMayoCor",gsub(".*_(.*)","\\1",cnd)),
                                                 gsub("^(.*)_.*","\\1",cnd))]

bar.plt <- ggplot(dat, aes(x=tested.condition, y = counts))+
    geom_bar(stat="identity", position = "dodge", color = "black", aes(fill = data.layer))+
    theme_bw()+
    labs(x=element_blank(), y = "sig. rxns", title = "Biopsy")+
    theme(axis.text.x= element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(bar.plt, file = file.path(res.dr, "SigRxns_Barplot.pdf"))

# other barplot
dats <- lapply(files.dir, fread)
sig.rxns <- lapply(dats, function(x) x[padj <0.05,unique(rxn)])
names(sig.rxns) <- gsub("(.*)\\.[G]?LMM\\.(.*)_coefs.csv","\\1_\\2",c(files.new, files.old))
covariates <- c("Baseline","14d","14d0dDiff","TimeSeq","intervention","")
name_base <- gsub("\\.$","",c(paste("Remission",covariates, sep ="."),paste("Response",covariates, sep =".")))
bar.dat <- lapply(name_base, function(x) {
           sel <- names(sig.rxns) %in% paste0(c("FVA_","PA_","rxnExpr_"), x)
           length(unique(unlist(sig.rxns[sel])))
})
names(bar.dat) <- name_base
bar.dat <- data.table(count = unlist(bar.dat),
                      outcome = sapply(names(bar.dat), function(x) strsplit(x,split = "\\.")[[1]][1]),
                      covariate = sapply(names(bar.dat), function(x) strsplit(x,split = "\\.")[[1]][2]))
bar.dat[,covariate := gsub("14d","d14",covariate)]
bar.dat[,covariate := gsub("Baseline","d0",covariate)]
bar.dat[,covariate := gsub("intervention","HBMayoCor",covariate)]
bar.dat[is.na(covariate),covariate :="None"]


bar.plt2 <- ggplot(bar.dat, aes(x=covariate, y = count))+
    geom_bar(stat="identity",color = "black")+
    theme_bw()+
    facet_grid(~outcome)+
    labs(x="covariate/subset/interaction", y = "sig. rxns", title = "Biopsy")+
    theme(axis.text.x= element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(bar.plt2, file = file.path(res.dr, "SigRxns_Barplot2.pdf"), width = 4, height=4)

# create a vennDiagram
sel <- grep("(Remission|Response)($|\\.intervention)",name_base, value = TRUE)
venn.lst <- list()
for(nm in sel){
    sel2 <- grep(paste0(nm,"$"), names(sig.rxns))
    venn.lst[[nm]] <- unique(unlist(sig.rxns[sel2]))
}


venn.plts <- list()
for(outcm in c("Response","Remission")){
    venn.dat <- venn.lst[names(venn.lst) %in% paste0(outcm, c("",".intervention"))]
    names(venn.dat) <- gsub(outcm, "None", gsub(paste0(outcm, ".intervention"), "HBMayoCor", names(venn.dat)))
    venn.plts[[outcm]] <- ggvenn(venn.dat)+
        labs(title = outcm)
}

venn.plt <- cowplot::plot_grid(plotlist = venn.plts)

ggsave(venn.plt, file = file.path(res.dr, "VennPlot_sigRxns.pdf"), width = 12.5, height = 4.5)

#plot the permutation results
dat2 <- lapply(file.path(old.res.dr, files.old), fread)
sig.counts<- lapply(dat2, function(x) length(x[padj <0.05,unique(rxn)]))
names(sig.counts) <- gsub("(.*)\\.[G]?LMM\\.(.*)_coefs.csv","\\1_\\2",c(files.old))
dat2 <- data.table(cnd = names(sig.counts),
                  counts = unlist(sig.counts))
dat2[,c("tested.condition", "data.layer") := list(gsub(".intervention","",gsub(".*_(.*)","\\1",cnd)),
                                                 gsub("^(.*)_.*","\\1",cnd))]

#read permutation tests
perm_dat <- data.table()
for(i in 1:nrow(dat2)){
    lyr <- dat2[i, data.layer]
    cnd <- dat2[i, tested.condition]
    sig.rxns <- scan(paste0("../perm_results/",lyr,".", cnd,".sigRxns.txt"))
    perm_dat <- rbind(perm_dat, data.table(perm.sigs = sig.rxns,
               data.layer = lyr,
               tested.condition = cnd))
}
#dat2 <- merge(perm_dat, dat2, by = c("data.layer","tested.condition"))

plt <- ggplot(perm_dat, aes(y=perm.sigs, x = data.layer))+
    geom_boxplot(aes(fill = tested.condition))+
    theme_bw()+
    labs(x="Data layer",fill = "Therapy response", y = "No. of sig. rxn in permutation", title = "Label permutation")
ggsave(plt, file = "../results/SigRxns_afterPermutation.pdf", width = 5, height = 4)



# compare FVA results for remission
woCorrection <- fread(file.path(res.dr, "FVA.LMM.Remission_coefs.csv"))
oldCorrected <- fread(file.path(old.res.dr, "FVA.GLMM.Remission.intervention_coefs.csv"))
HBMayo <- fread(file.path(old.res.dr, "FVA.LMM.HBMayo.etiology_coefs.csv"))
sigRxns <- list(Remission_uncorrected = woCorrection[padj <0.05 , unique(rxn)],
                Remission_HBMayo_corr = oldCorrected[padj <0.05, unique(rxn)],
                HBMayo_only = HBMayo[padj <0.05, unique(rxn)])
ggvenn(sigRxns)


