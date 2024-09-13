# libraries
## @knitr loadLibraries
require(data.table)
require(foreach)
require(parallel)
require(doMC)
require(clusterProfiler)


# load data
## @knitr loadData_stat
subsystems <- fread("../dat/subsystems.csv")

hgeos<- list()
for(dat.set in c("rxnExpr","PA","FVA")){
    # load the data
    cluster <- fread(paste0("../results/",dat.set,"_DBSCAN_Cluster.csv"))
    stats <- fread(paste0("../results/",dat.set,".LMM.Remission_coefs.csv"))
    #  expand the cluster
    stats.exp <- merge(stats, cluster, by.x = "rxn", by.y = "rep.rxn")
    sig.rxns <- unique(stats.exp[padj<0.05,rxn.y])
    # make the enrichment
    hygeo <- enricher(sig.rxns,
                      TERM2GENE = subsystems)
    if(!is.null(hygeo)){
        hygeo@result <- cbind(hygeo@result, data.set = dat.set)
    }
    hgeos[[dat.set]] <- hygeo
}

hgeos <- lapply(hgeos, function(x) data.table(x@result[,-8]))
hgeos <- do.call(rbind, hgeos)
hgeos[,p.adjust := p.adjust(pvalue, method = "BH")]
fwrite(hgeos, "../results/HygeoTests_all.csv")
