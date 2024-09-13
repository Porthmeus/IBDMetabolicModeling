# Porthmeus
# 25.10.23

# the idea is to count the how often metabolites are affected in the sig. reactions and whether this is more often than expected

require(data.table)
require(ggplot2)

testMetAbundance <- function(sig.rxns, mets2rxns, nperm = 999){
    # sig.rxns - a vector of all significant reactions
    # met2rxns - a data.frame with the first column containing the metabolites, the second column containing the reactions. Reaction names of sig.rxns and met2rxn must comply.
    
    # reformat mets2rxns to have expected behavior
    mets2rxns <- data.table(met = mets2rxns[[1]],
                            rxn = mets2rxns[[2]])
    # start with counting the metabolites affected
    sig.mets.sub <- table(mets2rxns[rxn %in% sig.rxns,met])
    
    # generate a dummy vector for the sampling
    nmets <- length(mets2rxns[,unique(met)])
    dummy.mets <- rep(0, nmets)
    names(dummy.mets) <- sort(mets2rxns[,unique(met)])
    # add the 0 to the sig metabolites
    sig.mets <- dummy.mets
    sig.mets[names(sig.mets.sub)] <- sig.mets.sub
    
    # get the unique rxns in the network
    all.rxns <- mets2rxns[,unique(rxn)]
    n <- length(sig.rxns)
    # create the vector to store the p values
    p.vec <- dummy.mets
    # do the testing
    for(i in 1:nperm){
        cat("\r", round(i/nperm, digits = 2)*100, "%               ")
        # reset the the test vector
        test.mets <- dummy.mets
        test.rxns <- sample(all.rxns, n)   
        test.mets.sub <- table(mets2rxns[rxn %in% test.rxns,met])
        test.mets[names(test.mets.sub)] <- test.mets.sub
        p.vec <- p.vec + (sig.mets > test.mets)
    }
    
    p.vec[p.vec == 0] <- 1
    p.vec <- 1-p.vec/nperm
    return(p.vec[names(sig.mets.sub)])
}

rsrc.dir <- file.path("..","data","Host")

# get the list of metabolites in the rections
met2rxn <- fread(file.path(rsrc.dir, "Mets2Subsystems.csv"))

# some variable to rely on
prfx <- c("PA","rxnExpr","FVA")

final.dat <- data.table()
# load the significant reactions
for(tissue in c("biopsy","blood")){
    #tissue <- "biopsy"
    
    # load the clustering data
    files <- sapply(prfx, function(x) list.files(file.path(rsrc.dir, tissue), pattern = paste0(x,".*Cluster.csv")))
    cluster <- lapply(names(files), function(x) cbind(data = x, fread(file.path(rsrc.dir, tissue, files[[x]]))))
    names(cluster) <- prfx

    for(set in c("HBMayo","Remission","Response")){
        #set <- "HBMayo"
        files <- sapply(prfx, function(x) list.files(file.path(rsrc.dir, tissue), pattern = paste0(x,".*",set)))
        coefs <- lapply(names(files), function(x) cbind(data = x, fread(file.path(rsrc.dir, tissue, files[[x]]))))
        names(coefs) <- prfx

        # expand the clusters
        for(prf in prfx){
            coefs[[prf]] <- merge(cluster[[prf]], coefs[[prf]], by.y = "rxn", by.x ="rep.rxn")
        }

        # get the significant reactions
        sig.rxn <- unique(unlist(sapply(coefs, function(x) x[padj < 0.05, rxn])))
        #sig.rxns <- sig.rxn
        # remove the compartment annotation from met2rxn
        met2rxn[,met.id.s := gsub("\\[[a-z]\\]","",metabolite)]

        sig.mets <- testMetAbundance(sig.rxn = sig.rxn, 
                                     mets2rxns = met2rxn[,.(met.id.s, reaction)],
                                     nperm = 9999)

        stats <- data.table(tissue = tissue,
                            set = set,
                            met = names(sig.mets),
                            p.value = sig.mets,
                            p.adj = p.adjust(sig.mets, method = "BH"))
        final.dat <- rbind(final.dat, stats)
    }
}
write.csv(final.dat, file = file.path(rsrc.dir, "metaboliteCountingStats.csv"), row.names = FALSE)
