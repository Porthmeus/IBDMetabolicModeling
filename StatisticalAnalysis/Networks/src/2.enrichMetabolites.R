# Porthmeus
# 26.10.23

# load the list of metabolites which are more often present in the significant data than expected and plot a venn across conditions

require(data.table)
require(ggplot2)
require(ggvenn)
require(RColorBrewer)
require(cowplot)
require(clusterProfiler)

# some variables
rsrc.dir <- file.path("..","resources")
net.dir <- file.path("..","results")

# rerun the analysis if the file does not exist already
#if(!file.exists(file.path(rsrc.dir, "metaboliteCountingStats.csv"))){
#    print("This might take a while...")
#    source("testMetabolitCounting.R")
#} 
#met.stats <- fread(file.path(rsrc.dir, "metaboliteCountingStats.csv"))

all_edges <- fread(file.path(net.dir,"baseEdges.csv"))
# remove compartments, exchange and transport reactions
sbst <- "\\[[a-z]\\]$"
rm.mets <- c("h","h2o","o2","pi","co2")
all_edges[,c("to.id","from.id") := list(gsub(sbst, "", to.id), gsub(sbst, "",from.id))]
all_edges <- all_edges[from.id != to.id,]
all_edges <- all_edges[!grepl("Transport|Exchange", subsystem),]
all_edges <- all_edges[!(from.id %in% rm.mets) & !(to.id %in% rm.mets),]
# melt down
met2rxn <- melt(all_edges[,.(from.id, to.id, reaction)], id.vars = "reaction", value.name = "met")
met2rxn[,variable := NULL]
met2rxn <- unique(met2rxn)

# get the network files and make an hgeo enrichment
tissues <- c("biopsy","blood")
sets <- c("HBMayo","Response","Remission")
network_all <- data.table()
hygeos <- list()
for(set in sets){
    for(tiss in tissues){
        #tiss <- "biopsy"
        dat <- fread(file.path(net.dir,paste0("basic.edge.table.1.", set,".",tiss,".csv")))
        network_all <- rbind(network_all, dat)
        hygeo <- enricher(gene = dat[,unique(rxn.id)],
                 TERM2GENE = met2rxn[,.(met,reaction)],
                 minGSSize = 3,
                 maxGSSize = 5000)
        hygeos[[paste(set,tiss, sep = ".")]] <- data.table(cbind("set" = set,
                                                    "tissue" = tiss,
                                                    data.table(hygeo@result)))
    }
}
hygeos <- do.call(rbind, hygeos)
hygeos[,p.adj := p.adjust(pvalue, method = "BH")]
hygeos[,sig := FALSE]
hygeos[p.adj < 0.05, sig := TRUE]

fwrite(hygeos, file = "../results/enrichedMetabolites.csv")

