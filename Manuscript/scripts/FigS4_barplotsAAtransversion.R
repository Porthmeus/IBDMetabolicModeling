# Porthmeus
# 15.02.24

require(data.table)
require(ggplot2)
require(igraph)
require(ggraph)
require(ggvenn)

dd.network <- file.path("..","data","Networks")
dd.host <- file.path("..","data","Host")
dd.figs <- file.path("..", "figures_raw")
# load base edges
base.edges <- fread(file.path(dd.network, "baseEdges.csv"))
base.edges[,c("from.id","to.id") := list(gsub("\\[[a-z]\\]$","",from.id),
                                         gsub("\\[[a-z]\\]$","",to.id))]

# get all amino acid transversion reactions
aas <- c("arg_L","his_L","lys_L","asp_L","glu_L","ser_L","thr_L","asn_L","gln_L","cys_L","gly_L","pro_L","ala_L","val_L","ile_L","leu_L","met_L","phe_L","tyr_L","trp_L")
rxns.dat <- base.edges[from.id %in% aas & to.id %in% aas & from.id != to.id, ]
rxns <- rxns.dat[c(-grep("Protein formation", subsystem), -grep("Transport", subsystem)),unique(reaction)]

# make a small enrichment test
## load data from associations
allRxns <- base.edges[,unique(reaction)]
for(tis in c("biopsy","blood")){
rxnExpr <- fread(file.path(dd.host, tis, "rxnExpr.LMM.HBMayo.etiology_coefs.csv"))
rxnExpr.cluster <- fread(file.path(dd.host, tis, "rxnExpr_DBSCAN_Cluster.csv"))
rxnExpr <- merge(rxnExpr.cluster[,.(rxn, rep.rxn,sign.cor)], rxnExpr, by.x = "rep.rxn", by.y = "rxn")
rxnExpr[,Estimate := Estimate * sign.cor]
PA <- fread(file.path(dd.host, tis, "PA.LMM.HBMayo.etiology_coefs.csv"))
PA.cluster <- fread(file.path(dd.host, tis, "PA_DBSCAN_Cluster.csv"))
PA <- merge(PA.cluster[,.(rxn, rep.rxn,sign.cor)], PA, by.x = "rep.rxn", by.y = "rxn")
PA[,Estimate := Estimate * sign.cor]
FVA <- fread(file.path(dd.host, tis, "FVA.LMM.HBMayo.etiology_coefs.csv"))
FVA.cluster <- fread(file.path(dd.host, tis, "FVA_DBSCAN_Cluster.csv"))
FVA <- merge(FVA.cluster[,.(rxn, rep.rxn,sign.cor.range, sign.cor.center)], FVA, by.x = "rep.rxn", by.y = "rxn")
FVA[coef == "range",Estimate := Estimate * sign.cor.range]
FVA[coef == "center",Estimate := Estimate * sign.cor.center]
# get all and all significant reactions for the tissue
sigRxns <- unique(c(rxnExpr[padj <0.05,rxn],
                 PA[padj <0.05, rxn],
                 FVA[padj <0.05, rxn]))
# create the contigency table
#                   |significant|not significant|
#  in AAtransversion|    val1   |      val2     |
# not AAtransversion|    val3   |      val4     |
val1 <- sum(rxns %in% sigRxns)
val2 <- length(rxns) - val1
val3 <- sum(!(sigRxns %in% rxns))
val4 <- length(allRxns) - sum(c(val1,val2,val3))
# sanity checks:
stopifnot(sum(val1,val2,val3,val4) == length(allRxns))
cont.tab <- matrix(c(val1,val2,val3,val4),ncol = 2, byrow=TRUE)
res <- fisher.test(cont.tab, alternative = "greater")
print(tis)
print(res)
}

st <- "HBMayo"
aaTrans.rxns <- data.table()
for(tis in c("blood","biopsy")){
    edges <- fread(file.path(dd.network, paste0("basic.edge.table.1.",st,".",tis,".csv")))
    aaTrans.rxns <-rbind(aaTrans.rxns, edges[rxn.id %in% rxns &
                      from.id %in% aas &
                      to.id %in% aas,
                  .(tissue = tis, from_to = paste0(paste(unique(gsub("_L$","",from.id)), collapse = " + "),
                                     " -> ",
                                     paste(unique(gsub("_L$","",to.id)), collapse = " + "))),
                  by=.(rxn.id, reac.name,  EffSize)])
}

# sort data 
aaTrans.rxns[,from_to_rxn := paste0(from_to, " (",rxn.id,")")]
aaTrans.od <- aaTrans.rxns[,.(od = mean(EffSize)), by = from_to_rxn]
aaTrans.od <- aaTrans.od[order(od), from_to_rxn]
aaTrans.rxns[,from_to_rxn := factor(from_to_rxn, levels = aaTrans.od)]

bar.plt <- ggplot(aaTrans.rxns, aes(y= EffSize, x = from_to_rxn)) +
    geom_bar(stat= "identity", aes(fill = EffSize), color = "black") +
    scale_fill_gradient2(low = "#c51b7d", mid = "white", high = "#4d9221")+
    geom_hline(yintercept = 0)+
    facet_grid(tools::toTitleCase(tissue)~.) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "left")+
    labs(y= "Mean effect size", x = "Reaction", fill ="Eff. Size")
#bar.plt
ggsave(bar.plt, file = file.path(dd.figs, "AAtransversion.pdf"),
       width = 7,
       height =5
       )

#venn.dat <- list("All" = unique(rxns),
#                 "Biopsy" = unique(aaTrans.rxns[tissue == "biopsy", rxn.id]),
#                 "Blood" = unique(aaTrans.rxns[tissue == "blood", rxn.id])
#                 )
#ggvenn(venn.dat)

# create a similar plot for NAD metabolism
mets <- c("quln","nad","nadp","nicrns","nmn","ncam","nicrnt","nac","rnam","nadh","dnad","1mncam","nadph")
nad.edges <- base.edges[grepl("NAD metabolism", subsystem) & from.id %in% mets & to.id %in% mets,]
rxns <- nad.edges[,unique(reaction)]

st <- "HBMayo"
nad.rxns <- data.table()
for(tis in c("blood","biopsy")){
    edges <- fread(file.path(dd.network, paste0("basic.edge.table.1.",st,".",tis,".csv")))
    nad.rxns <-rbind(nad.rxns, edges[rxn.id %in% rxns &
                      from.id %in% mets &
                      to.id %in% mets,
                  .(tissue = tis, from_to = paste0(paste(unique(gsub("_L$","",from.id)), collapse = " + "),
                                     " -> ",
                                     paste(unique(gsub("_L$","",to.id)), collapse = " + "))),
                  by=.(rxn.id, reac.name,  EffSize)])
}


nad.rxns[,from_to_rxn := paste0(from_to, " (",rxn.id,")")]
nad.od <- nad.rxns[,.(od = mean(EffSize)), by = from_to_rxn]
nad.od <- nad.od[order(od), from_to_rxn]
nad.rxns[,from_to_rxn := factor(from_to_rxn, levels = rev(nad.od))]

bar.plt <- ggplot(nad.rxns, aes(y= EffSize, x = from_to_rxn)) +
    geom_bar(stat= "identity", aes(fill = EffSize), color = "black") +
    scale_fill_gradient2(low = "#c51b7d", mid = "white", high = "#4d9221")+
    geom_hline(yintercept = 0)+
    facet_grid(tools::toTitleCase(tissue)~.) +
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
    labs(y= "Mean effect size", x = "Reaction")
bar.plt

# create a graph for the NAD metabolism
od <- unique(c("from.id","to.id", colnames(nad.edges)))
nad.graph <- graph_from_data_frame(nad.edges[,..od])
plot(nad.graph)
