# Porthmeus
# 9.02.23

require(data.table)
require(ggplot2)
require(ComplexUpset)
require(cowplot)
require(clusterProfiler)
require(RColorBrewer)

# plot the association of serum metabolites with the disease activity, response and remission
ddr <- file.path("..","data","Metabolomics")
ddrHost <- file.path("..","data","Host")

# Annotations
## metabolites general information from gapseq
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites_edited.tsv")
setkey(metabolites, "id")
## metabolomics to Recon3D tanslation
met.meta <- fread(file.path(ddr,"Metabolomics2Recon3D.csv"))
## metabolomics category information
met.categories <- fread(file.path(ddr, "metabolome_names.csv"))
## combine all information
met.meta <- merge(met.meta, met.categories[,.(Abbreviation, Category2,Category)], by = "Abbreviation", all.x = TRUE)
met.meta[,colAbbr := gsub("(^\\d)","X\\1",gsub("[\\(|\\.|:| |\n |\\)|/|-]",".",Abbreviation))]
met.meta <- met.meta[!duplicated(colAbbr),]
met.meta[,Name_metabolomics := Name]
## get information on which metabolite is present in which subsystem of Recon3D
met2subsystem <- fread(file.path(ddrHost, "Mets2Subsystems.csv"))
met2subsystem[,metabolite.woComp := gsub("\\[[a-z]\\]","",metabolite)]
## get information which subsystem was significanctly enriched in the host reaction analysis
sig.subs.df <- fread(file.path(ddrHost, "EnrichedSubsystems.csv"))

# expand the set by the clustered reactions
cluster <- fread(file.path(ddr,"metabolites_cluster.csv"))


coef.fls <- list.files(ddr, pattern = "_coefs.csv.*$")
coef.fls <- coef.fls[-grep("crosslinks", coef.fls)]

## create the heatap
### load the data
coefs <- lapply(coef.fls, function(set){
                set.title <- strsplit(set,
                             split = "\\.")[[1]][3]
                coef.tb <- fread(file.path(ddr, set))[,.(metabolite, Estimate, padj)]
                coef.tb <- cbind(set = set.title,
                                 coef.tb)
})
coefs <- do.call(rbind, coefs)
coefs.exp <- merge(coefs, cluster[,.(met, rep.met, sign.cor)], by.x= "metabolite", by.y = "rep.met", suffix = c(".rep",".all"), allow.cartesian=TRUE)
# correct the sign of the estimate by the sign of the correlation in the cluster
coefs.exp[,Estimate := Estimate*sign.cor]
coefs.exp <- merge(coefs.exp[padj < 0.05,], met.meta[,.(colAbbr,Recon3D)], by.x = "metabolite", by.y = "colAbbr", all.x = TRUE)
coefs.exp <- merge(coefs.exp, unique(met2subsystem[,.(subsystem,metabolite.woComp)]), by.x = "Recon3D", by.y= "metabolite.woComp")

### make hygeo enrichments of the metabolites in the subsystems
sub2mets <- fread(file.path("..","data","Networks","baseEdges.csv"))
sub2mets <- melt(sub2mets[,.(from.id,to.id,subsystem)], id.vars = "subsystem", value.name = "met.id")[,.(subsystem,met.id)]
sub2mets[,met := gsub("(.*)\\[[a-z]\\]$","\\1",met.id)]
sub2mets <- unique(sub2mets[,.(subsystem, met)])
enrichments <- list()
for(ss in coefs.exp[,unique(set)]){
    sig.mets <- unique(coefs.exp[set == ss,Recon3D])
    hygeo <- enricher(sig.mets,
                      TERM2GENE = sub2mets)
    result <- data.table(hygeo@result)
    result <- result[p.adjust <0.05, .(subsystem = ID, rxns = geneID,GeneRatio, BgRatio, p.adjust, pvalue, set = ss)]
    enrichments[[ss]] <- result
}
enrichments <- do.call(rbind, enrichments)


### count the metabolites in all subsystems and which are significantly changed
subsys.sizes <- unique(met2subsystem[,.(subsystem, metabolite.woComp)])[,.(subsys.size = .N), by = subsystem]
coefs.exp.plt <- unique(coefs.exp[,.(Recon3D,subsystem,set)])[,.(metsInSubsys = .N), by = .(subsystem, set)]
coefs.exp.plt <- merge(coefs.exp.plt, subsys.sizes)
coefs.exp.plt[,affected := metsInSubsys/subsys.size]

### cluster the subsystems hierarchically for a nicer plot
order.tb <- dcast(coefs.exp.plt, subsystem~set, value.var = "affected", fill = 0)
od <- hclust(dist(as.matrix(order.tb[,-1])))[["order"]]
od <- order.tb[od,subsystem]
coefs.exp.plt[,subsystem := factor(subsystem, levels = od)]

# add mock tissue
coefs.exp.plt[,tissue := "blood"]
coefs.exp.plt.blood <-  copy(coefs.exp.plt)
coefs.exp.plt <- rbind(coefs.exp.plt[,tissue := "biopsy"], coefs.exp.plt.blood)

### add the information on which subsystem was significantly enriched in Recon3D
coefs.exp.plt <- merge(coefs.exp.plt,
                       sig.subs.df[,.(subsystem, set, tissue, enriched.host = 1)],
                       by = c("set","subsystem","tissue"),
                       all.x = TRUE)
coefs.exp.plt[is.na(enriched.host),enriched.host := 0]

### add the subsystem enrichment of the metabolomics data
coefs.exp.plt <- merge(coefs.exp.plt, enrichments[,.(subsystem, set, enriched.mtblmcs = 2)], by =c("subsystem","set"), all = TRUE)
coefs.exp.plt[is.na(enriched.mtblmcs),enriched.mtblmcs := 0]

coefs.exp.plt[, enriched.sum := enriched.host + enriched.mtblmcs+1]
coefs.exp.plt[, enriched := c("None", "Host","Metabolomics","Both")[enriched.sum]]
coefs.exp.plt[set == "HBMayo",set := "HBI/Mayo"]
coefs.exp.plt[,set := factor(set, levels = c("HBI/Mayo","Response","Remission"))]

### create the plot
plt.heatmap <- ggplot(coefs.exp.plt, aes(y=subsystem, x = set,  size = affected))+
    geom_point( alpha = 0.5,shape = 21, aes(fill = enriched))+
    labs(x = "",
         y = "",
         size ="Ratio affected
met. in subsys.",
         fill = "Enriched in") +
    guides(#size = "none", # remove size legend
           y.sec = "axis", y = "none" # move y-tick labels to right
           )+
    #scale_color_gradient(low = "gray90", high = "red", limits = c(0,NA))+
    scale_fill_brewer(palette = "Dark2")+
    scale_shape_manual(values = c(19,4))+
    facet_grid(~tissue)+
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
#plt.heatmap

plts <- list()
rel.heights <- c()
for(set in rev(coef.fls)){
    # read coefficients and create a data.frame for plotting
    coefs <- fread(file.path(ddr, set))
    coefs.exp <- merge(coefs, cluster[,.(met, rep.met)], by.x= "metabolite", by.y = "rep.met", suffix = c(".rep",".all"))
    coefs.exp.plt <- merge(coefs.exp[padj < 0.05,], met.meta[,.(colAbbr,Recon3D)], by.x = "metabolite", by.y = "colAbbr", all.x = TRUE)
    coefs.exp.plt <- merge(coefs.exp.plt, unique(met2subsystem[,.(subsystem,metabolite.woComp)]), by.x = "Recon3D", by.y= "metabolite.woComp")
    # order for subsystems
    orderCat <- coefs.exp.plt[,.(odr = median(Estimate)), by = subsystem][order(odr), subsystem]
    coefs.exp.plt[,subsystem:= factor(subsystem, levels = orderCat)]
    # add enriched subsystems
    coefs.exp.plt[subsystem %in% enrichments[set == set,unique(subsystem)], enriched := "Metabolomics"] 
    coefs.exp.plt[is.na(enriched),enriched := "None"]
    # sanity check, if any metabolites are still available
    if(nrow(coefs.exp.plt) >0){

        # adjust plot to data set
        set.title <- strsplit(set,
                                 split = "\\.")[[1]][3]
        set.title <- gsub("HBMayo","HBI/Mayo score", set.title)
        if(set.title != "HBI/Mayo score"){
            coefs.exp.plt[,Estimate := sign(Estimate)* sqrt(abs(Estimate))]
            coefs.exp.plt[,`Std. Error` := sign(`Std. Error`)* sqrt(abs(`Std. Error`))]
            xlab <- "sqrt(Estimate)"
        } else {
            xlab <- "Estimate"
        }
        rel.heights <- c(rel.heights, 2+length(unique(coefs.exp.plt[,subsystem])))

        # create the actual plot and store in a list
        plts[[set]] <- 
            ggplot(coefs.exp.plt, aes(x = Estimate, y = subsystem)) +
                geom_vline(xintercept = 0, linetype = 2, color = "red")+
                geom_boxplot(alpha = 0.5, aes(fill = enriched)) +
                scale_fill_manual(values = brewer.pal(name = "Dark2", n = 4)[3:4])+
                labs(x = xlab, y = "", title = set.title)+
                theme_bw() +
                theme(legend.position = "none") 
    }
}



plt.left  <- plot_grid(plotlist= plts, labels = "AUTO", ncol = 1, rel_heights = rel.heights)
plt <- plot_grid(plt.left, plt.heatmap, labels = c(NA,"D"), ncol = 2)

ggsave(plt, 
       file = file.path("..","figures_raw", "FigS1_metabolomicsDiseaseAssociationSubsystems.pdf"),
       height = 12,
       width = 12)
