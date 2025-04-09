# Porthmeus
# 12.03.25

# libraries
require(clusterProfiler)
require(data.table)
require(ggplot2)

# directories
dat.dir.net <- file.path("..","data","Networks")
dat.dir.host <- file.path("..","data","Host")
res.dir <- file.path("..","figures_raw")

# load data
subsystems <- fread(file.path(dat.dir.host, "subsystems.csv"))
plts <- list()
for(st in c("Response","HBMayo")){
    st_title <- gsub("HBMayo","HBI/Mayo index",st)
    for(tis in c("biopsy","blood")){
        shortest_paths <- fread(file.path(dat.dir.net, paste("shortPath_005quant",tis,st,"csv",sep = ".")))
        enr <- enricher(gene = shortest_paths[,unique(rxn)],
                        TERM2GENE = subsystems)
        plts[[paste(st,tis,sep = ".")]] <- dotplot(enr) + ggtitle(tools::toTitleCase(tis))
    }
}

response.plts <- cowplot::plot_grid(plotlist = plts[grep("Response",names(plts))], labels = "AUTO")
ggsave(response.plts, file =file.path(res.dir,"FigSX_enrichmentShortestPath.Response.pdf"), height = 4, width =12)
HBMayo.plts <- cowplot::plot_grid(plotlist = plts[grep("HBMayo",names(plts))], labels = "AUTO")
ggsave(HBMayo.plts, file =file.path(res.dir,"FigS6_enrichmentShortestPath.HBMayo.pdf"), height = 4.5, width =12)
#cowplot::plot_grid(plotlist = rev(plts),ncol = 2)

