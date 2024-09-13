# Porthmeus
# 07.09.23

require(data.table)
require(limma)
require(cowplot)
require(RColorBrewer)
require(ggvenn)

rsrc.dir <- file.path("..","tables")
res.dir <- file.path("..","figures_raw")

# load data
coefs <- fread(file.path(rsrc.dir, "TabS1_StatisticsHost.csv"))

# create the venn mat
types <- coefs[,unique(Data.Type)]
cndts <- coefs[,unique(Condition)]
rxns <- coefs[,unique(rxn)]
venns <- list()
for(tis in c("blood","biopsy")){
    for(cndt in cndts){
        venn <- list()
        for(typ in types){
            sig.rxns <- coefs[padj < 0.05 & Data.Type == typ & Tissue == tis & Condition == cndt, unique(rxn)]
       #     mat[sig.rxns,typ] <- 1 
            venn[[typ]] <- sig.rxns
        }
        cndt <- gsub("HBMayo", "HBI/Mayo", cndt)
        venns[[paste(cndt,tools::toTitleCase(tis))]] <- venn
    }
}

plts <- list()
for(tis in sort(names(venns))){
    plt <- ggvenn(venns[[tis]])+
        scale_fill_brewer(palette = "Dark2") +
        ggtitle(tis)
    plt
    plts[[tis]] <- plt
}

plt <- plot_grid(plotlist = plts, labels = "AUTO", ncol = 2)

ggsave(plt,
       file = file.path(res.dir, "FigSX_VennSignificantRxnsHost.pdf"),
       width = 12,
       height = 18)
