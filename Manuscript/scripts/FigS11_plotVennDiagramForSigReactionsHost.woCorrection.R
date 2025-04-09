# Porthmeus
# 07.09.23

require(data.table)
require(limma)
require(cowplot)
require(RColorBrewer)
require(ggvenn)

rsrc.dir <- file.path("..","data", "Host")
res.dir <- file.path("..","figures_raw")

# load data

# get the network files and make an hgeo enrichment
tissues <- c("biopsy","blood")
sets <- c("Response","Remission")
prfx <- c("rxnExpr","PA","FVA")
coefs.all <- data.table()
for(tiss in tissues){
    for(set in sets){
        sig.rxns <- c()
        for(lr in prfx){
            # load clusters to expand them
            clus <- fread(file.path(rsrc.dir, tiss, paste0(lr, "_DBSCAN_Cluster.csv")))
            dat <- fread(file.path(rsrc.dir,tiss, paste0(lr, ".LMM.", set,"_coefs.csv")))
            if(lr == "FVA"){
                # merge cluster and coef dat 
                dat <- merge(clus[,.(rep.rxn,rxn, sign.cor.range, sign.cor.center)],dat, by.x = "rep.rxn", by.y = "rxn")
                dat[coef == "center", Estimate := Estimate * sign.cor.center]
                dat[coef == "range", Estimate := Estimate * sign.cor.range]
                coefs.all<- rbind(coefs.all, dat[padj < 0.05,.(layer=lr, set = set, tissue = tiss, rxn, Estimate)])
            } else {
                dat <- merge(clus[,.(rep.rxn,rxn, sign.cor)],dat, by.x = "rep.rxn", by.y = "rxn")
                dat[, Estimate := Estimate * sign.cor]
                coefs.all <- rbind(coefs.all, dat[padj < 0.05,.(layer=lr, set = set, tissue = tiss, rxn, Estimate)])
            }
        }
    }
}


# create the venn mat
types <- coefs.all[,unique(layer)]
cndts <- coefs.all[,unique(set)]
rxns <- coefs.all[,unique(rxn)]
venns <- list()
for(tis in c("blood","biopsy")){
    for(cndt in cndts){
        venn <- list()
        for(typ in types){
            sig.rxns <- coefs.all[layer == typ & tissue == tis & set == cndt, unique(rxn)]
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
       file = file.path(res.dir, "FigS11_VennSignificantRxnsHost.woCorrection.pdf"),
       width = 12,
       height = 12)
