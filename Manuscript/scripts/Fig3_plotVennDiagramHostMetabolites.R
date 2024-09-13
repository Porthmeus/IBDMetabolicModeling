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
rsrc.dir <- file.path("..","data","Host")
net.dir <- file.path("..","data","Networks")

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

# create a venn data frame
venn.dat <- dcast(hygeos[sig != FALSE], ID+set~tissue, value.var = "sig", fill = FALSE)

venn.plts <- list()
for(st in sets){
    venn.dat.plt <- venn.dat[set == st,]
    if(st == "HBMayo"){
        ttl <- "HBI/Mayo score"
    } else {
        ttl <- st
    }
    venn.plt <- ggplot(venn.dat.plt, aes(A=biopsy, B=blood)) +
        geom_venn(fill_color = brewer.pal(n = 3, "Set1")[c(2,1)]) +
        theme_void() +
        coord_fixed()+
        ggtitle(ttl)+
        theme(strip.text.x = element_text(size = 18))
    venn.plts[[st]] <- venn.plt
}
venn.plt <- plot_grid(plotlist = venn.plts, nrow = 1)


# load the sig. reactions and expand the cluster
prfx <- c("PA","rxnExpr","FVA")
coefs.all <- data.table()
for(tissue in c("biopsy","blood")){
    #tissue <- "biopsy"
    
    # load the clustering data
    files <- sapply(prfx, function(x) list.files(file.path(rsrc.dir, tissue), pattern = paste0(x,".*Cluster.csv")))
    cluster <- lapply(names(files), function(x) cbind(data = x, fread(file.path(rsrc.dir, tissue, files[[x]]))))
    names(cluster) <- prfx

    for(set in c("HBMayo","Remission","Response")){
        #set <- "HBMayo"
        files <- sapply(prfx, function(x) list.files(file.path(rsrc.dir, tissue), pattern = paste0(x,"\\..*\\.",set,"\\.")))
        coefs <- lapply(names(files), function(x) cbind(data = x, fread(file.path(rsrc.dir, tissue, files[[x]]))))
        names(coefs) <- prfx

        # expand the clusters
        for(prf in prfx){
            coefs[[prf]] <- merge(cluster[[prf]], coefs[[prf]], by.y = "rxn", by.x ="rep.rxn")
            if(prf != "FVA"){
                coefs[[prf]][, c("Estimate") := list(sign.cor * Estimate)]
            } else {
                coefs[[prf]][coef == "range", c("Estimate") := list(sign.cor.range * Estimate)]
                coefs[[prf]][coef == "center", c("Estimate") := list(sign.cor.center * Estimate)]
            }
            coefs.all <- rbind(coefs.all, cbind(tissue = tissue, set = set, coefs[[prf]][,.(rxn,Estimate,padj, dat = data.y)]))
        }
    }
}

# get the list of metabolites in the rections
met2rxn2ss <- fread(file.path(rsrc.dir, "Mets2Subsystems.csv"))

### create a barplot 
met2rxn2ss[,met.id.s := gsub("\\[[a-z]\\]","",metabolite)]
coefs.all.sig <- coefs.all[padj <0.05,]
coefs.all.sig <- coefs.all.sig[,.(Estimate = mean(Estimate)), by =.(tissue, set, rxn)]
# remove transport and exchange reactions
coefs.all.sig <- coefs.all.sig[rxn %in% network_all[,unique(rxn.id)],]
mets.coefs <- merge(coefs.all.sig,met2rxn2ss[,.(met.id.s, reaction,  met.name)], by.x = "rxn", by.y = "reaction", allow.cartesian=TRUE)
mets.coefs <- merge(mets.coefs, hygeos[sig != FALSE,.(tissue, set, met = ID, pvalue, p.adj)], by.x = c("tissue","set","met.id.s"), by.y = c("tissue","set","met"))
# count the up and down reactions per metabolite
mets.coefs[, UpDown := sign(Estimate)]
mets.cnt <- mets.coefs[,.(up = sum(UpDown == 1), down = sum(UpDown == -1), total = .N), by = .(tissue, set, met.name,met.id.s)]
# create NAs for missing values
mets.cnt[, c("rel.up","rel.down","rel.diff") := list(up/total, down/total, (up-down)/total)]
#mets.cnt.l <- dcast(mets.cnt, tissue+set~met.name, value.var = "total", fill = NA)
#mets.cnt.m <- melt(mets.cnt.l, id.vars = c("tissue","set"), variable.name = "met.name", value.name = "total")
#met.cnt <- merge(mets.cnt, mets.cnt.m[is.na(total),], by = c("tissue","set","met.name"),all = TRUE, suffix = c("",".y"))
#met.cnt[is.na(total), c("rel.up","rel.down","rel.diff","total") := list(0,0,0,0)]
#met.cnt[,total.y := NULL]
# select only the inner circles of the venns
sel <- venn.dat[biopsy == TRUE & blood == TRUE,.(met=ID,set)]
mets.cnt.sel <- merge(sel, mets.cnt, by.x = c("met","set"),by.y=c("met.id.s","set"))

# some reformatting for the plotting
mets.cnt.sel[,tissue := factor(tissue, levels = c("blood","biopsy"))]
mets.cnt.sel[,set := factor(gsub("HBMayo","HBI/Mayo score", set), levels = c("HBI/Mayo score","Response","Remission"))]
mets.cnt.sel[,met.name := gsub("Coenzyme A", "CoA",met.name)]
mets.cnt.sel[,met.name := gsub("Nicotinamide Adenine Dinucleotide Phosphate", "NADP",met.name)]
mets.cnt.sel[,met.name := gsub("Nicotinamide Adenine Dinucleotide", "NAD+",met.name)]
mets.cnt.sel[,met.name := gsub("Flavin Adenine Dinucleotide", "FAD",met.name)]
mets.cnt.sel[,met.name := gsub("Adenosine Monophosphate", "AMP",met.name)]
mets.cnt.sel[,met.name := gsub("Adenosine Diphosphate", "ADP",met.name)]
mets.cnt.sel[,met.name := gsub("Adenosine Triphosphate", "ATP",met.name)]

nudge <- rep(c(0.25,-0.25),length(unique(mets.cnt.sel[order(set, met.name, tissue),paste0( met.name, set, tissue)]))/2)
bar.plt <- ggplot(mets.cnt.sel, aes(x=met.name, y=rel.up, fill = tissue )) +
    geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
    geom_hline(yintercept = -0.5, color = "black", linetype = "dashed")+
    geom_hline(yintercept = 0, color = "black")+
    geom_bar(position =position_dodge(preserve = "single"), stat = "identity", alpha = 0.7, color = "black") +
    geom_bar(position =position_dodge(preserve = "single"), stat = "identity", alpha = 0.7, color = "black", aes(y=-rel.down)) +
    geom_text(aes(label = total,
                  hjust = (-sign(rel.diff)+1)/2,
              y = sign(rel.diff)*0.01),
            nudge_x = nudge, angle = 90, size = 3)+
#    geom_bar(position = "dodge", stat = "identity", alpha = 0.7, aes(y = -1*rel.down))+
    theme_bw() +
    scale_fill_brewer(palette = "Set1")+
    facet_grid(.~set, scale = "free", space = "free_x")+
    labs( y= "Rel. no. of pos./neg. reactions", x = element_blank(), fill = "Tissue")+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
bar.plt

plt <- plot_grid(venn.plt, bar.plt, ncol = 1)
ggsave(plt, file = file.path("..","figures_raw","FigSX_VennBarMetabolitesHost.pdf"),
       width = 10, height = 7)

# create a plot to seperate between substrate and product in reactions
nwk.melt <- unique(melt(network_all[,.(rxn.id, from.id, to.id, arrows, set, tissue, EffSize)],
                 id.vars = c("rxn.id","arrows","set","tissue","EffSize")))

nwk.cnt <- nwk.melt[, .(substrate.up = sum(EffSize > 0 & variable == "from.id" & arrows == "to"),
             substrate.down = -sum(EffSize < 0 & variable == "from.id" & arrows == "to"),
             product.up = sum(EffSize > 0 & variable == "to.id" & arrows == "to"),
             product.down = -sum(EffSize < 0 & variable == "to.id" & arrows == "to"),
             both.up = sum(EffSize > 0 &  arrows == "to;from"),
             both.down = -sum(EffSize < 0  & arrows == "to;from"),
             total = .N),
        by = .(value, set, tissue)]
mets.nwk.cnt <- merge(mets.cnt, nwk.cnt, by.x = c("tissue","set","met.id.s", "total"), by.y = c("tissue","set","value", "total"))
met.nwk.cnt.melt <- melt(mets.nwk.cnt[,.(tissue, set, met.id.s, met.name,total, substrate.up, substrate.down, product.up, product.down, both.up, both.down)], id.vars = c("tissue", "set", "met.id.s", "met.name","total"), value.name = "count")
met.nwk.cnt.melt[,rel.count := count/total]
met.nwk.cnt.melt[,Metabolite := gsub("\\.up|\\.down","", variable)]

met.nwk.cnt.melt[,tissue := factor(tissue, levels = c("blood","biopsy"))]
met.nwk.cnt.melt[,set := factor(gsub("HBMayo","HBI/Mayo score", set), levels = c("HBI/Mayo score","Response","Remission"))]
met.nwk.cnt.melt[,met.name := gsub("Coenzyme A", "CoA",met.name)]
met.nwk.cnt.melt[,met.name := gsub("Nicotinamide Adenine Dinucleotide Phosphate", "NADP",met.name)]
met.nwk.cnt.melt[,met.name := gsub("Nicotinamide Adenine Dinucleotide", "NAD",met.name)]
met.nwk.cnt.melt[,met.name := gsub("Flavin Adenine Dinucleotide", "FAD",met.name)]
met.nwk.cnt.melt[,met.name := gsub("Adenosine Monophosphate", "AMP",met.name)]
met.nwk.cnt.melt[,met.name := gsub("Adenosine Diphosphate", "ADP",met.name)]
met.nwk.cnt.melt[,met.name := gsub("Adenosine Triphosphate", "ATP",met.name)]

bar2.plt <- ggplot(met.nwk.cnt.melt, aes(y=rel.count, x = met.name))+
    geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
    geom_hline(yintercept = -0.5, color = "black", linetype = "dashed")+
    geom_hline(yintercept = 0, color = "black")+
    geom_bar(stat = "identity", aes(fill = Metabolite), alpha = 0.7, color = "black") + 
    facet_wrap(set~tissue, scale = "free", ncol = 2)+
    theme_bw()+
    scale_fill_brewer(palette="Paired")+
    labs( y= "Rel. no. of pos./neg. reactions", x = element_blank())+
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
bar2.plt
# create a hbmayo_only version

met.nwk.cnt.melt[,tissue := tools::toTitleCase(as.character(tissue))]
met.nwk.cnt.melt[,total.label := abs(total)]
met.nwk.cnt.melt[duplicated(total.label) & duplicated(paste(met.name, tissue, set)), total.label := NA]
# get the direction on which side of the 0 the number should be plotted
updown.label <- met.nwk.cnt.melt[,.(updown = sign(sum(rel.count))),by = .(tissue, set, met.id.s)]
met.nwk.cnt.melt <- merge(met.nwk.cnt.melt, updown.label, by = c("tissue","set","met.id.s"))
met.nwk.cnt.melt[updown == 0, updown :=1]

# remove very long metabolite names to exchange them manualy if needed
met.nwk.cnt.melt[nchar(met.name) >30,met.name := met.id.s]
for(st in unique(met.nwk.cnt.melt[,set])){
    bar2.plt.hbmayo <- ggplot(met.nwk.cnt.melt[set == st,], aes(y=rel.count, x = met.name))+
        geom_hline(yintercept = 0.5, color = "black", linetype = "dashed")+
        geom_hline(yintercept = -0.5, color = "black", linetype = "dashed")+
        geom_hline(yintercept = 0, color = "black")+
        geom_bar(stat = "identity", aes(fill = Metabolite), alpha = 0.7, color = "black") + 
        facet_grid(tissue~., scale = "free")+#, ncol = 2)+
        theme_bw()+
        scale_fill_brewer(palette="Paired")+
        geom_text(aes(label = total.label, y = updown*0.2),
                  angle = 90, size = 3)+
        labs( y= "Rel. no. of pos./neg. reactions", x = element_blank())+
        theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
    bar2.plt.hbmayo
    ggsave(bar2.plt.hbmayo, file = paste0("../figures_raw/Fig1_MetabolitesBarplot.",gsub("/| ","",st)  ,".pdf"),
           width = 8, height = 5)
}


##################################################
#
#### create the venn diagram
## create a true/false vector for significance
#met.stats[, sig:= FALSE]
#met.stats[p.adj <0.025, sig:= TRUE] # its a one sided test, thus reduced alpha values
## recast into plotable version
#venn.dat <- dcast(met.stats[sig != FALSE,], met+set~tissue, value.var = "sig", fill = FALSE)
#
#venn.dat[,set2 := factor(gsub("HBMayo","HBI/Mayo score",set), levels = c("HBI/Mayo score", "Response","Remission"))]
#
#venn.plt <- ggplot(venn.dat, aes(A=biopsy, B=blood)) +
#    geom_venn(fill_color = brewer.pal(n = 3, "Set1")[c(2,1)]) +
#    theme_void() +
#    coord_fixed()+
#    facet_wrap(~set2)+
#    theme(strip.text.x = element_text(size = 18))
##venn.plt
#
#
#
