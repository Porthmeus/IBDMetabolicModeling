# Porthmeus
# 06.06.23

require(data.table)
require(ggplot2)
require(clusterProfiler)
require(cowplot)
require(stringr)
require(fpc)
require(ggvenn)

# set some variables
dat.dir <- file.path("..","data", "Microbiome")
set.seed(682968)

phenotypes <- c(`HBI/Mayo score` = "HBMayo",
          Response = "Response",
          Remission = "Remission")

sim.sets <- c("MicrobiomeGS","Bacarena")


datasets <- data.table()
for(s.set in sim.sets){
    flux.cluster <- fread(file.path(dat.dir,s.set, "flux_cluster.csv"))
    for(i in 1:length(phenotypes)){
        # get the right file
        ph.type <- phenotypes[i]
        fl <- list.files(file.path(dat.dir,s.set), pattern = paste0("^flux", "\\.G?LMM\\.", ph.type, ".*", "_coefs.csv"))
        dat <- fread(file.path(dat.dir,s.set, fl))
        ##
        # get the significant reactions and expand the cluster
        dat <- dat[padj < 0.05,]
        dat <- merge(flux.cluster[,.(rep.rxn, rxn,sign.cor)], dat, by.x = "rep.rxn", by.y = "rxn")
        ##
        # add the information about simulation and phenotype
        dat <- cbind(sim.set = s.set, phenotype = names(phenotypes)[i], dat)
        ##
        # do some tinkering to make columns conform and put everything in one table
        if("rxn_id" %in% colnames(dat)){
            dat[,rxn := rxn_id]
            dat <- dat[,-"rxn_id"]
        }
        dat[,rxn := gsub("_[a-z]0$","",rxn)]
        colnames(dat)[grep("[t|z] value", colnames(dat))] <- "EffSize"
        colnames(dat)[grep("Pr(>|[t|z]|)",colnames(dat))] <- "p.value"
        # correct effect size and Estimates from the model
        dat[,Estimate := Estimate * sign.cor]
        dat[,EffSize := EffSize * sign.cor]
        if(ncol(datasets) >1){
            cols <- intersect(colnames(datasets), colnames(dat))
            datasets <- rbind(datasets[,..cols], dat[,..cols])
        } else {
            datasets <- rbind(datasets, dat)
        }

    }
}

dataset.flux <- datasets
# scale the estimates
dataset.flux[phenotype != "HBMayo", `Std. Error`:= sign(`Std. Error`)*sqrt(abs(`Std. Error`)), by = "phenotype"]
dataset.flux[phenotype != "HBMayo", Estimate := sign(Estimate)*sqrt(abs(Estimate)), by = "phenotype"]

# add reaction names
rxn2name <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_reactions_corrected.tsv")[,.(id, abbr_rxn = abbreviation, name_rxn= name)]
dataset.flux <- merge(dataset.flux, rxn2name, by.x = "rxn", by.y="id", all.x = TRUE)
dataset.flux[is.na(name_rxn),name_rxn := rxn]

# shorten them
dataset.flux[,name_rxn:= 
             gsub("UDP-GlcNAc:undecaprenyl-diph.-N-acetylmuramoyl-L-alanyl-gamma-D-glutamyl-meso-2,6-diaminopimeloyl-D-alanyl-D-alanine 4-beta-N-acetylglucosaminlytransfer.","ph.-n-ac.muramoyl-pentapeptide transfer.", fixed = TRUE,
         gsub("TECA4S", "Minor teichoic acid synthesis n30",
         gsub("TRANS-RXN-108F.cp","Transport of deoxyuridine",
         gsub("hydroxy","OH",
         gsub("Hydroxy","OH",
         gsub("6,7-Dimethyl-8-(1-D-ribityl)lumazine:6,7-dimethyl-8-(1-D-ribityl)lumazine 2,3-butanediyltransfer.","Riboflavin synth.", fixed = TRUE,
         gsub("2,5-Diamino-6-hydroxy-4-(5-ph.ribosylamino)-pyrimidine 2-aminohydrol.","5-amino-6-(5-ph.ribosylamino)uracil reduct.", fixed = TRUE,
         gsub("UDP-N-acetyl-D-mannosamine:N-acetyl-beta-D-glucosaminyldiph.undecaprenol beta-1,4-N-acetylmannosaminyltransfer.","N-acetylmannosaminyltransferase",
         gsub("TRANS-RXNBWI-115637.ce.maizeexp.L-ASPARTATE_L-ASPARTATE","Asparatate transport via proton symport",
         gsub("TRANS-RXNBWI-115637.ce.maizeexp.ASN_ASN","Asparagine transport via proton symport",
         gsub("TRANS-RXNAVI-26524.ce.brachyexp.SUCROSE_SUCROSE","Transport of Sucrose",
         gsub("acyl-carrier protein","ACP",
         gsub("acyl-carrier-protein","ACP",
         gsub("phospho","ph.",
         gsub("phospha","ph.",
         gsub("phosphate","ph.",
         gsub("ase",".",
         gsub(" $","",
         gsub("alpha","a", 
         gsub("^-","", 
         gsub("substituted","sub.",
         gsub("synthesis","synth.",
         gsub("lipoteichoic acid","LTA",
         gsub("N-acetylglucosamine","GlcNAc",
         gsub("N-acetyl-D-glucosamine","GlcNAc",
         gsub("(.*)\\(.*\\)$","\\1",
         gsub("UDP-2,3-bis[(3R)-3-hydroxymyristoyl]-alpha-D-glucosamine 2,3-bis[(3R)-3-hydroxymyristoyl]-beta-D-glucosaminyl 1-phosphate phosphohydrolase", "udp-2,3-diacylglucosamine hydrolase", fixed = TRUE,
         gsub("UDP-2,3-bis(3-hydroxytetradecanoyl)glucosamine:2,3-bis-(3-hydroxytetradecanoyl)-alpha-D-glucosaminyl-1-phosphate 2,3-bis(3-hydroxytetradecanoyl)-glucosaminyltransferase", "lipid-a-disaccharide synthase", fixed = TRUE,
         gsub("2,5-diamino-6-(5-phospho-D-ribitylamino)pyrimidin-4(3H)-one:NAD+ 1'-oxidoreductase","2,5-diamino-6-ribosylamino-4(3H)-pyrimidinone 5'-phosphate reductase", fixed = TRUE,
         gsub("UDP-N-acetylmuramoyl-L-alanyl-gamma-D-glutamyl-meso-2,6-diaminopimeloyl-D-alanyl-D-alanine:undecaprenyl-phosphate phospho-N-acetylmuramoyl-pentapeptide-transferase","phospho-n-acetylmuramoyl-pentapeptide transferase", fixed = TRUE,
              name_rxn))))))))))))))))))))))))))))))]
unique(dataset.flux[,name_rxn])

# create a plot per reaction
rxn.plts <- list()
rel.height <- c()
for(phntp in dataset.flux[,unique(phenotype)]){
    dat.plt <- dataset.flux[phenotype == phntp,]
    ## sort the names
    dat.plt[, name_rxn := factor(name_rxn, levels = unique(name_rxn[order(Estimate)]))]
    rel.height <- c(rel.height, length(dat.plt[,unique(name_rxn)]))
    inner.pl <- ggplot(dat.plt, aes(x = Estimate, y =name_rxn ,color = sim.set)) +
        geom_vline(xintercept = 0, linetype = 2, color = "red") +
        geom_point() +
        scale_color_brewer(palette = "Paired")+
        geom_errorbar(aes(xmin = Estimate - `Std. Error`, xmax = Estimate + `Std. Error`))+
        theme_bw() +
    #    theme(axis.text.y = element_text(angle =30, vjust = 1))+
        labs(x = "Estimate/sqrt(Estimate)", y = "Metabolite", title = phntp)
    rxn.plts[[phntp]] <- inner.pl
}

rxn.plt <- plot_grid(plotlist = rxn.plts, ncol = 1, rel_heights = rel.height)
ggsave(rxn.plt,
       limitsize = FALSE,
       height = 60,
       width = 12,
       file = file.path("..","figures_raw","FigSXX_internalRxns.Microbiome.Estimates.pdf"))

# create venn diagram for the overlap between both modeling approaches
venn.plts <- list()
venn.dats <- list()
for(phntp in dataset.flux[,unique(phenotype)]){
    venn.dat <- list()
    for(sim in dataset.flux[,unique(sim.set)]){
        venn.dat[[sim]] <- dataset.flux[sim.set == sim & phenotype==phntp, unique(rxn)]
    }
    venn.dats[[phntp]] <- venn.dat
    venn.plt <- ggvenn(venn.dats[[phntp]]) +
        scale_fill_brewer(palette = "Paired") +
        labs(title =phntp)
    venn.plts[[phntp]] <- venn.plt
}

venn.plt <- plot_grid(plotlist = venn.plts, nrow=1)
ggsave(venn.plt,
       file = file.path("..","figures_raw", "FigSXX_Microbiome.internalRxnsVenn.pdf"),
       width = 10, height = 3)


## make pathway enrichment
pwytb <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/meta_rea_pwy-gapseq.tbl")
pwytb[,pathway := gsub("\\|","",pathway)]
rxn2pwy <- pwytb[,.(rxn = unlist(strsplit(gapseq, split = " "))), by = .(pathway, pathway.name)]
pwy2name <- unique(rxn2pwy[,.(pathway, pathway.name)])
rxn2pwy <- unique(rxn2pwy[,.(rxn, pathway)])

# add the names again
dataset.flux <- merge(dataset.flux, rxn2pwy, by = "rxn")

gseas <- list()
for(ph.type in names(phenotypes)){
    dat <- dataset.flux[phenotype == ph.type,]
    vec.tab <- dat[,.(val = sum(EffSize)), by = rxn]
    vec <- vec.tab[,val]
    names(vec) <- vec.tab[,rxn]
    vec <- sort(vec, decreasing = TRUE)
    gsea <- GSEA(geneList = vec,
                 TERM2GENE = rxn2pwy[,.(pathway, rxn)],
                 TERM2NAME = pwy2name)
    gseas[ph.type] <- gsea
}

hygeos <- list()
for(ph.type in names(phenotypes)){
    dat <- dataset.flux[phenotype == ph.type,]
    vec <- unique(dat[padj <0.05,rxn])
    hygeo <- enricher(gene = vec,
                 TERM2GENE = rxn2pwy[,.(pathway, rxn)],
                 TERM2NAME = pwy2name)
    hygeos[[ph.type]] <- hygeo
}

# combine the enrichment results
dat.nrch <- data.table()
for(n in names(phenotypes)){
    dat1 <- gseas[[n]]@result[,c("ID","Description","pvalue","setSize")]
    colnames(dat1)[4] <- "Count"
    dat1 <- dat1[dat1$Count > 2,]
    if(nrow(dat1) >0){
        dat1 <- cbind(dat1, method = "gsea")
    }
    dat2 <- hygeos[[n]]@result[,c("ID","Description","pvalue","Count")]
    dat2 <- dat2[dat2$Count > 2,]
    if(nrow(dat2) >0){
        dat2 <- cbind(dat2, method = "hygeo")
    }
    dat <- data.table(cbind(rbind(dat1,dat2), phenotype = n))
    dat.nrch <- rbind(dat.nrch, dat)
}

dat.nrch[,padj := p.adjust(pvalue, method = "BH")]
dat.nrch <- dat.nrch[padj  <0.05,]
# sort the pathways for those relevant for bacteria
dat.nrch[,flag1 := grepl("bacteria", Description)]
dat.nrch[,flag2 := grepl("fungi|archea|plants|eukaryotes|mitochondrial|engineered", Description)]
dat.nrch <- dat.nrch[flag1==flag2,]
subs <- list("\\(S\\)-" = "",
             hydrogen = "H+",
             biosynthesis = "synth.",
             fermentation = "ferm.",
             deoxyribonucleotide = "dN",
             ribonucleotide = "N",
             tranformation = "transform.",
             anaerobic = "anaer.",
             aerobic = "aer.")
for(i in 1:length(subs)){
    old = names(subs)[i]
    nw = subs[[i]]
    dat.nrch[,Description:= gsub(old,nw, Description)]
}

# plot the results
## expand the pathway ID to corresponding reactions
rxn2nrch <- merge(rxn2pwy[,.(rxn, pathway)], dat.nrch, by.x = "pathway", by.y = "ID")
plt.dat.nrch <- merge(dataset.flux[,.(rxn, sim.set, phenotype, Estimate, padj, pathway = pathway.y)],rxn2nrch,
                      by = c("rxn", "phenotype","pathway"),
                      suffix = c(".flux",".pathway"))
#plt.dat.nrch <- unique(plt.dat.nrch)
# order the pathways and remove duplicates
plt.dat.ord <- plt.dat.nrch[,.(od = median(Estimate),mn = mean(Estimate),sel = paste0(phenotype,median(Estimate),length(rxn))), by = .(pathway, phenotype)]
#plt.dat.ord <- plt.dat.ord[!duplicated(sel),]
plt.dat.nrch <- merge(plt.dat.ord, plt.dat.nrch)
#plt.dat.nrch[,sel := paste0(sel,pathway,rxn)]
#plt.dat.nrch <-plt.dat.nrch[!duplicated(sel),]

pwy.plts <- list()
rel.height <- c()
dat.final <- data.table()
for(n in names(phenotypes)){
    dat <- plt.dat.nrch[phenotype == n,]
    # group the subsystems
    a <- matrix(0, ncol = length(dat[,unique(pathway)]), nrow = length(dat[,unique(rxn)]),
                dimnames = list(dat[,unique(rxn)], dat[,unique(pathway)]))
    for(pwy in dat[,unique(pathway)]){
        a[dat[pathway == pwy,rxn],pwy] <- 1
    }
    a.dist <- dist(t(a), method ="binary")
    cluster <- dbscan(a.dist, method = "dist", MinPts = 1, eps = 0.9)
    tbl.cluster <- data.table( cluster = cluster[["cluster"]],
                              eps = cluster[["eps"]],
                              MinPts = cluster[["MinPts"]],
                              isseed = cluster[["isseed"]],
                              pathway = colnames(a))
    dat <- merge(dat, tbl.cluster, by = "pathway")
    #dat[,Description := factor(Description, levels = unique(Description[order(mn)]))]
    dat[,group := paste0("Group ",cluster)] 
    dat.plt <- unique(dat[,.(Estimate, group, sim.set)])
    # add the manual found group names
    group2name <- file.path("..","tables","TabSX_Goups2ManualPathways.csv")
    if(file.exists(group2name)){
        group2name <- fread(group2name)
        dat.plt <- merge(dat.plt, group2name[phenotype == n,], by ="group")
        dat.plt[, cluster := group]
        dat.plt[, group := Description2]
    }
    # remove glucosylglycerol synth. as it is only known in Cyanos
    dat.plt <- dat.plt[group != "Glucosylglycerol synth.",]
    dat.ord <- dat.plt[,.(od = median(Estimate),mn = mean(Estimate)), by = group]
    dat.plt <- merge(dat.plt, dat.ord, by =  "group")
    dat.plt[,group := factor(group, levels = unique(group[order(od)]))]
    #dat[,cluster := factor(paste0("Group",cluster), levels = unique(paste0("Group",cluster)[order(mn)]))]
    if(n != "HBI/Mayo score"){
        xlab = "sqrt(Estimate)"
    } else {
        xlab = "Estimate"
    }
    # censor some outlier
    dat.plt[abs(Estimate) > 5, Estimate := sign(Estimate)*5]
    plt <- ggplot(dat.plt, aes(x=Estimate, y =group)) +
        geom_vline(xintercept = 0, linetype = 2, color = "red" )+
        geom_jitter(width = 0, height = 0.25, alpha = 1, aes(color = sim.set))+
        geom_point(shape = 24,alpha = 1, fill = "#33a02c",size = 2,data=dat.ord, aes( x =od))+
        labs(title = n, y = "Pathway", x =xlab )+
        scale_color_brewer(palette = "Paired")+
       # scale_y_discrete(labels = function(x) str_wrap(x, width = 50))+
        theme_bw() 
#        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    pwy.plts[[n]] <- plt
    rel.height <- c(rel.height, length(dat[,unique(group)]))
    dat.final <- rbind(dat.final, dat)
}
x <- 6.5
plt.all<- cowplot::plot_grid(plotlist= pwy.plts, ncol = 1, rel_heights = (rel.height/sum(rel.height))+(1/x), labels = c("D","E","F"))
ggsave(plt.all,
       file = file.path("..","figures_raw","Fig1_EnrichedPathwaysMicrobeRxn.pdf"),
       height = x,
       width = 6)
write.csv(dat.final, row.names = FALSE,
          file = file.path("..","temp","EnrichedMicPathways.csv"))
