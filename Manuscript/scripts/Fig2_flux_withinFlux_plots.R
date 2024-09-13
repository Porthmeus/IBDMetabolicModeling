# Porthmeus
# 17.05.23

require(data.table)
require(ggplot2)
require(clusterProfiler)

dat.dir <- file.path("..","results")

phenotypes <- c(`HBI/Mayo score` = "HBMayo",
          Response = "Response",
          Remission = "Remission")

flux.sets <- c("flux","withinflux")

# read datasets
datasets <- data.table()
for(f.set in flux.sets){
    for(i in 1:length(phenotypes)){
        ph.type <- phenotypes[i]
        fl <- list.files(dat.dir, pattern = paste0("^",f.set, ".*", ph.type, ".*", "_coefs.csv"))
        dat <- fread(file.path(dat.dir, fl))
        dat <- dat[,-"df"]
        colnames(dat)[6:7] <- c("EffSize","p.value")
        dat <- cbind(flux.set = f.set, phenotype = names(phenotypes)[i], dat)
        datasets <- rbind(datasets, dat)
    }
}
datasets[,phenotype := factor(phenotype, levels = names(phenotypes))]

# scale the estimates
datasets.orig <- datasets
datasets <- datasets.orig[padj <0.05,]
datasets[flux.set == "withinflux" & phenotype != "HBMayo", `Std. Error`:= sign(`Std. Error`)*sqrt(abs(`Std. Error`)), by = "phenotype"]
datasets[flux.set == "withinflux" & phenotype != "HBMayo", Estimate := sign(Estimate)*sqrt(abs(Estimate)), by = "phenotype"]
datasets[, name := factor(name, levels = name[order(Estimate)])]

# plot the within metabolite exchange of the microbiome
datasets[,phenotype := factor(phenotype, levels = names(phenotypes))]
within.pl <- ggplot(datasets[flux.set == "withinflux",], aes(x = Estimate, y =name )) +
    geom_vline(xintercept = 0, linetype = 2, color = "red") +
    geom_point() +
    geom_errorbar(aes(xmin = Estimate - `Std. Error`, xmax = Estimate + `Std. Error`))+
    facet_wrap(~phenotype, scale = "free", ncol = 1) +
    theme_bw() +
    labs(x = "Estimate/sqrt(Estimate)", y = "Metabolite", title = "Exchange within Microbiome")
ggsave(within.pl, file = "withinFlux.pdf", height = 6, width = 4)


# expand the cluster
flux.cluster <- fread(file.path(dat.dir, "flux_cluster.csv"))

dataset.flux <- datasets.orig[flux.set == "flux" & padj < 0.05,]
dataset.flux <- merge(flux.cluster[,.(rxn, rep.rxn)], dataset.flux, by.x = "rep.rxn", by.y = "rxn", allow.cartesian=TRUE, all.y = TRUE)
# scale the estimates
dataset.flux[phenotype != "HBMayo", `Std. Error`:= sign(`Std. Error`)*sqrt(abs(`Std. Error`)), by = "phenotype"]
dataset.flux[phenotype != "HBMayo", Estimate := sign(Estimate)*sqrt(abs(Estimate)), by = "phenotype"]
dataset.flux[is.na(name_rxn), name_rxn := rxn]

# add the names again and shorten them
rxn2name <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_reactions_corrected.tsv")[,.(id, abbr_rxn = abbreviation, name_rxn= name)]
dataset.flux<- dataset.flux[,-"name_rxn"]
dataset.flux <- merge(dataset.flux, rxn2name, by.x = "rxn", by.y="id", all.x = TRUE)

dataset.flux[,name_rxn:= 
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
         gsub("(.*)\\(.*\\)$","\\1",name_rxn)))))))))))))))))))))))))]
dataset.flux[, name_rxn := factor(name_rxn, levels = unique(name_rxn[order(Estimate)]))]
dataset.flux[, abbr_rxn := factor(abbr_rxn, levels = unique(abbr_rxn[order(Estimate)]))]
dataset.flux[,unique(name_rxn)]
# plot the estimates of the reactions of the internal fluxes
dataset.flux[, name_rxn := factor(name_rxn, levels = unique(name_rxn[order(Estimate)]))]
dataset.flux[,phenotype := factor(phenotype, levels = names(phenotypes))]
inner.pl <- ggplot(dataset.flux, aes(x = Estimate, y =name_rxn )) +
    geom_vline(xintercept = 0, linetype = 2, color = "red") +
    geom_point() +
    geom_errorbar(aes(xmin = Estimate - `Std. Error`, xmax = Estimate + `Std. Error`))+
    facet_wrap(~phenotype, scale = "free") +
    theme_bw() +
#    theme(axis.text.y = element_text(angle =30, vjust = 1))+
    labs(x = "Estimate/sqrt(Estimate)", y = "Metabolite", title = "Flux of micr. reactions")
ggsave(inner.pl, file = "reactionFlux.pdf", height = 14, width = 16)


# create a GSEA for the reactions
rxn2pwy <- datasets.orig[flux.set == "flux",
                                 .(pathway = unlist(strsplit(gsub("\\|","",pathway), split = ";"))
                                                     ),
                               by = rxn]
rxn2pwy <- rxn2pwy[!is.na(pathway),]
# expand reactions
rxn2pwy <- merge(flux.cluster[,.(rxn, rep.rxn)], rxn2pwy, by.x = "rep.rxn", by.y = "rxn", allow.cartesian=TRUE)
pwy2name <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/meta_pwy.tbl")[,.(id = gsub("\\|","",id),name)]

rxn2pwy<- rxn2pwy[(rxn2pwy[,pathway] %in% pwy2name[,id]),]

gseas <- list()
for(ph.type in names(phenotypes)){
    dat <- dataset.flux[phenotype == ph.type,]
    vec <- dat[,EffSize]
    names(vec) <- dat[,rxn]
    vec <- sort(vec, decreasing = TRUE)
    gsea <- GSEA(geneList = vec,
                 TERM2GENE = rxn2pwy[,.(pathway, rxn)],
                 TERM2NAME = pwy2name)
    gseas[ph.type] <- gsea
}

hygeos <- list()
for(ph.type in names(phenotypes)){
    dat <- dataset.flux[phenotype == ph.type,]
    vec <- dat[padj <0.05,rxn]
    hygeo <- enricher(gene = vec,
                 TERM2GENE = rxn2pwy[,.(pathway, rxn)],
                 TERM2NAME = pwy2name)
    hygeos[[ph.type]] <- hygeo
}

# combine the enrichment results
dat.nrch <- data.table()
for(n in names(phenotypes)){
    dat1 <- gseas[[n]]@result[,c("ID","Description","pvalue")]
    if(nrow(dat1) >0){
        dat1 <- cbind(dat1, method = "gsea")
    }
    dat2 <- hygeos[[n]]@result[,c("ID","Description","pvalue")]
    if(nrow(dat2) >0){
        dat2 <- cbind(dat2, method = "hygeo")
    }
    dat <- data.table(cbind(rbind(dat1,dat2), phenotype = n))
    dat.nrch <- rbind(dat.nrch, dat)
}
dat.nrch <- dat.nrch[pvalue <0.05,]

# plot the results
## expand the pathway ID to corresponding reactions
rxn2nrch <- merge(rxn2pwy[,.(rxn, pathway)], dat.nrch, by.x = "pathway", by.y = "ID")
plt.dat.nrch <- merge(dataset.flux,rxn2nrch,
                      by = c("rxn", "phenotype"))
plt.dat.nrch <- unique(plt.dat.nrch[,.(rxn, phenotype, Estimate, `Std. Error`, EffSize, padj,  pathway = pathway.y, Description)])
# order the pathways and remove duplicates
plt.dat.ord <- plt.dat.nrch[,.(od = median(Estimate),mn = mean(Estimate),sel = paste0(phenotype,median(Estimate)*length(rxn))), by = .(pathway, phenotype)]
plt.dat.ord <- plt.dat.ord[!duplicated(sel),]
plt.dat.nrch <- merge(plt.dat.ord, plt.dat.nrch)

pwy.plts <- list()
for(n in names(phenotypes)){
    dat <- plt.dat.nrch[phenotype == n,]
    dat[,Description := factor(Description, levels = unique(Description[order(od)]))]
    if(n != "HBI/Mayo score"){
        xlab = "sqrt(Estimate)"
    } else {
        xlab = "Estimate"
    }
    plt <- ggplot(dat, aes(x=Estimate, y = Description)) +
        geom_vline(xintercept = 0, linetype = 2, color = "red")+
        geom_jitter(width = 0, height = 0.25, alpha = 0.5)+
        geom_point(shape = 24, fill = "blue",size = 2, aes( x =mn))+
        labs(title = n, y = "Pathway", x =xlab )+
        theme_bw()
    pwy.plts[[n]] <- plt
}
plt.all<- cowplot::plot_grid(plotlist= pwy.plts, ncol = 1)
ggsave(plt.all,
       file = "EnrichedPathwaysInternalRxn.pdf",
       height = 6,
       width = 4)

