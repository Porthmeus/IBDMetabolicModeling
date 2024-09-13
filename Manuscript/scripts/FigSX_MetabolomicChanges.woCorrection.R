# Porthmeus
# 30.01.23

require(data.table)
require(ggplot2)
require(ComplexUpset)
require(cowplot)
require(ggsankey)

# plot the association of serum metabolites with the disease activity, response and remission
ddr <- file.path("..","data","Metabolomics")
ddrHost <- file.path("..","data","Host")
ddrFig <- file.path("..","figures_raw")

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
sig.subs.df <- fread(file.path(ddrHost, "EnrichedSubsystems.woCorrection.csv"))

# expand the set by the clustered reactions
cluster <- fread(file.path(ddr,"metabolites_cluster.csv"))

coef.fls <- list.files(ddr, pattern = "_coefs.csv.*$")
coef.fls <- coef.fls[-grep("LMM", coef.fls)]

### load the data
coefs <- lapply(coef.fls, function(set){
                set.title <- gsub("_coefs","",strsplit(set,
                             split = "\\.")[[1]][2])
                coef.tb <- fread(file.path(ddr, set))[,.(metabolite, Estimate,`Std. Error`, padj)]
                coef.tb <- cbind(set = set.title,
                                 coef.tb)
})
coefs <- do.call(rbind, coefs)
coefs.exp <- merge(coefs, cluster[,.(met, rep.met, sign.cor)], by.x= "metabolite", by.y = "rep.met", suffix = c(".rep",".all"), allow.cartesian=TRUE)
# correct estimates
coefs.exp[, Estimate := sign.cor* Estimate]
# get the corrected names and the associated subsystems
coefs.exp <- merge(coefs.exp[padj < 0.05,], met.meta[,.(colAbbr,Recon3D,Name_metabolomics, Category2)], by.x = "metabolite", by.y = "colAbbr", all.x = TRUE)
coefs.exp <- merge(coefs.exp, unique(met2subsystem[,.(subsystem,metabolite.woComp)]), by.x = "Recon3D", by.y= "metabolite.woComp")
subsystems.metabolomics <- coefs.exp[,.(.N), by =.(set,subsystem)][,.(.N), by = set]

# load the enriched subsystems
enriched <- fread(file.path(ddrHost, "EnrichedSubsystems.woCorrection.csv"))
subsystems.all <- fread(file.path(ddrHost, "subsystems.csv"))

for(st in c("Response","Remission")){
#    st <- "HBMayo"
    # create a plot for the estimates and a barplot for the number of detected metabolites in the different categories
    set.title <- gsub("HBMayo","HBI/Mayo score", st)
    dat.plt <- unique(coefs.exp[set == st,.(Name_metabolomics, Category2, Estimate, `Std. Error`)])
    dat.bar <- dat.plt[,.(N = .N, total = nrow(dat.plt),x = ""), by = Category2]

    plt.bar <- ggplot(dat.bar, aes(y = N,x = x, fill = Category2))+
        geom_bar(stat = "identity", color = "black")+
        labs(x="", y = "No. of metabolites", fill = "Met. Class", title = paste(set.title, "- metabolomics assoc.")) +
        geom_text(aes(label = paste0(round(N/total*100, 2), "%")), position = position_stack(vjust = 0.5), size = 2)+
        scale_fill_brewer(palette = "Dark2")+
        theme_bw()+
        theme(legend.position = "right")
    plt.bar

    # some formatting to get the right plot
    if(set.title != "HBI/Mayo score"){
        dat.plt[,Estimate := sign(Estimate)* sqrt(abs(Estimate))]
        dat.plt[,`Std. Error` := sign(`Std. Error`)* sqrt(abs(`Std. Error`))]
        xlab <- "sqrt(Estimate)"
    } else {
        xlab <- "Estimate"
    }

    # sort the metabolites after estimate
    met.levels <- unique(dat.plt[order(Estimate), Name_metabolomics])
    dat.plt[,Name_metabolomics := factor(Name_metabolomics, levels = met.levels)]

    estimate.plt <- ggplot(dat.plt, aes(x=Estimate, y = Name_metabolomics )) +
                geom_vline(xintercept = 0, linetype = 2, color = "red")+
                geom_errorbar(aes(xmin = Estimate - `Std. Error`, xmax = Estimate + `Std. Error`, color = Category2))+
                geom_point(size = 3, shape = 21, aes(fill  = Category2)) +
                labs(x = xlab, y = "Blood metabolite")+
                scale_color_brewer(palette = "Dark2")+
                scale_fill_brewer(palette = "Dark2")+
                theme_bw() +
                theme(legend.position = "none") 
    estimate.plt


    # create a sankey plot for the 
    dat.plt <- unique(coefs.exp[set == st,.(Metabolite = Name_metabolomics, subsystem,Category2)])
    dat.plt.count <- unique(dat.plt[,.(Metabolite, subsystem)])
    # filter the data to the enriched subsystems
    dat.plt.count <- dat.plt.count[subsystem %in% enriched[set == st, subsystem],]
    # add number of metabolites in each subsystem
    dat.plt.count[,Subsystem := paste0(subsystem, " (", .N, ")"), by = subsystem]

    # bring the data into the right form
    dat.plt.count <- make_long(dat.plt.count, Metabolite, Subsystem)
    # add colors 
    dat.plt.count <- merge(data.table(dat.plt.count), unique(dat.plt[,.(Metabolite, Category2)]), by.x = "node", by.y = "Metabolite", all.x = TRUE)

    sankey_plot <- ggplot(dat.plt.count,
           aes(node = node,
                 next_node = next_node,
                 x=x,
                 next_x = next_x,
                 fill = Category2, label = node)) +
        geom_sankey(color = "black", node.fill = "white", alpha = 0.7)+
        labs(fill = "Met. class", x = "")+
        scale_fill_brewer(palette = "Dark2")+
        geom_sankey_label(size = 3, fill = "white")+
        theme_sankey()+
        theme(legend.position = "none") 
    sankey_plot

    # merge the plots

    plt.up <- plot_grid(estimate.plt,plt.bar, nrow = 1, rel_widths = c(2,1),labels = "AUTO")
    plt.all <- plot_grid(plt.up,  sankey_plot, ncol = 1, labels = c(NA,"C"), rel_heights = c(1,2))
    ggsave(plt.all, file = file.path(ddrFig, paste0("Fig2_MetabolomicChanges.",st,".woCorrection.pdf")),
           height = 3*length(met.levels)*0.1 +3,
           width = 12)

    # create a small test if there where more metabolites found significant in the enriched subsystems than expected by chance
    perm <- 9999
    all.subs  <- unique(subsystems.all[,subsystem])
    sig.subs <- unique(enriched[set==st,subsystem])
    sig.subs.mets <- unique(coefs.exp[set == st,subsystem])
    n.sig.subs.mets <- length(sig.subs.mets)
    n.sig.overlap <- length(intersect(sig.subs, sig.subs.mets))
    
    test.smaller <- rep(NA, perm)
    for(i in 1:perm){
        test.subs <- sample(all.subs, n.sig.subs.mets)
        n.test.overlap<- length(intersect(sig.subs, test.subs))
        test.smaller[i] <- n.test.overlap < n.sig.overlap
    }
    print(paste(st,"p =", 1-sum(test.smaller)/perm))
}
