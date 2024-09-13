# Porthmeus
# 30.01.23

require(data.table)
require(ggplot2)
require(ComplexUpset)
require(cowplot)

# plot the association of serum metabolites with the disease activity, response and remission
ddr <- file.path("..","data","Metabolomics")
tmp.dir <- file.path("..","temp")

# Annotations
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites_edited.tsv")
setkey(metabolites, "id")
met.meta <- fread(file.path(ddr,"Metabolomics2Recon3D.csv"))
met.categories <- fread(file.path(ddr, "metabolome_names.csv"))
met.meta <- merge(met.meta, met.categories[,.(Abbreviation, Category2)], by = "Abbreviation", all.x = TRUE)
met.meta[,colAbbr := gsub("(^\\d)","X\\1",gsub("[\\(|\\.|:| |\n |\\)|/|-]",".",Abbreviation))]
met.meta <- met.meta[!duplicated(colAbbr),.(colAbbr, Name, Category2)]
met.meta[,Name_metabolomics := Name]

# expand the set by the clustered reactions
cluster <- fread(file.path(ddr,"metabolites_cluster.csv"))

coef.fls <- list.files(ddr, pattern = "_coefs.csv.*$")
coef.fls <- coef.fls[-grep("crosslinks", coef.fls)]

plts <- list()
rel.heights <- c()
for(set in rev(coef.fls)){

    # read coefficients and create a data.frame for plotting
    coefs <- fread(file.path(ddr,set))
    coefs.exp <- merge(coefs, cluster[,.(met, rep.met,sign.cor)], by.x= "metabolite", by.y = "rep.met", suffix = c(".rep",".all"))
    # correct estimate by sign of the cluster correlation
    coefs.exp[, Estimate := Estimate*sign.cor]
    coefs.exp.plt <- merge(coefs.exp[padj < 0.05,], met.meta[,.(colAbbr,Name_metabolomics,Category2)], by.x = "metabolite", by.y = "colAbbr", all.x = TRUE)
    coefs.exp.plt[,Name_metabolomics:= factor(Name_metabolomics, levels = unique(Name_metabolomics[order(Estimate)]))]
    write.csv(file = file.path(tmp.dir,paste("NAin_",set, sep ="")),
              coefs.exp.plt[is.na(Name_metabolomics),])

    # check if there are any metabolomics associations
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
        rel.heights <- c(rel.heights, nrow(coefs.exp.plt))

        # create the actual plot and store in a list
        plts[[set]] <- 
            ggplot(coefs.exp.plt,aes(x=Estimate, y = Name_metabolomics )) +
            geom_vline(xintercept = 0, linetype = 2, color = "red")+
            geom_point(size = 3, aes(color = Category2)) +
            labs(x = xlab, y = "Blood metabolite", title = paste(set.title, "association to blood metabolomics"))+
            geom_errorbar(aes(xmin = Estimate - `Std. Error`, xmax = Estimate + `Std. Error`, color = Category2))+
            scale_color_brewer(palette = "Paired")+
            theme_bw() +
            theme(legend.position = "none") 
        plts[[paste0(set,"_bar")]] <- 
            ggplot(coefs.exp.plt, aes( y = Category2)) +
            geom_bar(aes(fill = Category2)) +
            scale_fill_brewer(palette = "Paired")+
            guides(fill = "none")+
            theme_bw() +
            labs(y = "Compound class", fill = "Compound class",title = paste("#",  "classes"))

        }
}


plt <- plot_grid(plotlist= plts, labels = "AUTO", ncol = 2, rel_heights = rel.heights, rel_widths = c(rep(c(2,1),3)))
plt

ggsave(plt, 
       file = file.path("..","figures_raw", "FigS1_metabolomicsDiseaseAssociationSingle.pdf"),
       height = 20,
       width = 10)
