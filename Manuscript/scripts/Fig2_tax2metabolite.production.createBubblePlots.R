# Porthmeus
# 08.12.23

require(data.table)
require(ggplot2)

dat.dr <- file.path("..","data","Microbiome", "MetaboliteProductionPerOrganism")
res.dr <- file.path("..","figures_raw")

# load the plotting data
met.prod <- unique(fread(file.path(dat.dr, "tax2metabolite.csv.gz")))

# filter for non-significant effect sizes and add metabolite names
met.prod.sig <- met.prod[sign(LCL) == sign(UCL),]
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")
met.prod.sig <- merge(met.prod.sig, metabolites[,.(id,name)], by.x = "cpd", by.y = "id", all.x = TRUE)
met.prod.sig[,name := gsub("indol","Indol",name)]
met.prod.sig[,name := gsub("isocholate","Isocholate",name)]

# rename micorbiomeGS to MGS to fit the plot
met.prod.sig[sim.set == "MicrobiomeGS", sim.set := "MGS"]
# correct indol spelling
met.prod.sig[,name := gsub("ndol$","ndole", name)]
# add the information for neg/pos assocation to the original data
met.prod.sig[t.value.ori.assoc > 0,pos.neg.asso := "high infl."]
met.prod.sig[t.value.ori.assoc < 0,pos.neg.asso := "low infl."]

# create a plot per tax level
tax.levels <- met.prod.sig[,unique(tax.levels)]


for(tx.lv in tax.levels){
    # seperate host metabolites and internal microbial metabolites
    for(flux in c("outflux","withinFlux")){

        # create the correct title
        if(flux == "outflux"){
            ttl <- "Host exchanged metabolites"
        } else{
            ttl <- "Microbial exchanged metabolites"
        }
        # cluster taxa for plotting
        dat.plt <- met.prod.sig[tax.levels == tx.lv & target.set == flux,]
        dat.cluster <- dcast(dat.plt, new.level~name, value.var = "emmean", fill = 0,fun.aggregate = "mean")
        clst <- hclust(dist(dat.cluster[,-1]))
        lvls <- dat.cluster[[1]][clst$order]
        dat.plt[,new.level2 := factor(new.level, levels = lvls)]
        # cluster metabolites
        clst <- hclust(dist(t(dat.cluster[,-1])))
        lvls <- colnames(dat.cluster)[-1][clst$order]
        dat.plt[,name := factor(name, levels = lvls)]
        
        # creat the actual plot
        pl <- ggplot(dat.plt, aes(x = name, y = new.level2)) +
            geom_point(alpha = 0.5,shape = 21, aes(size = log10(abs(emmean)), fill = factor(sign(emmean))))+
            theme_bw() +
            facet_grid( ~ sim.set+pos.neg.asso, scale = "free", space = "free")+
            scale_fill_discrete(type = c("#52a0da","#f37a76"), labels=c("uptake","production"))+
            scale_size_continuous(name = "log10(|ES|)")+
            theme(axis.text.x = element_text(angle = 45, hjust = 1))+
            labs(x=element_blank(), y = tx.lv, fill = "Direction", title = ttl)
        pl
        ggsave(pl, 
               file = file.path(res.dr, paste("tax2metabolite",tx.lv,flux,"pdf", sep = ".")),
               width = 0.2* length(dat.plt[,unique(paste0(sim.set, name))])+3,
               height = 0.12*length(levels(dat.plt[,new.level2]))+1.5)

    }
}


