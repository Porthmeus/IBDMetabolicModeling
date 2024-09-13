# Porthmeus
# 07.06.23

require(data.table)
require(ggplot2)
require(cowplot)

# set some variables
dat.dir <- file.path("..","data", "Microbiome")


phenotypes <- c(`HBI/Mayo score` = "HBMayo",
          Response = "Response",
          Remission = "Remission")

sim.sets <- c("MicrobiomeGS","Bacarena")


datasets <- data.table()
for(s.set in sim.sets){
    flux.cluster <- fread(file.path(dat.dir,s.set, "withinFlux_cluster.csv"))
    for(i in 1:length(phenotypes)){
        # get the right file
        ph.type <- phenotypes[i]
        fl <- list.files(file.path(dat.dir,s.set), pattern = paste0("^withinFlux", "\\.G?LMM\\.", ph.type, ".*", "_coefs.csv"))
        print(file.path(s.set,fl))
        dat <- fread(file.path(dat.dir,s.set, fl))
        ##
        # get the significant reactions and expand the cluster
        dat <- dat[padj < 0.05,]
        dat <- merge(flux.cluster[,.(rep.rxn, rxn, sign.cor)], dat, by.x = "rep.rxn", by.y = "rxn")
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
        # correct the sign of the estimates and the effect size by the correlation results within each cluster
        dat[,Estimate := Estimate*sign.cor]
        dat[,EffSize := Estimate*sign.cor]
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


dataset.flux[,phenotype := factor(phenotype, levels = names(phenotypes))]
dataset.flux[,name := factor(name, levels = unique(name[order(Estimate)]))]
dataset.flux[,xmin := Estimate -`Std. Error`]
dataset.flux[,xmax := Estimate +`Std. Error`]

# capitalize metabolite names
capwords <- function(s) {
    cap <-paste0(toupper(substring(s, 1, 1)),substring(s, 2))
    return(cap)
}
dataset.flux[!is.na(name),name := sapply(name,capwords)]


plts <- list()
rel.height <- c()
for(phntp in dataset.flux[,unique(phenotype)]){
    dat.plt <- dataset.flux[phenotype == phntp,]
    od <- dat.plt[,.(od = mean(Estimate)), by = name]
    dat.plt[,name := factor(name, levels = od[order(od),name])]
    rel.height <- c(rel.height, length(dat.plt[,unique(name)]))
    if(phntp == "HBI/Mayo score"){
        xlab = "Estimate"
    } else {
        xlab = "sqrt(Estimate)"
    }
    within.pl <- ggplot(dat.plt, aes(x = Estimate, y =name, color = sim.set )) +
        scale_color_brewer(palette = "Paired")+
        geom_vline(xintercept = 0, linetype = 2, color = "red") +
        geom_point() +
        geom_errorbar(aes(xmin = xmin, xmax = xmax))+
        theme_bw() +
        labs(x = xlab, y = "Metabolite", title = phntp, color = "Simulation")
    plts[[phntp]] <- within.pl
}
rel.height <- (rel.height/sum(rel.height))+(1/length(plts))
plt.all <- plot_grid(plotlist = plts, ncol = 1, rel_heights = rel.height)

ggsave(plt.all,
       file = file.path("..","figures_raw","FigSX_withinFluxMicrobiom.pdf"),
       height = 5,
       width = 4)
