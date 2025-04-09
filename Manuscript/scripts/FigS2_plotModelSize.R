# Porthmeus
# 13.03.25

require(data.table)
require(ggplot2)
require(GGally)
require(cowplot)
#require(lme4)
#require(lmerTest)


# directories
rsrc.biopsy <- file.path("..","data","Host","biopsy")
rsrc.blood <- file.path("..","data","Host","blood")
rsrc.mic <- file.path("..","data","Microbiome")
res.dir <- file.path("..","figures_raw")


# files 
dat.biopsy <- fread(file.path(rsrc.biopsy, "ModelSize.biopsy.csv"))
dat.blood <- fread(file.path(rsrc.blood, "ModelSize.blood.csv"))
dat.mic <- fread(file.path(rsrc.mic, "community_stats.csv"))


# combine files and plot host data
dat<-rbind(dat.biopsy, dat.blood)

hist.plt <- ggplot(dat, aes( x=modelSize))+
    geom_density(fill = "grey", alpha = 0.7)+
   # geom_histogram(bins = 30, alpha =0.7)+
    theme_bw()+
    facet_grid(~Tissue)+
    geom_vline(data = dat[,.(meanModelSize = mean(modelSize)), by = Tissue],aes(xintercept = meanModelSize), color = "red", linetype = 2)+
    labs(x = "Model size", y = "Count")
hist.plt



dat[,.(meanModelSize = mean(modelSize),
       sdModelSize = sd(modelSize)), by = Tissue]

# plot microbial data
facet_names <-c("Community_ID" = "Community_ID",
                    "Number_of_bacterial_models" = "No. bacterial models",
                    "Genra" = "Genera",
                    "Flux_to_host" = "No. of met. produced",
                    "Flux_from_diet" = "No. of met. consumed",
                    "Flux_within_microbiome" = "No. of cross feeding met.",
                    "Internal_reactions" = "Community model size")
colnames(dat.mic) <- facet_names[colnames(dat.mic)]
mic.dat <- ggpairs(dat.mic, columns = 2:ncol(dat.mic), aes(alpha = 0.7)) +
    theme_bw()

plotAll<- plot_grid(hist.plt, ggmatrix_gtable(mic.dat),ncol = 1, rel_heights = c(1,4),labels = "AUTO")
ggsave(plotAll, file = file.path(res.dir, "FigS2_ModelSizes.pdf"), width = 10, height =12)

dat.mic[,.(No.models = mean(`No. bacterial models`),
           No.genera = mean(Genera),
           No.Flux_to_host = mean(`No. of met. produced`, na.rm = TRUE),
           No.Flux_from_diet = mean(`No. of met. consumed`, na.rm = TRUE),
           No.Flux_within = mean(`No. of cross feeding met.`, na.rm = TRUE),
           Internal = mean(`Community model size`, na.rm = TRUE),
           No.models.sd = sd(`No. bacterial models`, na.rm = TRUE),
           No.genera.sd = sd(Genera, na.rm = TRUE),
           No.Flux_to_host.sd = sd(`No. of met. produced`, na.rm = TRUE),
           No.Flux_from_diet.sd = sd(`No. of met. consumed`, na.rm = TRUE),
           No.Flux_within.sd = sd(`No. of cross feeding met.`, na.rm = TRUE),
           Internal.sd = sd(`Community model size`, na.rm = TRUE))]
