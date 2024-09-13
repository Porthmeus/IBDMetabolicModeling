# Porthmeus
# 01.12.23

# plot the values for HB/Mayo over time

require(data.table)
require(ggplot2)
require(cowplot)

# load the data
dat.dir <- file.path("..","data","Host","Cohort")
plt.dir <- file.path("..","figures_raw")
clinic <- fread(file.path(dat.dir, "ClinicalDataMergeLongImputeCor.csv"))

# add information about imputation, cohort and normalize
clinic[,c("imputed","Cohort"):= list(factor(c("no","yes")[1+(is.na(HB_Mayo))], levels = c("yes","no")),
                                     gsub("^(.*)_.*","\\1",PatientID))]
clinic[,HB_Mayo_norm := HB_Mayo_impu/max(HB_Mayo_impu), by = "Diagnosis"]


resp.plt <- ggplot(clinic, aes(x = Time_seq, y = HB_Mayo_norm)) +
    geom_smooth(method = "lm") +
    geom_point(aes(color = imputed)) +
    theme_bw() +
    scale_color_brewer(palette = "Paired")+
    facet_grid(Diagnosis ~ Responder) +
    labs(x = "Time [d]", y = "normalized HBI/Mayo index", title = "Response against Diagnosis")


rem.plt <- ggplot(clinic, aes(x = Time_seq, y = HB_Mayo_norm)) +
    geom_smooth(method = "lm") +
    geom_point(aes(color = imputed)) +
    theme_bw() +
    scale_color_brewer(palette = "Paired")+
    facet_grid(Diagnosis ~ Remision) +
    labs(x = "Time [d]", y = "normalized HBI/Mayo index", title = "Remission against Diagnosis")


cohort.plt <- ggplot(clinic, aes(x = Time_seq, y = HB_Mayo_norm)) +
    geom_smooth(method = "lm",aes(color = Remision)) +
    geom_point(aes(color = Remision),
               #size = 0.75,
               alpha =0.5,
               stroke = NA
               ) +
    theme_bw() +
    ylim(0,1)+
    scale_color_brewer(palette = "Dark2")+
    facet_grid(Diagnosis ~ Cohort) +
    labs(x = "Time [d]",
         y = "normalized HBI/Mayo index",
        # title = "Cohort against Diagnosis",
         color = "Remission")

cohort.plt2 <- ggplot(clinic, aes(x = Time_seq, y = HB_Mayo_norm)) +
    geom_smooth(method = "lm",aes(color = Responder)) +
    geom_point(aes(color = Responder),
             #  size = 0.75,
               alpha = 0.5,
              # shape = 21,
               stroke = NA) +
    theme_bw() +
    ylim(0,1)+
    scale_color_brewer(palette = "Dark2")+
    facet_grid(Diagnosis ~ Cohort) +
    labs(x = "Time [d]",
         y = "normalized HBI/Mayo index",
         #title = "Cohort against Diagnosis",
         color = "Response")
cohort.plt2

plt.all <- plot_grid( cohort.plt2,cohort.plt,labels = "AUTO")

ggsave(plt.all,
       file = file.path(plt.dir, "FigSX_HBMayoVsTime.pdf"),
       width = 10, height = 4)
