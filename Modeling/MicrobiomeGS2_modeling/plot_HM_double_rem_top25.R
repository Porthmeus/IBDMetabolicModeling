rm(list=ls())
require(data.table)
require(ggplot2)
setwd("/home/samer/projects/intervention4/figs/")
dou <- fread("/home/samer/projects/intervention4/results/double_effect_inter_log2.csv")
rem <- fread("/home/samer/projects/intervention4/results/rem_effect_inter_log2.csv")

dou.melt <- melt(dou, id.vars = "V1", variable.name = "Intervention", value.name = "log2FC")
dou.melt[,Inter_type := "Doubling"]
rem.melt <- melt(rem, id.vars = "V1", variable.name = "Intervention", value.name = "log2FC")
rem.melt[,Inter_type := "Removal"]
inter <- rbind(dou.melt, rem.melt)
colnames(inter)[1] <- "Target"
inter <- inter[abs(log2FC) > 0.05,]
inf_ind <- which(abs(inter$log2FC)==Inf)
inter[abs(log2FC) > 4&abs(log2FC)!=Inf, log2FC := sign(log2FC)*4]
#inter <- inter[grepl("_o$", Target),]
#inter[, Target := gsub("_o$","", Target)]

# get the significant results
sigs <- fread("/home/samer/projects/intervention4/data/all_targets.csv")

# create desired effect plot
inter.sigs <- merge(inter, sigs, by.x = "Target", by.y = "cpd", all.x = TRUE, allow.cartesian=TRUE)
inter.sigs[,desired := "no"]
inter.sigs[sign(log2FC) == -1*sign(estimate), desired := "yes"]
inf_ind <- which(abs(inter.sigs$log2FC)==Inf)
inter.sigs[inf_ind, log2FC := sign(log2FC)*3]
inter.sigs[inf_ind, desired := "Switch"]

# create heatmap of positive interventions
sel.tab <- inter.sigs[,.(yes_r =sum(desired=="yes"& Inter_type=="Removal"),
                         no_r = sum(desired=="no"& Inter_type=="Removal"),
                         yes_d =sum(desired=="yes"& Inter_type=="Doubling"),
                         no_d = sum(desired=="no"& Inter_type=="Doubling")),
                      by = .(Intervention)]
sel.tab <- sel.tab[order(yes_r + yes_d)]
sel <- sel.tab[(yes_r> no_r)|(yes_d>no_d), unique(Intervention)]
sel <- sel.tab[, unique(Intervention)]
sel <- sel.tab[(yes_r> 0)|(yes_d>0), unique(Intervention)]

metabolites_b <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites_edited.tsv")
#metabolites <- metabolites[,.(id,name = tools::toTitleCase( name))]
metabolites_b <- metabolites_b[, c("id", "name")]
metabolites_h <- metabolites_b
metabolites_h$id <- paste0(metabolites_h$id, "_o") 
metabolites <- rbind(metabolites_b, metabolites_h)
metabolites$name <- ifelse(nchar(metabolites$name) > 16, paste(substr(metabolites$name, 1, 16),
                                                               ".", sep = ""), metabolites$name)

inter.sigs <- merge(inter.sigs, metabolites, by.x = "Intervention", by.y = "id", all.x = TRUE)
inter.sigs <- merge(inter.sigs, metabolites, by.x = "Target", by.y = "id",
                    all.x = TRUE, suffix = c(".intervention",".target"))

inter.sigs.plt <- inter.sigs[Intervention %in% sel & !is.na(estimate),]

inter.sigs.plt <- unique(inter.sigs.plt[,-"estimate"])
inter.sigs.plt <- inter.sigs.plt[!name.intervention %in% c("O2","Arsenate","N2", NA),]
inter.sigs.plt <- inter.sigs.plt[!name.target %in% c("CO2","H2O", NA),]
inter.sigs.plt$order <- match(inter.sigs.plt$Intervention, sel)
inter.sigs.plt <- inter.sigs.plt[order(order)]
inter.sigs.plt[,Direction := "Increase"]
inter.sigs.plt[log2FC<0, Direction := "Decrease"]
inter_same_met_name <- inter.sigs.plt[, .SD[uniqueN(Intervention) > 1], by = name.intervention]
to_remove <- unique(inter_same_met_name$Intervention)[2]
inter.sigs.plt <- inter.sigs.plt[Intervention != to_remove]
host_ind <- grepl ("_o", inter.sigs.plt$Target)
inter.sigs.plt.h <- inter.sigs.plt[host_ind,]
h25 <- tail(unique(inter.sigs.plt.h$name.intervention), 25)
inter.sigs.plt.h <- inter.sigs.plt.h[inter.sigs.plt.h$name.intervention%in%h25,]
inter.sigs.plt.b <- inter.sigs.plt[!host_ind,]
b25 <- tail(unique(inter.sigs.plt.b$name.intervention), 25)
inter.sigs.plt.b <- inter.sigs.plt.b[inter.sigs.plt.b$name.intervention%in%h25,]

#write.csv(inter.sigs.plt.h, "/home/samer/projects/intervention4/results/top25_interv_host.csv")
#write.csv(inter.sigs.plt.b, "/home/samer/projects/intervention4/results/top25_interv_mic.csv")


plt2 <- ggplot(inter.sigs.plt.h, aes(y = factor(name.intervention, levels = unique(inter.sigs.plt.h$name.intervention)),
                                     x = name.target,
                                     fill = desired,
                                     size = abs(log2FC))) +
  geom_point(shape = 21, alpha = 0.8) +
  facet_grid(.~Inter_type) +
  theme_bw() +
  labs(x= "Target",
       y = "Intervention",
       fill = "Desired effect",
       size = "log2(FC)")+
  scale_fill_manual(values=c("#F05954","#9911DC", "#2788D1"))+
  theme(axis.text.x = element_text(angle = 45, hjust =1, size=6),
        axis.text.y = element_text(size=6))

ggsave(
  filename = "/home/samer/projects/intervention4/figs/HMs/rem_double/rem_doub_host_top25.png",
  plot = plt2,
  device = "png",
  height = 4,
  width = 8,
  units = "in",
  dpi = 300  # Adjust the DPI (resolution) as needed
)

plt2 <- ggplot(inter.sigs.plt.b, aes(y = factor(name.intervention, levels = unique(inter.sigs.plt.b$name.intervention)),
                                     x = name.target,
                                     fill = desired,
                                     size = abs(log2FC))) +
  geom_point(shape = 21, alpha = 0.8) +
  facet_grid(.~Inter_type) +
  theme_bw() +
  labs(x= "Target",
       y = "Intervention",
       fill = "Desired effect",
       size = "log2(FC)")+
  scale_fill_manual(values=c("#F05954", "#2788D1"))+
  #    scale_fill_gradient2(low = "blue",mid = "white", high = "red")+
  theme(axis.text.x = element_text(angle = 45, hjust =1, size=6),
        axis.text.y = element_text(size=6))

ggsave(
  filename = "/home/samer/projects/intervention4/figs/HMs/rem_double/rem_doub_bac_top25.png",
  plot = plt2,
  device = "png",
  height = 4,
  width = 8,
  units = "in",
  dpi = 300  # Adjust the DPI (resolution) as needed
)
