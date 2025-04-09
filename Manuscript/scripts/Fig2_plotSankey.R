# Porthmeus
# 20.05.24

require(data.table)
require(ggplot2)
require(ggsankey)
require(cowplot)

ddrHost <- file.path("..","data","Host")
ddrFig <- file.path("..","figures_raw")

# read data
met2sub <- fread(file.path(ddrHost, "MicrobiomeMets2HostSubsystems.csv"))
met2sub[,Host_subsystem := paste0(Host_subsystem, " (", .N,")"), by = Host_subsystem]
colnames(met2sub) <- gsub("Mic_met", "Micr. metabolite", gsub("Host_subsystem","Host subsystem", colnames(met2sub)))

#### Data table looks like this ###
#     Micr. metabolite           Host subsystem     Category2
#  1:        L-Leucine   Peptide metabolism (1)     Amino acids and related
#  2:        L-Leucine        Miscellaneous (2)     Amino acids and related

# the rest is bringing the data into the right form
met2sub.long <- make_long(met2sub, `Micr. metabolite`, `Host subsystem`)
met2sub.long <- merge(data.table(met2sub.long), unique(met2sub[,.(`Micr. metabolite`, Category2)]), by.x = "node", by.y = "Micr. metabolite", all.x = TRUE)

sankey_plot <- ggplot(met2sub.long,
       aes(node = node,
             next_node = next_node,
             x=x,
             next_x = next_x,
             fill = Category2,
             label = node)) +
    geom_sankey(color = "black", node.fill = "white", alpha = 0.7)+
    labs(fill = "Metabolite class", x = "")+
    scale_fill_brewer(palette = "Dark2")+
    geom_sankey_label(size = 3, fill = "white")+
    theme_sankey()
ggsave(sankey_plot, file = file.path(ddrFig, "Fig2_sankeyplot.pdf"), height = 5, width = 8)
