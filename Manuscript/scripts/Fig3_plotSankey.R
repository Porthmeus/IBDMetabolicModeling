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
met2sub.long <- make_long(met2sub, `Micr. metabolite`, `Host subsystem`)
met2sub.long <- merge(data.table(met2sub.long), unique(met2sub[,.(`Micr. metabolite`, Category2)]), by.x = "node", by.y = "Micr. metabolite", all.x = TRUE)
# adjust names
#met2sub.long[, c("x","next_x") := list(factor(gsub("Mic_met", "Micr. metabolite", gsub("Host_subsystem","Host subsystem", x))),
#                                       factor(gsub("Mic_met", "Micr. metabolite", gsub("Host_subsystem","Host subsystem", next_x))))]
#

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
ggsave(sankey_plot, file = file.path(ddrFig, "Fig3_sankeyplot.pdf"), height = 5, width = 8)
