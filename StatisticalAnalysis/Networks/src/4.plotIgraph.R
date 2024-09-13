# Porthmeus
# 29.01.24

# create consistent images for the networks

require(data.table)
require(ggplot2)
require(ggraph)
require(igraph)
require(RColorBrewer)
require(cowplot)


dat.dir <- file.path("..","results")
for(st in c("HBMayo","Response","Remission")){
    plt.list <- list()
    for(tis in c("biopsy","blood")){
        if(st == "HBMayo"){
            st.title <- "HBI/Mayo score"
        } else {
            st.title <- st
        }
        #st <- "HBMayo"
        #tis <- "biopsy"
        nodes <- fread(file.path(dat.dir, paste0("basic.nodes.table.1.",st,".",tis,".csv")))
        edges <- fread(file.path(dat.dir, paste0("basic.edge.table.1.",st,".",tis,".csv")))

        # add the information of the enriched metabolites
        enriched.mets <- fread(file.path(dat.dir, "enrichedMetabolites.csv"))
        nodes[met.id %in% enriched.mets[set == st & tissue == tis & sig== TRUE, ID],ES.enrichedMets := TRUE]

        # remove some transporters, which are still in there
        edges <- edges[from.id != to.id,]

        # get only the largest cluster if network is large
        if(nrow(edges) > 5000){
            grph <- graph_from_data_frame(edges[,.(from.id, to.id)])
            clsts <- clusters(grph)
            largest.clst <- which.max(clsts[["csize"]])
            mmbs <- names(clsts[["membership"]])[clsts[["membership"]]==largest.clst]
            nodes <- nodes[met.id %in% mmbs,]
            edges <- edges[from.id %in% mmbs | to.id %in% mmbs,]
        }

        # calculate the degree of metabolites
        grph <- graph_from_data_frame(edges[,.(from.id, to.id)])
        dgr.grph <- degree(grph)
        nodes[,degree := dgr.grph[met.id]]
        # make distances by degree of the metabolite
        edges[,distance := dgr.grph[to.id] + dgr.grph[from.id]]

        # calculate neg/pos ratio of associated reactions
        mets.asso <- unique(edges[,.(rxn.id, EffSize)])
        met2rxn <- unique(melt(edges[,.(rxn.id, to.id, from.id)], id.var = "rxn.id", value.name = "met.id")[,.(rxn.id,met.id)])
        mets.asso <- merge(mets.asso, met2rxn)
        mets.asso <- mets.asso[,.(pos = sum(sign(EffSize) == 1),
                          neg = sum(sign(EffSize) == -1),
                          total = .N), by = met.id]
        mets.asso[,pos.neg.ratio := (pos-neg)/total]
        nodes <- merge(nodes, mets.asso[,.(met.id, pos.neg.ratio)], by = "met.id", all.x = TRUE)

        # get the top 20 metabolites in the degree an give it name in the plot
        # create a list to shorten down the names
        subs <- list("Uridine-5'-Diphosphate" ="UDP",
                    "Uridine Diphosphate" ="UDP",
                    "Nicotinamide Adenine Dinucleotide Phosphate - Reduced" = "NADPH",
                    "Nicotinamide Adenine Dinucleotide Phosphate"="NADP+",
                    "Nicotinamide Adenine Dinucleotide - Reduced" = "NADH",
                    "Nicotinamide Adenine Dinucleotide"="NAD+",
                    "Flavin Adenine Dinucleotide Reduced" = "FADH2",
                    "Flavin Adenine Dinucleotide Oxidized" = "FAD",
                    "Coenzyme A" = "CoA",
                    "Adenosine Triphosphate" = "ATP",
                    "Adenosine Diphosphate" = "ADP",
                    "Adenosine Monophosphate" = "AMP",
                    "Hydrogen Peroxide" = "H2O2",
                    "Diphosphate" = "PP",
                    "\\(.*\\)$" = ""
                     )
        nodes[!is.na(ES.enrichedMets), name_plt := met.name]
        #nodes[order(degree, decreasing=TRUE)[1:20],name_plt :=met.name]
        # do the name shortening
        for(ori in names(subs)){
            nodes[,name_plt:= gsub(ori, subs[[ori]], name_plt)]
        }
        #nodes[order(degree, decreasing=TRUE)[1:20],name_plt]

        # target subsystems ### NOT USED
        target <- c("NAD metabolism")
        edges[subsystem %in% target, target_plt := subsystem]
        ###

        # colors
        cols <- brewer.pal(name = "Dark2", n=6)
        red <- cols[4]
        green <- cols[5]
        blue <- cols[3]
        orange <- cols[2]
        # do some specific stuff for the metabolites
        #nodes.metmic <- copy(nodes)
        #edges.metmic <- copy(edges)

        nodes[,met.mics:= ""]
        nodes[!is.na(ES.Bacarena) | !is.na(ES.MicrobiomeGS), met.mics := paste0(met.mics,";Microbiome")]
        nodes[!is.na(ES.Metabolomics), met.mics := paste0(met.mics, ";Metabolomics")]
        nodes[,met.mics:=gsub("^;","",met.mics)]
        # add a name for the metabolites
        nodes[met.mics != "", name_plt.metmic := met.name]
        metmic.mets <- nodes[met.mics != "", met.id]
        edges[, is.metmic := FALSE]
        edges[from.id %in% metmic.mets | to.id %in% metmic.mets, is.metmic := TRUE]

        # sort the index for the edges and nodes.metmic
        id_first <- which(colnames(edges) %in% c("from.id","to.id"))
        col_idx.edge <- c(colnames(edges)[id_first], colnames(edges)[-id_first])
        id_first <- grep("met.id", colnames(nodes))
        col_idx.node <- c(colnames(nodes)[id_first], colnames(nodes)[-id_first])
        # load the data as graph and add the degree information
        graph <- graph_from_data_frame(d = edges[,mget(col_idx.edge)],
                                       vertices = nodes[,mget(col_idx.node)])

        # create an overview plot
        plt.base <- ggraph(graph,
                      #layout = "linear",
                      #circular = TRUE,
                      layout =  "igraph",
                      algorithm = "with_fr",
                      weights = 1/E(graph)$distance,
                      )

        plt <- plt.base+
            geom_edge_fan(aes(colour = factor(sign(EffSize))),
                          alpha = 0.15)+
           # geom_edge_density(aes(fill = subsystems_plt))+
            scale_edge_color_manual(values = c(red,green))+
            geom_node_point(aes(size = degree,
                                fill = pos.neg.ratio,
                                shape = c("no","yes")[1 + !is.na(name_plt)]),
                            color = "gray"
                            ) +
            geom_node_text(aes(label = name_plt),
                           #color = "gray30",
                           repel = TRUE,
                           min.segment.length = unit(0, "lines"))+
            scale_fill_gradient2(low = red, mid = "white",high = green, midpoint =0) + # adjust colors 
            scale_shape_manual(values = c(21,22))+
            labs(edge_color = "up/down regulation", fill = "pos./neg. assoc. ratio", size = "Degree",title = paste(st.title, tis), shape = "Enriched Metabolite")+
            theme(panel.background = element_blank(),
                  legend.key = element_rect(fill = "transparent")
            )
        plt.list[[paste0(tis, ".net")]] <- plt

        # create highlighting metabolomics/microbiome metabolites

        plt.metmic <-  plt.base+
           # geom_edge_density(aes(fill = subsystems_plt))+
            geom_edge_fan(aes(colour = factor(sign(EffSize) * is.metmic), alpha = 0.01+is.metmic*0.99))+
            scale_edge_color_manual(values = c(red,"gray",green))+
            geom_node_point(aes(size = degree,
                                fill = met.mics,
                                alpha = 0.4+(met.mics!="")*0.6),
                            shape = 21,
                            color = "gray"
                            ) +
            geom_node_text(aes(label = name_plt.metmic),
                           #color = "gray30",
                           repel = TRUE,
                           min.segment.length = unit(0, "lines"))+
            #scale_fill_continuous(low = red, high = green) + # adjust colors 
            scale_fill_manual(values = c("gray",brewer.pal(name = "Dark2", n = 6)))+
            labs(edge_color = "up/down regulation", fill = "Sig. association in", size = "Degree",title = paste(st.title, tis, "micr./met. assoc."))+
            guides(alpha = "none", edge_alpha = "none")+
            theme(panel.background = element_blank(),
                  legend.key = element_rect(fill = "transparent")
            )
        plt.list[[paste0(tis, ".metmic")]] <- plt.metmic
    }
    plt.all <- plot_grid(plotlist = plt.list,  ncol = 2)
    ggsave(plt.all, file = file.path(dat.dir, paste0("ggraphGridPlot.",st,".jpg")),dpi=600,
           width = 24, height = 16, )
}
