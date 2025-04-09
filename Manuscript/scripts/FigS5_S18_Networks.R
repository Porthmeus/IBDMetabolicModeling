# Porthmeus
# 29.01.24

# create consistent images for the networks

require(data.table)
require(ggplot2)
require(ggraph)
require(igraph)
require(RColorBrewer)
require(cowplot)


# directories
dat.dir <- file.path("..","data","Networks")
res.dir <- file.path("..","figures_raw")

for(st in c("Response","Remission","HBMayo")){
    plt.list.networks <- list()
    plt.list.paths <- list()
    counter <- -2
    for(tis in c("biopsy","blood")){
        counter <- counter +2
        if(st == "HBMayo"){
            st.title <- "HBI/Mayo score"
        } else {
            st.title <- st
        }
        #st <- "HBMayo"
        #tis <- "blood"
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
                    "Diphosphate" = "PPi",
                    "D-Galactosyl"= "D-Gal",
                    "Digalactosyl" = "Di-Gal",
                    "Phosphatidylethanolamine" = "PE",
                    "\\(.*\\)$" = "",
                    "(Sm.*), Sphingomyelin"= "\\1",
                    "^.*, Lysopc"="LysoPC"
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
        # correct some lengthy names
        for(sn in names(subs)){
            nodes[,name_plt.metmic := gsub(sn, subs[[sn]], name_plt.metmic)]
        }
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

        # calculate path distances between all nodes
        pathdist <- distances(graph, weights = log(E(graph)$distance), mode = "in")
        # reformat

        pathdist <- melt(data.table(pathdist, keep.rownames = TRUE),
                         id.var ="rn",
                         value.name = "distance",
                         variable.name = "to")
        colnames(pathdist)[1] <- "from"
        # remove self distances and impossible paths
        pathdist <- pathdist[distance < Inf & from != to,]

        # calculate the lower 5% cutoff and look for enrichment of low distances paths between different categories
        cutoff <- 0.05
        nperm <- 100
        thresh <- quantile(pathdist[from != to & distance <Inf,distance],cutoff)
        n.all <- nrow(pathdist)
        n.low <- nrow(pathdist[distance <= thresh, ]) # will need that number later for the calculation
        met.sel <- list("Hub" = nodes[ES.enrichedMets== TRUE,met.id],
                        "Mic" = c(nodes[!is.na(ES.Bacarena) | !is.na(ES.MicrobiomeGS), met.id]),
                        "Met" = nodes[!is.na(ES.Metabolomics), met.id]
                        )
        cmprs.pairs <- combn(names(met.sel),2)

        stats <- data.table()
        dists.sel<- data.table()
        for(i in 1:ncol(cmprs.pairs)){
            met.sel1 <- met.sel[[cmprs.pairs[1,i]]]
            met.sel2 <- met.sel[[cmprs.pairs[2,i]]]
            pathdist.sel <- pathdist[(from %in% met.sel1 & to %in% met.sel2) | (from %in% met.sel2 & to %in% met.sel1),]
            non.pathdist.sel <- pathdist[!((from %in% met.sel1 & to %in% met.sel2) | (from %in% met.sel2 & to %in% met.sel1)),]
            dists.sel <- unique(rbind(dists.sel, 
                                      cbind(pathdist.sel, con.sets = paste0(cmprs.pairs[1,i],"-",cmprs.pairs[2,i])),
                                      cbind(non.pathdist.sel, con.sets = "all")
                                      )
                        )
            if(nrow(pathdist.sel)>3){
                t.stat <- t.test(pathdist.sel[,distance], non.pathdist.sel[,distance], alternative = "less")
                stats <- rbind(stats,
                               data.table(pair1 = cmprs.pairs[1,i],
                                          pair2 = cmprs.pairs[2,i],
                                          t.val = t.stat[["statistic"]],
                                          log2FC = log2(t.stat[["estimate"]][1]/t.stat[["estimate"]][2]),
                                          conf.int = abs(t.stat[["conf.int"]])[1]-abs(diff(t.stat[["estimate"]])),
                                          std.err = t.stat[["stderr"]],
                                          p.value = t.stat[["p.value"]]))
            } else {
                stats <- rbind(stats,
                               data.table(pair1 = cmprs.pairs[1,i],
                                          pair2 = cmprs.pairs[2,i],
                                          t.val = 0,
                                          log2FC = 0,
                                          conf.int = Inf,
                                          std.err = Inf,
                                          p.value = 1)) 
            }
        }
        stats[,padj := p.adjust(p.value, method = "BH")]
        fwrite(stats, file = file.path("..","temp", paste0(st,".",tis, ".NetworkSetDistancesStats.csv")))

        stats[,c("con.sets", "distance") := list(paste0(pair1,"-",pair2), max(dists.sel[,distance]*0.9))]
        stats[,label2 := ""]
        stats[padj < 0.05,label2 := "*"]

        # create a small boxplot comparing the distances

        dist.plt <- ggplot(dists.sel, aes(x=con.sets, y = distance))+
            geom_boxplot() +
            geom_text(data=stats,
                      #angle = 45,
                      size = 5,
                      aes(label =label2))+# paste0("p = ",formatC(padj, format = "e", digits = 2))))+
            labs(x="", y = "Path length", title = tools::toTitleCase(tis))+
            theme_bw()+
            theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))

        # calculate the shortest paths between the low distance values
        dists.sel.low <- dists.sel[con.sets != "all" & distance < thresh,]
        # add the common names
        dists.sel.low <- merge(dists.sel.low, nodes[,.(met.id, met.name)], by.x = "from", by.y = "met.id")
        dists.sel.low <- merge(dists.sel.low, nodes[,.(met.id, met.name)], by.x = "to", by.y = "met.id", suffix = c(".from",".to"))

        for(sn in names(subs)){
            dists.sel.low[,c("met.name.from","met.name.to") := list(gsub(sn,subs[[sn]], met.name.from),
                                                                   gsub(sn,subs[[sn]], met.name.to))]
        }
        # sort by clustering
        if(nrow(dists.sel.low) >3){
            dist.sel.low.mat <- dcast(data=dists.sel.low, met.name.from~met.name.to, value.var = "distance",fun.aggregate = mean, fill = 0)
            dists.sel.low.rn <- dist.sel.low.mat[[1]]
            dists.sel.low.cn <- colnames(dist.sel.low.mat)[-1]
            dist.sel.low.mat <- as.matrix(dist.sel.low.mat[,-1])
            od.from <- dists.sel.low.rn[hclust(dist(dist.sel.low.mat))[["order"]]]
            od.to <- dists.sel.low.cn[hclust(dist(t(dist.sel.low.mat)))[["order"]]]
            dists.sel.low[,c("met.name.from","met.name.to") := list(factor(met.name.from, levels = od.from),
                                                  factor(met.name.to, levels = od.to))]
        
            fromto.plt <- ggplot(dists.sel.low, aes(x = met.name.from, y = met.name.to)) +
                   geom_point(aes(size = 1/distance, color = distance), alpha = 0.7)+ 
                   facet_grid(~con.sets, scale = "free", space = "free") +
                   scale_color_gradient(low = "red", high = "gray")+
                   guides(size = "none")+
                   theme_bw() +
                   labs(x = "From",y = "To", title = paste(tools::toTitleCase(tis), "path length < 5% quant."))+
                   theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))

        } else {
            fromto.plt <-ggplot(dists.sel.low)
        }

       # table(c(a[,from],a[,as.character(to)]))

        

        # create an overview plot
        plt.base <- ggraph(graph,
                      #layout = "linear",
                      #circular = TRUE,
                      layout =  "igraph",
                      algorithm = "with_fr",
                      weights = log(E(graph)$distance),
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
                           max.overlaps= 20,
                           min.segment.length = unit(0, "lines"))+
            scale_fill_gradient2(low = red, mid = "white",high = green, midpoint =0) + # adjust colors 
            scale_shape_manual(values = c(21,22))+
            labs(edge_color = "up/down regulation", fill = "pos./neg. assoc. ratio", size = "Degree",title = paste("   ",st.title, tis), shape = "Enriched Metabolite")+
            theme(panel.background = element_blank(),
                  legend.key = element_rect(fill = "transparent")
            )
        #plt

        plt.list.networks[[paste0(tis, ".net")]] <- plt

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
            labs(edge_color = "up/down regulation", fill = "Sig. association in", size = "Degree",title = paste("   ",st.title, tis, "micr./met. assoc."))+
            guides(alpha = "none", edge_alpha = "none")+
            theme(panel.background = element_blank(),
                  legend.key = element_rect(fill = "transparent")
            )
            #plt.metmic

        plt.list.networks[[paste0(tis, ".metmic")]] <- plt.metmic
        # add the other two plots for the path distances
        plt.dists <- plot_grid(dist.plt,fromto.plt, rel_widths = c(1,5),nrow = 1,labels = LETTERS[c(1,2)+counter])

        plt.list.paths[[paste0(tis, ".dists")]] <- plt.dists
    }
    plt.all.networks <- plot_grid(plotlist = plt.list.networks, byrow = FALSE, ncol = 2,labels = "AUTO")
    plt.all.paths <- plot_grid(plotlist = plt.list.paths, ncol = 1)
   # ggsave(plt.all.paths, file = file.path(res.dir, paste0("Fig5_pathwayPlots.",st,".pdf")),
   #        width = 14, height = 12)
    ggsave(plt.all.networks, file = file.path(res.dir, paste0("FigS5_ggraphGridPlot.",st,".pdf")),#dpi=600,
           width = 30, height = 16)
}
