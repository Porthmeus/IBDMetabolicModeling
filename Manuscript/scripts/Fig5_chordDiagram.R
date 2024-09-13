require(data.table)
require(ggplot2)
require(ggraph)
require(igraph)
require(tidygraph)
require(RColorBrewer)
require(cowplot)


# directories
dat.dir <- file.path("..","data","Networks")
res.dir <- file.path("..","figures_raw")

circle.plts <- list()
for(st in c("Response","Remission","HBMayo")){
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
        #tis <- "biopsy"
        nodes <- fread(file.path(dat.dir, paste0("basic.nodes.table.1.",st,".",tis,".csv")))
        edges <- fread(file.path(dat.dir, paste0("basic.edge.table.1.",st,".",tis,".csv")))

        # add the information of the enriched metabolites
        enriched.mets <- fread(file.path(dat.dir, "enrichedMetabolites.csv"))
        nodes[met.id %in% enriched.mets[set == st & tissue == tis & sig== TRUE, ID],ES.enrichedMets := TRUE]

        # calculate the degree of metabolites
        grph <- graph_from_data_frame(edges[,.(from.id, to.id)])
        dgr.grph <- degree(grph)
        nodes[,degree := dgr.grph[met.id]]
        # make distances by degree of the metabolite
        edges[,distance := dgr.grph[to.id] + dgr.grph[from.id]]
        
        # define metabolites which are from different pools of association
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
        
        # get those pathways below a certain distance
        stats <- data.table()
        dists.sel<- data.table()
        for(i in 1:ncol(cmprs.pairs)){
            met.sel1 <- met.sel[[cmprs.pairs[1,i]]]
            met.sel2 <- met.sel[[cmprs.pairs[2,i]]]
            pathdist.sel <- pathdist[(from %in% met.sel1 & to %in% met.sel2) | (from %in% met.sel2 & to %in% met.sel1),]
            dists.sel <- unique(rbind(dists.sel, 
                                      cbind(pathdist.sel, con.sets = paste0(cmprs.pairs[1,i],"-",cmprs.pairs[2,i])))
                        )
        }
        dists.sel <- dists.sel[distance <= thresh,]
        if(nrow(dists.sel) > 2){
            # create another annotation data frame for the metabolites
            met.sel.dt <- data.table()
            for(n in names(met.sel)){
                met.sel.dt <-rbind(met.sel.dt, 
                                   data.table(group = n,
                                    metabolite = met.sel[[n]]))
            }

            # add the information
            dists.sel <- merge(dists.sel, met.sel.dt, by.x = "from",by.y = "metabolite", allow.cartesian = TRUE)
            dists.sel <- merge(dists.sel, met.sel.dt, by.x = "to",by.y = "metabolite",suffix = c(".from",".to"), allow.cartesian = TRUE)
            dists.sel[,c("to","from"):= list(paste0(group.to, "_",to), paste0(group.from,"_",from))]
            dists.sel <- dists.sel[group.from != group.to,]
            met.sel.dt[,group.met := paste0(group, "_", metabolite)]
            met.sel.dt[,.(group.met, group, metabolite)]

            # add full names of metabolites
            met.sel.dt <- merge(met.sel.dt, nodes[, .(met.id, met.name)], by.x = "metabolite", by.y = "met.id")
            # shorten common names
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
            for(ori in names(subs)){
                met.sel.dt[,met.name:= gsub(ori, subs[[ori]], met.name)]
            }

            # sort some things out to get the correct layout
            met.sel.dt[, group.met:= factor(group.met)]
            met.sel.dt <- met.sel.dt[order(group.met),]
            met.sel.dt[, rel_id := (1:nrow(met.sel.dt))/nrow(met.sel.dt)]
            met.sel.dt[, angle := (rel_id *-360)-90]
            met.sel.dt[angle>-270, angle := angle-180]
            met.sel.dt[,angle:= angle+360]
            met.sel.dt[, xnudge := sin((2*rel_id)*pi)]#*sign(rel_id*-360+180)]
            met.sel.dt[, ynudge := cos((2*rel_id)*pi)]
            met.sel.dt[, hjust := as.integer(xnudge <=1e-12)]#(-sign(xnudge)+1)/2]
            


            
            # create the layout
            # determine margin by length of strings
            met.lim <- met.sel.dt
            met.lim[,c("name_length","cos", "sin","tan") := list(lengths(strsplit(met.name, "")),
                                                                    cos(rel_id*2*pi),
                                                                    sin(rel_id*2*pi),
                                                                    tan(rel_id*2*pi))]

            lim_y_max <- met.lim[which.max(name_length * sign(cos)*cos^2),name_length]
            lim_y_min <- met.lim[which.max(name_length * -1*sign(cos)*cos^2),name_length]
            lim_x_max <- met.lim[which.max(name_length * sign(sin)*sin^2),name_length]
            lim_x_min <- met.lim[which.max(name_length *-1* sign(sin)*sin^2),name_length]
            add <- 1
            fact <- 18
            fact2 <- 36
            lim_y_max <- add+(lim_y_max)/(fact2+lim_y_max)
            lim_y_min <- add+(lim_y_min)/(fact+lim_y_min)
            lim_x_max <- add+(lim_x_max)/(fact+lim_x_max)
            lim_x_min <- add+(lim_x_min)/(fact2+lim_x_min)

            # create the chord diagram
            met.sel.dt[,group_long := gsub("Met","Metabolome",gsub("Mic","Microbiome",group))]
            sel <- unique(c("group.met", colnames(met.sel.dt)))
            grph <- graph_from_data_frame(dists.sel, directed = TRUE ,vertices = met.sel.dt[,..sel])
            grph <- as_tbl_graph(grph)

    #        grph.layout<- create_layout(grph,layout = "linear", circular = TRUE) 
            circle.plt <- ggraph(grph, layout = "linear", circular = TRUE) +
                scale_edge_alpha('Edge direction', guide = 'edge_direction')+
                geom_edge_arc(#alpha = 0.5,
                              aes(edge_colour = con.sets,
                                  alpha = after_stat(index)))+
                #scale_edge_color_brewer(palette = "Set1")+
                scale_edge_color_manual(values=brewer.pal(n=6, name = "Dark2")[-(1:3)])+
                geom_node_point(aes(fill=group_long), size = 4, shape = 21)+
                geom_node_text(size = 3, aes(label = met.name,
                                   angle = angle,
                                   x = x+xnudge*0.04,
                                   y = y+ynudge*0.04,
                                   hjust = hjust
                                   ))+
                ylim(-lim_y_min,lim_y_max)+
                xlim(-lim_x_min,lim_x_max)+
                scale_fill_brewer(palette="Dark2")+
                scale_color_brewer(palette="Dark2")+
                coord_fixed()+
                labs(title = tools::toTitleCase(tis), fill="Metabolite", edge_colour = "Metabolite connection")+
                theme(panel.background = element_blank(),
                      legend.key = element_rect(fill = "transparent"))
        circle.plt

        circle.plts[[paste(st, tis, sep = ".")]] <- circle.plt
        }
    }
}



boxplot.plts <- list()
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

        # do some specific stuff for the metabolites
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
                      size = 8,
                      aes(label =label2))+# paste0("p = ",formatC(padj, format = "e", digits = 2))))+
            labs(x="", y = "Path length", title = tools::toTitleCase(tis))+
            theme_bw()+
            theme(axis.text.x = element_text(angle=45, hjust = 1, vjust=1))
        boxplot.plts[[paste0(st, tis, ".")]] <- dist.plt
    }
}

# combine the plots
# HBMayo
hbmayo.plt <- plot_grid(boxplot.plts[[5]], circle.plts[[3]],boxplot.plts[[6]], circle.plts[[4]],
                        ncol = 2, rel_widths=c(1,5), labels = "AUTO")
ggsave(hbmayo.plt, file = file.path(res.dir, paste0("Fig5_pathwayPlots.HBMayo.pdf")),
           width = 14, height = 16)
response.plt <- plot_grid(boxplot.plts[[1]], circle.plts[[1]],boxplot.plts[[2]], circle.plts[[2]],
                        ncol = 2, rel_widths=c(1,5), labels = "AUTO")
ggsave(response.plt, file = file.path(res.dir, paste0("FigS10_pathwayPlots.Response.pdf")),
           width = 14, height = 16)
