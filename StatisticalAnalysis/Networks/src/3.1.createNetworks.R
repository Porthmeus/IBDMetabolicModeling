# Porthmeus
# 07.03.23

require(igraph)
require(visNetwork)
require(data.table)
require(qgraph)
source("./createVisNetworks.R")

#
dat.dir <- file.path("..","results")

# get the tested centrality measure
tst.cent <- "degree"
for(st in c("HBMayo","Response","Remission")){
    # split by tissue
    for(tiss in c("blood","biopsy")){
        # tiss <- "biopsy"
        #st <- "HBMayo"
        nodes <- fread(file.path(dat.dir, paste0("basic.nodes.table.1.",st,".",tiss,".csv")))
        edges <- fread(file.path(dat.dir, paste0("basic.edge.table.1.",st,".",tiss,".csv")))
        
        # remove some transport reactions
        edges <- edges[from.id != to.id,]

        # synchronize edge/nodes ids
        nodes.edges <- pruneNetwork(nodes,edges)
        nodes <- nodes.edges[["nodes"]]
        edges <- nodes.edges[["edges"]]

        # calculate degree of each node
        grph <- graph_from_data_frame(edges[,.(from.id, to.id)])
        dgr.grph <- degree(grph)
        nodes[,degree := dgr.grph[met.id]]

        # make distances by degree of the metabolite
        edges[,distance := dgr.grph[to.id] + dgr.grph[from.id]]
        #edges <- edges[,c("distance") := list( -1*abs(EffSize))]

        # calculate neg/pos ratio of associated reactions
        mets.asso <- unique(edges[,.(rxn.id, EffSize)])
        met2rxn <- unique(melt(edges[,.(rxn.id, to.id, from.id)], id.var = "rxn.id", value.name = "met.id")[,.(rxn.id,met.id)])
        mets.asso <- merge(mets.asso, met2rxn)
        mets.asso <- mets.asso[,.(pos = sum(sign(EffSize) == 1),
                              neg = sum(sign(EffSize) == -1),
                              total = .N), by = met.id]
        mets.asso[,pos.neg.ratio := (pos-neg)/total]
        nodes <- merge(nodes, mets.asso[,.(met.id, pos.neg.ratio)], by = "met.id", all.x = TRUE)

        # add tissue and layer of significance to the dropdown menu later on
        nodes[grepl("Metabolomics",src), subsystem:= paste0(subsystem,",Metabolomics")]
        nodes[grepl("Microbiome",src), subsystem:= paste0(subsystem,",Microbiome")]
        nodes[grepl("Centrality",src), subsystem:= paste0(subsystem,",Centrality")]
        
        # replace " " with NA for the src
        nodes[is.na(src), src := " "]

        # save the final nodes/edge table
        fwrite(nodes,
               file = file.path(dat.dir, paste0("final.nodes.table.1.", st, ".", tiss, ".csv")))
        fwrite(edges,
               file = file.path(dat.dir, paste0("final.edges.table.1.", st, ".", tiss, ".csv")))

        # create the nodes mapping
        nodes.final <- nodes[,
                        .( id,
                        label = met.id,
                        color = mapColorToScale(pos.neg.ratio,cols =brewer.pal(name = "PiYG",n=9)[c(1,9)], mid = 0),
#                        color = mapColorToScale(pos.neg.ratio, mid = 0),
                        src = src,
                        shape = mapShapeToFactor(factor(src), "outside"),
                        tissue = tissue,
                        size = mapSizeToScale(degree),
                        pos.neg.ratio,
                        subsystem = subsystem)]
        nodes.final.titel <- nodes[,.( title = paste0("<p><b>",met.name,"</b><br>",
                                       "Tissue: ",tissue,"<br>",
                                       "Changed in: ", src,"<br>",
                                       "Pos/Neg assoc.: ", round(pos.neg.ratio, digits = 2), "<br>",
                                       "Bacarena Eff. Size: ", round(ES.Bacarena, digits = 2), "<br>",
                                       "MicrobiomeGS Eff. Size: ", round(ES.MicrobiomeGS, digits = 2), "<br>",
                                       "Metabolomics Eff. Size: ", round(ES.Metabolomics, digits = 2), "<br>",
                                       "</p>")), by = id]
        if(!all(nodes.final[,id] == nodes.final.titel[,id])){
           stop("Check nodes table for duplicates, can joined data.tables if there are duplicates in the edge list")
        } else {
            nodes.final <- cbind(nodes.final, title = nodes.final.titel[,title])
        }

        # create a table for the legend
        nodes.legend <- data.table(color = c(nodes.final[c(which.max(pos.neg.ratio),which.min(pos.neg.ratio)), color]), size = 25, shape = "dot", label = c("pos. ass. of rxns","neg. ass. of rxns"))
        nodes.legend <- rbind(nodes.legend,
                              data.table(color = brewer.pal(name = "Dark2", n = 3)[1],
                                         label = c(paste("Degree: ",range(nodes[,degree]))),
                                         size = range(nodes.final[,size]),
                                         shape = "dot"
                                         ))
        nodes.legend <- rbind(nodes.legend,
                              cbind(color = brewer.pal(name = "Dark2", n = 3)[1],
                               unique(nodes.final[, .(label = src,size = 25, shape)])))
        
        # opacity value
        op <- "1A"
        edges.final <- edges[,.(
                             to = to,
                             from = from,
                             arrows = arrows,
#                             length = distance,
                             distance = distance,
                             weight = abs(EffSize),
                             subsystem = subsystem,
                             color_dat = mapColorToScale(sign(EffSize),cols =paste0(brewer.pal(name = "PiYG",n=9)[c(1,9)], op) ),
                             ES.dir = sign(EffSize)
                             )]
        # change highlight color
        edges.final[,highlight_dat := gsub(paste0(op, "$"),"FF",color_dat)]
        edges.final[,color := color_dat]
        edges.final.title <- edges[,.(title = paste0("<p><b>",rxn.id,"</b><br>",
                                            reac.name,"<br>",
                                            "Eff.Size: ",round(EffSize,2),"<br>",
                                            "Subsyst.: ",subsystem, "</p>")), by = .(to,from)]
        if(! all(c(edges.final.title[,from] == edges.final[,from], edges.final.title[,to] == edges.final[,to]))){
           stop("Check edges table for duplicates, can joined data.tables if there are duplicates in the edge list")
        } else {
            edges.final <- cbind(edges.final, title =edges.final.title[,title])
        }
#                                
#
        # create a table for the legend
        edge.legend <- unique(edges.final[,.(color= gsub(paste0(op, "$"),"FF",color_dat),
                                             label = c("Negative association","None","Positive association")[ES.dir+2],
                                             width = 3)])

        # get some nicer layouts
        #g <- graph_from_data_frame(edges.final[,.(from,to)])
        #l <- layout_(g, with_())#, component_wise())
    #    l <- qgraph.layout.fruchtermanreingold(edges.final[,.(from,to)], vcount = (nrow(nodes.final)^2)*8, repulse.rad = nrow(nodes.final)^3)
        #nodes.final[,x:= l[,1]]
        #nodes.final[,y:= l[,2]]


        vis.net <- visNetwork(nodes = nodes.final,
                              edges = edges.final,
                              main = paste(gsub("HBMayo","HBI/Mayo score",st),tiss),
                              width = "100%", height = "900px") %>%
            visIgraphLayout(randomSeed = 346, layout = "layout_with_fr", weights = edges.final[,log(distance)]) %>%
            #visIgraphLayout(randomSeed = 346) %>%
            visPhysics(enabled = FALSE) %>%
            visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover =FALSE, labelOnly = FALSE, hideColor = rgb(200,200,200,255*0.001,maxColorValue=255)),
                       nodesIdSelection = list(enabled = TRUE, values = nodes.final[src != "", id], main = "Sig. in other layer"),
                       selectedBy = list( variable = "subsystem", multiple = TRUE, main = "Subsystem/Sig. layers")) %>%
            visLegend(addEdges = edge.legend, addNodes = nodes.legend, position = "right", ncol = 1,width = 0.1, useGroup = FALSE) %>%
            visEdges(width = 2, color = list(color = edges.final[,color_dat],highlight =  edges.final[,highlight_dat]))
        visSave(vis.net,
                file =file.path(dat.dir, paste0(st,".",tiss,".network.1.html")),
                selfcontained = TRUE
        )

    }
}

