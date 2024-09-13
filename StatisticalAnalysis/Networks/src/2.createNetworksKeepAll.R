# Porthmeus
# 07.03.23

require(igraph)
require(visNetwork)
require(data.table)
source("./createVisNetworks.R")

dat.dir <- file.path("..","results")
nodes <- fread(file.path(dat.dir, "basic.nodes.table.csv"))
edges <- fread(file.path(dat.dir, "basic.edge.table.csv"))

# tested centrality measure
tst.cent <- "degree"
# get on network for each set
for(set in c("HBMayo","Response","Remission")){
    set.rxn <- paste0("EffectSize",set)

   # edges.id <- edges[!is.na(get(set.rxn)),]
    nodes.edges <- pruneNetwork(nodes = nodes, edges = edges)

    nodes.id <- nodes.edges[["nodes"]]
    edges.id <- nodes.edges[["edges"]]
    grph <- graph_from_data_frame(edges.id[,.(from,to,rxn.id)])
    dgr.mets <-degree(grph)
    nodes.id[,degree := dgr.mets[as.character(id)]]

    # include he central measure significance
    dg.pvals <- fread(file.path("..","results", paste(set,tst.cent,"significance.csv", sep = ".")))
    dg.pvals[,padj := p.adjust(p.value, method = "BH")]
    sig.dg <- dg.pvals[padj <0.05, V1]
    set(x= nodes.id,
        i = which(nodes.id[,met.id] %in% sig.dg),
        j = set, 
        value = nodes.id[met.id %in% sig.dg,gsub("^,","",paste(get(set),tst.cent,sep = ","))]) #get(set) := ]


    # create the nodes mapping
    nodes.final <- nodes.id[,
                    .(id =  id, 
                    label = met.id,
                    title = paste0("<p><b>",met.name,"</b><br>","Changed: ",get(set),"</p>"),
                    color = mapColorToFactor(factor(get(set)), cols = "Dark2"),
                    set = get(set),
                    size = mapSizeToScale(degree))]

    edg.subs <- unique(rbind(edges.id[,.(label = to.id,subsystem)],edges.id[,.(label = from.id, subsystem)]))#[!duplicated(label),]
    edg.subs <- edg.subs[,.(subsystem = paste(gsub(",","",subsystem), collapse = ", ")),by=label]
    nodes.final <- merge(nodes.final, edg.subs, by = "label")
    nodes.final[,subsystem := gsub(",$","",paste(subsystem,gsub("Sig. in $","",paste0("Sig. in ",set)), sep=","))]
    nodes.final[set != "",subsystem := gsub(",$","",paste(subsystem,"Sig. in other layer", sep=","))]


    # adjust Eff.size from NA (not significant) to 0
    set(edges.id,
        i = which(edges.id[,is.na(get(set.rxn))]),
        j = set.rxn,
        value = 0)

    edgcols <- brewer.pal(name = "PiYG",n=9)[c(1,9)]
    edgcols <- c(edgcols[1],"gray",edgcols[2])
    edges.final <- edges.id[,.(
                         from = from,
                         to = to,
                         arrows = arrows,
                         subsystem = subsystem,
                         #label = rxn.id,
                         color = mapColorToScale(sign(get(set.rxn)),cols =edgcols ),
                         title = paste0("<p><b>",rxn.id,"<br>",reac.name,"</b><br>","Eff.Size: ",round(get(set.rxn),2),"<br>",subsystem, "</p>"))]

    # create a table for the legend
    nodes.legend <- cbind(unique(nodes.final[,.(color,label = paste0("Changed:",set))]),
                               size = 25, shape = "dot")
    nodes.legend <- rbind(nodes.legend,
                          data.table(color = brewer.pal(name = "Dark2", n = 3)[1],
                                     label = c(paste("Degree: ",range(nodes.id[,degree]))),
                                     size = range(nodes.final[,size]),
                                     shape = "dot"
                                     ))
    edge.legend <- data.table(color =  brewer.pal(name = "PiYG", n=9)[c(1,9)],
                              label = c("Negative association", "Positive association"),
                              width = 3
                              )

    vis.net <- visNetwork(nodes = nodes.final,
                          edges = edges.final,
                          main = gsub("HBMayo","HBI/Mayo score",set),
                          width = "100%", height = "1500px") %>%
        visIgraphLayout(randomSeed = 346, layout = "layout_in_circle") %>%
        visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover =FALSE, labelOnly = FALSE, hideColor = rgb(200,200,200,255*0,maxColorValue=255)), # 200,200,200,255*0.3
                   nodesIdSelection = list(enabled = TRUE, values = nodes.final[set != "", id], main = "Sig. in other layer"),
                   selectedBy = list( variable = "subsystem", multiple = TRUE, main = "Subsystem/Sig. layers")) %>%
        visLegend(addEdges = edge.legend, addNodes = nodes.legend, position = "right", ncol = 1,width = 0.1, useGroup = FALSE) %>%
        visEdges(width = 3)
    visSave(vis.net,
            file =file.path(dat.dir, paste0(set,".networkKeepAll.html")),
            selfcontained = TRUE
    )
}
