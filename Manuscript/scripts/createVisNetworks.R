# Porthmeus
# 07.03.23

# This file contains functions to handle the network data better and create html files for visualizations of networks

require(RColorBrewer)
require(data.table)
require(scales)

mapColorToFactor <- function(x, cols = "Paired"){
    # maps a colors to a vector of entries
    x <- factor(x)
    # check if colors are part of ColorBrewer, else use the values given
    if(length(cols) == 1){
        if(cols %in% rownames(brewer.pal.info)){
            n.cols <- brewer.pal.info[cols, "maxcolors"]
            cols <- brewer.pal(n = n.cols, name = cols)
        } 
    } else {
        n.cols <- length(cols)

    }
    n.facs <- length(unique(x))
    if(n.facs > n.cols){
        warning(paste("The number of unique entries in your vector (",nfacs,") is larger than the number of possible values in your palette (",n.cols,"). Will recycle colors"))
        cols <- rep(cols, ceiling(n.facs/n.cols))
    }
    cols <- cols[1:n.facs]
    names(cols) <- sort(unique(x))
    return(cols[x])
}

mapColorToScale <- function(x, cols = c("blue","red"), mid = NULL){
    # maps a color gradient to values of x
    # separate values in mid and combine afterwards if mid!=NULL
    if(is.null(mid)){
        cols <- gradient_n_pal(colours = cols, values = range(x))(x)
    } else {
        cols <- gradient_n_pal(colours = cols, values = c(min(x),mid, max(x)))(x)
    }
    return(cols)
}
#        norm_x <- (x - min(x))/(max(x)-min(x))
#    } else {
#        i_l <- which(x<=mid)
#        i_h <- which(x>mid)
#        x_l <- x[i_l]
#        x_h <- x[i_h]
#        x_ld <- max(x_l-mid)
#        x_hd <- min(x_h-mid)
#        norm_x_l <- (x_l - min(x_l))/(max(x_l)-min(x_l))
#        norm_x_h <- (x_h - min(x_h))/(max(x_h)-min(x_h))
#        norm_x <- x
#        norm_x[i_l] <- c(norm_x_l+x_ld-1)
#        norm_x[i_h] <- c(norm_x_h+x_hd)
#        norm_x <- (norm_x - min(norm_x))/(max(norm_x)-min(norm_x))
#
#plot(1:10, col = colorRamp(c("blue","gray","red"),bias = 1)(10))
#
#
#
#    cols <- rgb(colorRamp(cols)(norm_x), maxColorValue = 255)
#    return(cols)
#}

mapShapeToFactor <- function(x, shapes = "all"){
    # shapes must be either "all","inside","outside" or a vector of available shapes given by visNetwork - see visNodes
    
    outside <- c("dot", "diamond",  "star", "triangle", "triangleDown", "hexagon", "square")
    inside <- c("ellipse", "circle", "database", "box", "text")

    x <- factor(x)
    n.facs <- length(levels(x))
    if(length(shapes) == 1){
        if(tolower(shapes) == "all"){
            shapes <- c(inside,outside)
        } else if(tolower(shapes) == "inside"){
            shapes <- inside
        } else if(tolower(shapes) == "outside"){
            shapes <- outside
        }
    }
    n.shapes <- length(shapes)
    if(n.facs > n.shapes){
        warning(paste("The number of unique entries in your vector (",nfacs,") is larger than the number of possible values in your palette (",n.cols,"). Will recycle colors"))
        shapes <- rep(shapes, ceiling(n.facs/n.shapes))
    }
    shapes <- shapes[1:n.facs]
    names(shapes) <- sort(unique(x))
    return(shapes[x])
}

mapSizeToScale <- function(x){
    # rescales everything between 5-45
    x <- (x-min(x))/(max(x)-min(x))
    x <- round((x*40)+5)
    return(x)
}

pruneNetwork <- function(nodes, edges){
    # takes a nodes and an edge data table and returns the same data.tables, but this time only with those entries which are present in the other one and added id, from and to columns
    nodes <- data.table(nodes)
    edges <- data.table(edges)
    nodes_in_edges <- unique(c(edges[,from.id],edges[,to.id]))
    nodes.id <- nodes[met.id %in% nodes_in_edges,]
    nodes.id[,id:= 1:nrow(nodes.id)]
    edges.id <- merge(edges, nodes.id[,.(met.id,id)], by.y = "met.id", by.x = "to.id", all.x = TRUE)
    edges.id <- merge(edges.id, nodes.id[,.(met.id,id)], by.y = "met.id", by.x = "from.id", suffix =c(".to_index",".from_index"), all.x = TRUE)
    colnames(edges.id)[grep("_index", colnames(edges.id))] <- c("to","from")
    edges.id <- edges.id[!is.na(to) & !is.na(from),]
    return(list(nodes = nodes.id, edges = edges.id))
}


