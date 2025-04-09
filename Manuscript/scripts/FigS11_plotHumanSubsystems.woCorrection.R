# Porthmeus
# 03.11.22

# plot the enriched subsystems of the different analysis

require(data.table)
require(ggplot2)
require(clusterProfiler)

# data dir
dat.dir <- file.path("..","data","Host")

# annotations
subsystems <- fread(file.path(dat.dir,"subsystems.csv"))

# create a list to store the final significant subsystems in
sig.subs.all <- list()
for(set in c("LMM.Response","LMM.Remission")){
    coefs.plt <- data.table()
    for(tissue in c("blood","biopsy")){
        # I need the cluster information to expand the estimates accordingly
        cluster.rxnExpr <- fread(file.path(dat.dir, tissue,"rxnExpr_DBSCAN_Cluster.csv"))
        cluster.PA <- fread(file.path(dat.dir, tissue,"PA_DBSCAN_Cluster.csv"))
        cluster.FVA <- fread(file.path(dat.dir, tissue,"FVA_DBSCAN_Cluster.csv"))

        # define Figure number 
        figNo <- 0
        figNo <- figNo + 1
        # get the coef table
        coefs.rxnExpr <- fread(file.path(dat.dir, tissue,paste0("rxnExpr.",set,"_coefs.csv")))
        coefs.PA <- fread(file.path(dat.dir, tissue,paste0("PA.", set, "_coefs.csv")))
        coefs.FVA <- fread(file.path(dat.dir, tissue,paste0("FVA.", set,"_coefs.csv")))

        # get the significant reactions
        sigs.rxnExpr <- coefs.rxnExpr[padj < 0.05,]
        sigs.PA <- coefs.PA[padj < 0.05,]
        sigs.FVA <- coefs.FVA[padj < 0.05,] 
        sigs <- list(rxnExpr = sigs.rxnExpr,
                     PA = sigs.PA,
                     FVA = sigs.FVA)

        # calculate hypergeometric test
        hgt <- enricher(unique(unlist(sapply(sigs, function(x) x[,rxn]))), TERM2GENE = subsystems)
    #    dotplot(hgt, x = "Count", showCategory = nrow(data.frame(hgt)))

        # do the gsea test

        # get the common columns of the coefs tabs
        clmns <- intersect(colnames(coefs.rxnExpr), colnames(coefs.PA))
        clmns <- intersect(colnames(coefs.FVA), clmns)
        # merge on 
        coefs.all <- rbind(
                           cbind(coefs.rxnExpr[,..clmns], data.set = "rxnExpr"),
                           cbind(coefs.PA[,..clmns], data.set = "PA"),
                           cbind(coefs.FVA[,..clmns], data.set = "FVA"))
        gseaTab <- coefs.all[,.(Estimate.median = median(Estimate)), by = rxn]
        gseaVector <- gseaTab[,Estimate.median]
        names(gseaVector) <- gseaTab[,rxn]
        gseaVector <- sort(gseaVector, decreasing = TRUE)
        gsea.all <- GSEA(gseaVector,
                     TERM2GENE = subsystems,
                     minGSSize = 3,
                     maxGSSize = 1E6,
                     verbose = FALSE,
                     pvalueCutoff = 0.05)


        # merge the cluster with the coef table to expand the list of reactions to the original one and correct the direction of the estimate
        coefs.rxnExpr <- merge(coefs.rxnExpr,
                               cluster.rxnExpr[,.(rxn,rep.rxn, sign.cor)],
                               by.x = "rxn",
                               by.y = "rep.rxn",
                               suffix = c("",".clustered"),
                               all.x = TRUE)
        coefs.rxnExpr[,Estimate := Estimate * sign.cor]
        coefs.PA <- merge(coefs.PA,
                               cluster.PA[,.(rxn,rep.rxn, sign.cor)],
                               by.x = "rxn",
                               by.y = "rep.rxn",
                               suffix = c("",".clustered"),
                               all.x = TRUE)
        coefs.PA[,Estimate := Estimate * sign.cor]
        coefs.FVA <- merge(coefs.FVA,
                               cluster.FVA[,.(rxn,rep.rxn, sign.cor.center, sign.cor.range)],
                               by.x = "rxn",
                               by.y = "rep.rxn",
                               suffix = c("",".clustered"),
                               all.x = TRUE)
        coefs.FVA[coef == "center",Estimate := Estimate * sign.cor.center]
        coefs.FVA[coef == "range",Estimate := Estimate * sign.cor.range]

        # combine the results again
        # get the common columns of the coefs tabs
        clmns <- intersect(colnames(coefs.rxnExpr), colnames(coefs.PA))
        clmns <- intersect(colnames(coefs.FVA), clmns)
        # merge on 
        coefs.all <- rbind(
                           cbind(coefs.rxnExpr[,..clmns], data.set = "rxnExpr"),
                           cbind(coefs.PA[,..clmns], data.set = "PA"),
                           cbind(coefs.FVA[,..clmns], data.set = "FVA"))

        # create the plots from  subsystems
        set.name <- strsplit(set, split = "\\.")[[1]][2]
        sig.subs<- unique(c(data.frame(hgt)[,"ID"], data.frame(gsea.all)[,"ID"]))
        sig.subs.all[[paste(set.name, tissue, sep = ".")]] <- sig.subs
        coefs.all.subs <- merge(coefs.all, subsystems, by.x = "rxn.clustered", by.y = "reaction", keep.x = TRUE)
        # scale the estimate for the GLMMs
        ttl <- paste("Host", set.name)
        xlab <- "sqrt(Estimate)"

        # remove 0 and INF odds ratios from the list, because these are awkward to plot
        coefs.all.subs  <- coefs.all.subs[exp(Estimate) !=0 & exp(Estimate) < Inf,]

        # combine blood and biopsy data sets
        coefs.plt <-rbind(coefs.plt,cbind(coefs.all.subs[subsystem %in% sig.subs & padj <0.05,.(Estimate, subsystem)],tissue = tissue))

    }
        sig.subs <- coefs.plt[!is.na(Estimate), unique(subsystem)]
        coefs.plt[is.na(Estimate), subsystem :=sig.subs[1]]
        coefs.plt <- unique(coefs.plt)
        sub.order <- coefs.plt[,.(od = median(Estimate)), by = subsystem]
        coefs.plt[,subsystem := factor(subsystem, level = sub.order[order(od), subsystem])]
        
        x_lim <- range(coefs.plt[,sqrt(abs(Estimate))*sign(Estimate)])
        if(x_lim[1] < -5){
            x_lim[1] <- -5
        }
        if(x_lim[2] > 5){
            x_lim[2] <-5
        }


        plt.sub.est.sig <- ggplot(coefs.plt, aes(x= sqrt(abs(Estimate))*sign(Estimate), y = subsystem)) +
            geom_vline(xintercept = 0, color = "red", linetype = 2)+
            geom_boxplot(fill = NA) +
            ggtitle(ttl)+
            facet_grid(~tissue, scale = "free_x")+
            xlim(x_lim[1],x_lim[2])+
            xlab(xlab)+
            ylab("")+
            guides(y.sec = "axis", y = "none") + # move labels to the right
            theme_bw()

        #browser()
        ggsave(plt.sub.est.sig, file = file.path("..","figures_raw",paste0("FigS11_HostSubsystemPlots",set,".woCorrection.pdf")),
               width = 8, height = 1+0.12*length(sig.subs))

}

    # save the enriched subsystems for later usage
    sig.subs.df <- do.call(rbind,
                           lapply(names(sig.subs.all), function(x) {
                                            if (length(sig.subs.all[[x]]) > 0){
                                               return(cbind(set = gsub("\\..*$","",x),
                                                 tissue = gsub("^.*\\.","",x),
                                                 subsystem = sig.subs.all[[x]]))
                               }}))
    write.csv(sig.subs.df, file = file.path(dat.dir, "EnrichedSubsystems.woCorrection.csv"), row.names = FALSE)
