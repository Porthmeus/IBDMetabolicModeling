# Porthmeus
# 23.08.22

require(data.table)
require(ggplot2)
require(ComplexUpset)
require(RColorBrewer)

dat.dir <- file.path("..","data","Microbiome")

# Annotations
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites_edited.tsv")
setkey(metabolites, "id")

coefs_summary <- list()
figno<- 0
for(set in c("LMM.HBMayo.etiology","GLMM.Response.intervention","GLMM.Remission.intervention")){
    figno <- figno+1
    # coefficients
    ba.coefs <- fread(file.path(dat.dir, paste0("Bacarena/outflux.", set, "_coefs.csv")))
    mgs.coefs <- fread(file.path(dat.dir, paste0("MicrobiomeGS/outflux.", set, "_coefs.csv")))

    # expand the set by the clustered reactions
    ba.cluster <- fread(file.path(dat.dir, "Bacarena","outflux_cluster.csv"))
    ba.coefs.exp <- merge(ba.coefs, ba.cluster[,.(rxn, rep.rxn,sign.cor)], by.x= "rxn", by.y = "rep.rxn", suffix = c(".rep",".all"))
    # correct the estimates and t/z-values for the expanded clusters
    ba.coefs.exp[,Estimate := Estimate * sign.cor]
    # find the effectsize column and rename 
    eff.col <- grep("z value|t value", colnames(ba.coefs.exp))
    colnames(ba.coefs.exp)[eff.col] <- "Effect.size"
    ba.coefs.exp[,Effect.size := Effect.size * sign.cor]
    ba.coefs.exp[,cpd := gsub(".*_(cpd[0-9]*)","\\1",rxn.all)]
    ba.coefs.exp[,name:= metabolites[cpd,name]]

    mgs.cluster <- fread(file.path(dat.dir, "MicrobiomeGS","outflux_cluster.csv"))
    mgs.coefs.exp <- merge(mgs.coefs, mgs.cluster[,.(rxn, rep.rxn,sign.cor)], by.x= "rxn", by.y = "rep.rxn", suffix = c(".rep",".all"))
    # correct the estimates and t-values for the expanded clusters
    mgs.coefs.exp[,Estimate := Estimate * sign.cor]
    # find the effectsize column and rename 
    eff.col <- grep("z value|t value", colnames(mgs.coefs.exp))
    colnames(mgs.coefs.exp)[eff.col] <- "Effect.size"
    mgs.coefs.exp[,Effect.size := Effect.size * sign.cor]
    mgs.coefs.exp[,cpd := gsub(".*_(cpd[0-9]*)_e0","\\1",rxn.all)]
    mgs.coefs.exp[,name:= metabolites[cpd,name]]

    # combine data sets in one table
    cols<- union(colnames(mgs.coefs.exp), colnames(ba.coefs.exp))
    coefs.exp <- data.table(rbind(
                                  cbind(mgs.coefs.exp[,..cols], dataset = "MicrobiomeGS"),
                                  cbind(ba.coefs.exp[,..cols], dataset = "Bacarena")
                                  ))

    # correct the metabolite names
    coefs.exp[is.na(name), name := cpd]
    coefs.exp[,name := tools::toTitleCase(gsub(" \\(.*\\)","",name))]
    coefs.exp <- coefs.exp[padj <0.05,]
    coefs.sort <- coefs.exp[,.(srt = mean(Estimate)), by = name]
    coefs.exp[,name:= factor(name, levels = coefs.sort[order(srt),name])]

    # read raw data
    mgs.flux <- read.csv(file.path(dat.dir, "MicrobiomeGS", "future.eMed_outflux.csv"), row.names = 1)
    ba.flux <- read.csv(file.path(dat.dir, "Bacarena", "future.eMed.bacarena_outflux.csv"), row.names = 1)
    mgs.flux.melt <- reshape2::melt(as.matrix(mgs.flux), varnames = c("Sample","Rxn"))
    ba.flux.melt <- reshape2::melt(as.matrix(ba.flux), varnames = c("Sample","Rxn"))
    flux.raw <- data.table(rbind(mgs.flux.melt, ba.flux.melt))

    flux.drctn <- flux.raw[,.(production = mean(value)>0) ,by=.(Rxn)]
    flux.drctn <- flux.drctn[production==T, by.microbiome := "produced"]
    flux.drctn <- flux.drctn[!production==T, by.microbiome := "consumed"]
    flux.drctn <- flux.drctn[,by.microbiome := factor(by.microbiome, levels = c("produced","consumed"))]
    coefs.drctn <- merge(coefs.exp, flux.drctn, by.x = "rxn", by.y = "Rxn")
    
    # scale the estimates for the GLMMs
    if(grepl("GLMM",set)){
        coefs.drctn[, Estimate := sign(Estimate)*sqrt(abs(Estimate))]
        coefs.drctn[, `Std. Error` := sign(`Std. Error`)*sqrt(abs(`Std. Error`))]
        xlab = "sqrt(Estimate)"
        title = paste("Microbiome", strsplit(set,split = "\\.")[[1]][2])
    } else { 
        xlab = "Estimate"
        title = paste("Microbiome","HBI/Mayo score")
    }
    
    # save the data for later usage
    coefs_summary[[set]] <- cbind(coefs.exp, testcondition = title)
    

    pl <- ggplot(coefs.drctn[padj <= 0.05,],aes(x=Estimate, y = name, color = dataset, shape = by.microbiome)) +
        geom_vline(xintercept = 0, linetype = 2, color = "red")+
        geom_point()+#(size = 3) +
        scale_color_brewer(palette = "Paired") +
        geom_errorbar(aes(xmin = Estimate - `Std. Error`, xmax = Estimate + `Std. Error`))+
        ylab("") +
        xlab(xlab)+
        ggtitle(title)+
        theme_bw()
    ggsave(pl, file = paste0(file.path("..","figures_raw","/Fig1"),LETTERS[figno],"_MicrobialMetabolites.", set, ".pdf"), height = 0.12*length(unique(coefs.drctn[padj <= 0.05, name])) +1 , width = 4)
}

coefs_summary <- lapply(coefs_summary, function(x) {
                            if("df" %in% colnames(x)){
                                x<- x[,-"df"]
                            }
                            colnames(x)[6:7] <- c("Effect.size","p.value")
       return(x)})

coefs_summary <- do.call(rbind, coefs_summary)
write.csv(file = file.path("..","temp","significantMetabolites.csv"),
          coefs_summary,
          row.names = FALSE)

