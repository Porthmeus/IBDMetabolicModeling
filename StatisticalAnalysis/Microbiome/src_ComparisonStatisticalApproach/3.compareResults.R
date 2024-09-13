# Porthmeus
# 08.07.24

require(data.table)
require(ggplot2)
require(ggvenn)
require(cowplot)

# define sets
sets <- c("Bacarena","MicrobiomeGS")

# small helper function
readData <- function(file.description, file.list, directory){
    dat <- fread(file.path(directory, file.list[file.description]))
    extra_info <- strsplit(file.description, split = "\\.")[[1]]
    flux <- extra_info[1]
    outcome <- extra_info[2]
    if(length(extra_info) > 2){
        covariate <- extra_info[3]
    } else {
        covariate <- "None"
    }
    dat <- dat[,c("flux","outcome","covariate","set") := list(flux, outcome,covariate,file.description)]
    return(dat)
}

dat.all <- data.table()
for(st in sets){
# load data

    set.dr <- file.path("..", paste0(st,"_therapy"),"src")
    res.dir <- file.path(set.dr, "..","results")
    res.dir.old <- file.path(set.dr,"..","..",st,"results")

    fls <- list.files(res.dir, "(Response|Remission).*_coefs.csv")
    names(fls) <- gsub("_coefs.csv","",fls)
    fls.old <- list.files(res.dir.old, "(Response|Remission).*_coefs.csv")
    names(fls.old) <- gsub("intervention","HBMayoCor",gsub("(GLMM\\.|_coefs.csv)","",fls.old))

    dat.new <- lapply(names(fls), readData, file.list = fls, directory = res.dir)
    dat.old <- lapply(names(fls.old), readData, file.list= fls.old, directory = res.dir.old)
    dat.temp <- append(dat.new, dat.old)
    dat.temp[["fill"]] = TRUE
    dat.temp <- do.call(rbind,dat.temp)
    dat.temp[,method := st]
    dat.all <- rbind(dat.all, dat.temp, fill = TRUE)
}
dat.all<- dat.all[!is.na(rxn),]

dat.bar <- dat.all[,.(value = sum(padj<0.05,na.rm=TRUE)), by = .(outcome, covariate,method)]

# get the barplot
plt.bar <- ggplot(dat.bar, aes(x=covariate, y = value)) +
    geom_bar(stat = "identity",aes(fill = method), color = "black")+
    facet_grid(~outcome)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust = 1, vjust = 1))+
    labs(title = "Microbiome", x = "covariate/subset/interaction", y = "sig. rxns", fill = "Method")
ggsave(plt.bar, file = file.path("..","results_summary","SigRxns_Barplot.pdf"), width =4.5, height=4)

# get the venn diagram


lst.response <- list()
lst.response[["HBMayoCor"]] <- dat.all[outcome == "Response" & covariate == "HBMayoCor" & padj <0.05, unique(rxn)]
lst.response[["None"]] <- dat.all[outcome == "Response" & covariate == "None" & padj <0.05, unique(rxn)]
venn.response <- ggvenn(lst.response) + ggtitle("Response")

lst.remission <- list()
lst.remission[["HBMayoCor"]] <- dat.all[outcome == "Remission" & covariate == "HBMayoCor" & padj <0.05, unique(rxn)]
lst.remission[["None"]] <- dat.all[outcome == "Remission" & covariate == "None" & padj <0.05, unique(rxn)]
venn.remission <- ggvenn(lst.remission)+ggtitle("Remission")

venn.diagram <- plot_grid(venn.response,venn.remission)
ggsave(venn.diagram, file = file.path("..","results_summary","VennPlot_sigRxns.pdf"), width = 12.5, height = 4.5)

