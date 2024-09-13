# Porthmeus
# 31.08.24

require(data.table)
require(ggplot2)
require(cowplot)
require(stringr)

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

# get the metadata
clinic <- fread("../../resources/ClinicalDataMergeLongImputeCor.csv")
idconv <- fread("../../resources/IdentifierConversion.csv")
metabolites <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_metabolites.tsv")
rxnsAnno <- fread("https://raw.githubusercontent.com/jotech/gapseq/master/dat/seed_reactions.tsv")

plt.lst <- list()
dat.sel <- dat.all[padj <0.05 & covariate %in% c("d0","d14","d14d0Diff"),]
for(i in 1:nrow(dat.sel)){
    fl <- file.path("..","..","resources", paste0("future.eMed.", tolower(dat.sel[i, method]),"_", dat.sel[i,flux],".csv"))
    fl <- gsub("withinFlux","withinOflux", fl)
    fl <- gsub(".microbiomegs_","_",fl)
    fluxes <- fread(fl)
    colnames(fluxes)[1] <-"sample"
    fluxes.melt <- melt(fluxes,
                         id.vars = "sample",
                         variable.name = "rxn",
                         value.name = "flux")
    ## rename the samples to match the expected form
    fluxes.melt[, SeqID := gsub("^X","",sample)]
    fluxes.melt <- merge(fluxes.melt, idconv, by.x ="SeqID", by.y = "sample")
    fluxes.melt[, PatientID_Time := paste(Cohort, Person.Numbercode, Time_seq, sep = "_")]
    fluxes.merged <- merge(fluxes.melt, clinic, by = "PatientID_Time")
    fluxes.merged <- fluxes.merged[rxn == dat.sel[i,rxn], .(PatientID, Time_seq.x, flux, seqtype, source, Response=Responder, Remission=Remision)]
    # create a variable which is easily readible
    fluxes.merged[,plt.value := flux]
    if(dat.sel[i, covariate] == "d0"){
        fluxes.merged <- fluxes.merged[Time_seq.x == 0,]
        y.lab <- "flux at 0d"
    } else if(dat.sel[i, covariate] == "d14"){
        fluxes.merged <- fluxes.merged[Time_seq.x == 14,]
        y.lab <- "flux at 14d"
    } else if(dat.sel[i, covariate] == "d14d0Diff"){
        fluxes.merged <- fluxes.merged[Time_seq.x %in% c(0,14),]
        fluxes.merged <- dcast(data = fluxes.merged, PatientID+seqtype+source+Response+Remission~paste0("d",Time_seq.x), value.var = "flux")
        fluxes.merged[,plt.value := d14-d0]
        fluxes.merged <- fluxes.merged[!is.na(plt.value),]
        y.lab <- "flux diff. 14d-0d"
    }
    # get the reaction name
    if(grepl("^EX",dat.sel[i, rxn])){
           met <- gsub("^EX_(cpd\\d*).*","\\1",dat.sel[i, rxn])
           met.name <-metabolites[id == met, name]
           rxn.name <- paste0("Export of ", tolower(met.name))
    } else {
        rxn.id <- gsub(".*(rxn\\d*).*","\\1", dat.sel[i,rxn])
        rxn.name <- rxnsAnno[id == rxn.id, name]
    }
    #
    ttl <- str_wrap(dat.sel[i,paste(rxn.name, method, sep=", ")], 25)
    plt <- ggplot(fluxes.merged[get(dat.sel[i,outcome]) != "C",], aes(x=get(dat.sel[i,outcome]), y = plt.value)) +
        geom_boxplot()+
        labs(x = dat.sel[i,outcome], y = y.lab, subtitle = ttl)+
        theme_bw()
    plt.lst[[ttl]] <- plt
}

plt <- plot_grid(plotlist = plt.lst, labels = "AUTO")
ggsave(plt, file = "../results_summary/singleTimePointPlot.pdf", width = 7, height = 3.5)


