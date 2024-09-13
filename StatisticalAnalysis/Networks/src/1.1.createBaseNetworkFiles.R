# Porthmeus
# 07.03.23

require(data.table)
require(igraph)
require(parallel)
require(doMC)
require(foreach)
require(visNetwork)

# define data dir
dat.dir <- file.path("..","results")
rsrc.dir <- file.path("..","resources")
sets <- c("HBMayo","Response","Remission")

# get the basic information on rxn2metabolite relationship
bas <- fread(file.path(rsrc.dir,"Host","Mets2Subsystems.csv"))
#rm.mets <- fread(file.path(rsrc.dir,"Host","cofactorsAnorganicCompounds.csv"))
rm.mets <- c("h","h2o","o2","pi","co2")
bas[,met.simple := gsub("\\[[a-z]\\]","",metabolite)]
bas <- bas[!met.simple %in% rm.mets,]
bas <- bas[!grepl("^Transport|^Exchange", subsystem),]


edg <- fread(file.path(dat.dir, "baseEdges.csv"))
# remove transport reactions and compartments
edg <- edg[!grepl("^Transport|^Exchange", subsystem),]
edg <- unique(edg[,c("from.id","to.id") := list(gsub("\\[[a-z]\\]$","",from.id),
                                                            gsub("\\[[a-z]\\]$","",to.id))])
colnames(edg)[3] <- "rxn.id"
edg <- edg[!from.id %in% rm.mets,]
edg <- edg[!to.id %in% rm.mets,]

# add back the information on the reaction names
edg <- merge(edg, unique(bas[,.( reac.name, reaction)]), by.x = "rxn.id", by.y = "reaction", all.x = TRUE)

# write basic edge.table
fwrite(edg, file.path(dat.dir, "basic.edge.table.1.all.csv"))

# add information on sig reactions

for(tissue in c("blood","biopsy")){
    coefs_all.host <- data.table()
    lst.dir <- file.path(rsrc.dir, "Host", tissue)
    for(st in sets){
        fls <- list.files(lst.dir, pattern = paste0(st, ".*_coefs.csv"))
        for(f in fls){
            coefs <- fread(file.path(lst.dir, f))
            if(st == "HBMayo"){
                coefs <- coefs[padj < 0.05,.(rxn, coef,
                                  EffSize = `t value`,
                                  set = st,
                                  tissue = tissue)]
            } else {
                coefs <- coefs[padj < 0.05,.(rxn, coef,
                                  EffSize = `z value`,
                                  set = st,
                                  tissue = tissue)]
            }
            # expand cluster and correct the sign of the effect size with the
            # sign of the correlation within the cluster representative
            f.dat <- strsplit(f, split="\\.")[[1]][1]
            f.clt <- fread(file.path(lst.dir,paste0(f.dat, "_DBSCAN_Cluster.csv")))
            if(f.dat != "FVA"){
                coefs <- merge(f.clt[,.(rep.rxn, rxn, sign.cor)], coefs, by.x = "rep.rxn", by.y = "rxn")
                coefs[,EffSize := EffSize * sign.cor]
            } else {
                coefs <- merge(f.clt[,.(rep.rxn, rxn, sign.cor.center, sign.cor.range)], coefs, by.x = "rep.rxn", by.y = "rxn")
                coefs[coef == "center",EffSize := EffSize * sign.cor.center]
                coefs[coef == "range",EffSize := EffSize * sign.cor.range]
            }
            # remove correlation columns
            cor.col.id <- grep("sign.cor", colnames(coefs), value = TRUE)
            coefs[, (cor.col.id) := NULL]
            coefs_all.host <- rbind(coefs_all.host, coefs)
        }
    }
    # get the sum of the EffSizes for each subset
    coefs_all.host <- coefs_all.host[,.(EffSize = mean(EffSize)), by = .(rxn,set,tissue)]

# add the information of the metabolomics data
lst.dir <- file.path(rsrc.dir, "Metabolomics")
coefs_all.metabolome <- data.table()
for(st in sets){
    fls <- list.files(lst.dir, pattern = paste0(st, ".*_coefs.csv"))
    for(f in fls){
        coefs <- fread(file.path(lst.dir, f))
        if(st == "HBMayo"){
            coefs <- coefs[padj < 0.05,.(metabolite,
                              EffSize = `t value`,
                              set = st,
                              src = "Metabolomics")]
        } else {
            coefs <- coefs[padj < 0.05,.(metabolite,
                              EffSize = `z value`,
                              set = st,
                              src = "Metabolomics")]
        }
        coefs_all.metabolome <- rbind(coefs_all.metabolome, coefs)
    }
}
# decluster the coefs
cluster.metabolome <- fread(file.path(lst.dir, "metabolites_cluster.csv"))
coefs_all.metabolome <- merge(cluster.metabolome[,.(rep.met, metabolite = met, sign.cor)], coefs_all.metabolome, by.x = "rep.met", by.y = "metabolite")
# correct EffSize to the sign of the correlation in the cluster
coefs_all.metabolome[,EffSize := EffSize* sign.cor]
# translate to met.ids
met2vmh <- fread(file.path(lst.dir, "Metabolomics2Recon3D.csv"))
coefs_all.metabolome <- merge(coefs_all.metabolome, met2vmh[,.(eukratis = make.names(Abbreviation), vmh = Recon3D)], by.x = "metabolite", by.y = "eukratis")

## add the metabolites from the bacteria simulations
lst.dir <- file.path(rsrc.dir, "Microbiome")
coefs_all.microbiome <- data.table()
for(sim in c("MicrobiomeGS","Bacarena")){
    for(st in sets){
        fls <- list.files(file.path(lst.dir,sim), pattern = paste0("outflux.*",st, ".*_coefs.csv"))
        for(f in fls){
            coefs <- fread(file.path(lst.dir,sim, f))
            if(st == "HBMayo"){
                coefs <- coefs[padj < 0.05,.(metabolite=cpd,
                                  EffSize = `t value`,
                                  set = st,
                                  sim = sim)]
            } else {
                coefs <- coefs[padj < 0.05,.(metabolite=cpd,
                                  EffSize = `z value`,
                                  set = st,
                                  sim = sim)]
            }
            coefs[,metabolite := gsub("^.*_(cpd\\d*)","\\1", metabolite)]
            # expand cluster
            f.clt <- fread(file.path(lst.dir,sim,"outflux_cluster.csv"))
            f.clt[,rep.met := gsub("^.*_(cpd\\d*).*$","\\1", rep.rxn)]
            f.clt[,met := gsub("^.*_(cpd\\d*).*$","\\1", rxn)]
            coefs <- merge(f.clt[,.(rep.met, metabolite = met), sign.cor], coefs, by.x = "rep.met", by.y = "metabolite")
            # correct EffSize sign with correlation sign in the cluster
            coefs[,EffSize := EffSize*sign.cor]
            coefs_all.microbiome <- rbind(coefs_all.microbiome, coefs)
        }
    }
}
coefs_all.microbiome[,metabolite := gsub("^.*_(cpd\\d*)","\\1", metabolite)]
coefs_all.microbiome[,src := "Microbiome"]
# translate the metabolites
met2vmh <- fread(file.path(lst.dir, "seed2vmh.csv"))
coefs_all.microbiome <- merge(coefs_all.microbiome, met2vmh[,.(seed =M2_SEED, vmh = M1_VMH)], by.x = "metabolite", by.y = "seed")

# create networks specific to each condition that was tested
for(st in sets){
    coefs.host <- coefs_all.host[set == st,]
    # create a common vector for biopsy and blood
    coefs.host.melt <- coefs.host[,.(EffSize = mean(EffSize, na.rm = TRUE)), by = .(rxn, set, tissue)]
#    coefs.host <-merge(coefs.host.melt,coefs.host, by = "rxn")
    edg.set <- merge(edg, coefs.host.melt, by.x = "rxn.id", by.y= "rxn")

    # create a nodes data frame
    nodes.id <- data.table(met.id = unique(c(edg.set[,from.id], edg.set[,to.id])))
    nodes.id[,id := 1:nrow(nodes.id)]
    nodes.id <- merge(nodes.id, unique(bas[,.(metabolite=gsub("\\[[a-z]\\]","",metabolite), met.name)]),by.x = "met.id",by.y = "metabolite", all.x = TRUE)
    # add the indication for blood and biopsy metabolites
    bas.set <- merge(bas, coefs.host.melt, by.x = "reaction",by.y = "rxn")
    bas.set <- bas.set[,.(reaction,  subsystem = paste(unique(gsub(",",";",subsystem)), collapse = ","), met.name), by = met.simple]
    nodes.id <- merge(nodes.id, unique(bas.set[!is.na(tissue),.(met.simple, tissue, subsystem)]), by.x = "met.id",by.y = "met.simple", all.x = TRUE)

    # add the metabolomics and the microbiome data
    coefs.microbiome <- coefs_all.microbiome[set == st,]
    coefs.metabolome <- coefs_all.metabolome[set == st,]
    coefs.metabolome[,sim:= "Metabolomics"] # just to keep some information, but still to simplify later
    clmns <- intersect(colnames(coefs.metabolome), colnames(coefs.microbiome))
    coefs.met.mic <- rbind(coefs.metabolome[,..clmns], coefs.microbiome[,..clmns])
    coefs.met.mic.melt <- 
        coefs.met.mic[,.(EffSize = mean(EffSize, na.rm = TRUE),
                  src = paste(unique(src), collapse = ";")), by = vmh]
    coefs.cast <- dcast(coefs.met.mic, vmh ~ paste0("ES.",sim), value.var = "EffSize", fun.aggregate = mean)
    coefs.met.mic <- merge(coefs.met.mic.melt, coefs.cast)
    # check if all casts are present
    col.checks <- c("ES.Bacarena","ES.Metabolomics","ES.MicrobiomeGS")
    col.checks <- col.checks[!(col.checks %in% colnames(coefs.met.mic))]
    if(length(col.checks) >0){
        for(i in col.checks){
            coefs.met.mic[[i]] <- NA
        }
    }
    nodes.id <- merge(nodes.id, coefs.met.mic, by.x = "met.id", by.y = "vmh", all.x = TRUE)
    
    # write the tables to disk
    write.csv(nodes.id,
              file = file.path(dat.dir, paste0("basic.nodes.table.1.",st,".",tissue,".csv")),
              row.names = FALSE)
    write.csv(edg.set,
              file = file.path(dat.dir, paste0("basic.edge.table.1.",st,".",tissue, ".csv")),
              row.names = FALSE)
}
}
