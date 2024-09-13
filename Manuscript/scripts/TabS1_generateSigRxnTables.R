# Porthmeus
# 07.09.23

require(data.table)

# directories
rsrc <- file.path("..","data","Host")
res.dir <- file.path("..","tables")

# variables to combine
tissue <- c("biopsy","blood")
types <- c("rxnExpr","PA","FVA")
cndtns <- c("HBMayo","Response","Remission")

coefs.all <- data.table()
for(tis in tissue){
    for(typ in types){
        # get the cluster file
        cls <- fread(file.path(rsrc, tis, paste0(typ,"_DBSCAN_Cluster.csv")))
        # reduce to the relevant columns
        suppressWarnings(cls[,c("V1", "cluster", "eps", "MinPts", "isseed","ClusterID", "clusterID") := NULL])
        for(cnd in cndtns){
            # get the coef file
            dr <- file.path(rsrc, tis)
            coefs <- fread( file.path(dr, list.files(dr,pattern = paste0(typ,".*\\.", cnd,".*_coefs.csv$"))))
            # renaming of the columns to harmonize namining
            cls.n <- c(grep("[t|z] value",colnames(coefs)),
                       grep("Pr\\(>\\|", colnames(coefs)))
            colnames(coefs)[cls.n] <- c("EffSize", "p.value")
            if(typ != "FVA"){
                coefs <- coefs[,.(Condition = cnd,
                                  Tissue = tis,
                                  Data.Type = typ,
                                  coef = typ,
                                  rxn, Estimate, `Std. Error`, EffSize, p.value, padj)]
            } else {
                coefs <- coefs[,.(Condition = cnd,
                                  Tissue = tis,
                                  Data.Type = typ,
                                  coef = coef,
                                  rxn, Estimate, `Std. Error`, EffSize, p.value, padj)]
            }
            # expand the clusters and correct the values
            coefs <- merge(cls, coefs, by.x = "rep.rxn", by.y = "rxn")
            if(typ != "FVA"){
                coefs[,Estimate := Estimate * sign.cor]
                coefs[,EffSize := EffSize * sign.cor]
            } else {
                coefs[coef == "center",Estimate := Estimate * sign.cor.center]
                coefs[coef == "center",EffSize := EffSize * sign.cor.center]
                coefs[coef == "range",Estimate := Estimate * sign.cor.range]
                coefs[coef == "range",EffSize := EffSize * sign.cor.range]
            }
            # combine everything into one table
            suppressWarnings(coefs[,c("sign.cor","sign.cor.range","sign.cor.center") := NULL])
            coefs.all <- rbind(coefs.all, coefs)
        }
    }
}

fwrite(coefs.all, file = file.path(res.dir, "TabS1_StatisticsHost.csv"), row.names = FALSE)
