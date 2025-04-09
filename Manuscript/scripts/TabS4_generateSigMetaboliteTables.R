# Porthmeus
# 10.02.23

require(data.table)

ddr <- file.path("..","data","Metabolomics")
fls <- list.files(ddr, pattern = "coefs.csv")
fls <- fls[-grep("crosslinks", fls)]

tab <- lapply(file.path(ddr,fls),function(x){
                  tb <- fread(x)[padj <0.05,]
                  mod.type <- strsplit(x, split = "\\.")[[1]][2+2]
                  set <- strsplit(x, split = "\\.")[[1]][3+2]
                  tb <- cbind(tb, model.type = mod.type, set = set)
                  if(mod.type == "LMM"){
                      tb <- tb[,-"df"]
                  }
                  colnames(tb)[5:6] <- c("t/z value","Pr(>|t/z|)")
                  return(tb)
})

tab <- do.call(rbind, tab)

# add the human readible names
met2recon <- fread(file.path("..","data","Metabolomics","Metabolomics2Recon3D.csv"))
met2recon[,colAbbr := gsub("(^\\d)","X\\1",gsub("[\\(|\\.|:| |\n |\\)|/|-]",".",Abbreviation))]
tab <- merge(tab, met2recon[,-c("Note","Ex_rxn")], by.x = "metabolite", by.y = "colAbbr", all.x = TRUE)


write.csv(tab, file.path("..","tables","TabS4_StatsMetabolomics.csv"), row.names = FALSE)
