rm(list=ls())
setwd("/home/samer/projects/intervention4/rscripts/")
library(stringr)
pheno_flux <- read.csv("/home/samer/projects/IBD3/IBDClinicalFluxFinal.csv", sep = "\t")
rownames(pheno_flux) <- pheno_flux$sample
cpd_ind <- grep("EX_cpd", names(pheno_flux))
cpds <- names(pheno_flux)[cpd_ind]
cpds <- str_remove_all(cpds, "EX_|_e0") ##all metabolic flux ids (within and between)
dirs <- list.dirs("/home/samer/projects/intervention4/interv_per_comm")
ids <- list.files("/home/samer/projects/intervention4/interv_per_comm")
dirs <- paste0(dirs[2:length(dirs)], "/")
interv_ids <- c() ## all interventions
for (dir in dirs){
  interv_ids <- union(interv_ids, list.files(dir[1]))
}
target_list <- list() ## a list of dataframes (intervention x sample) for each target
interv_ids <- str_remove(interv_ids, ".csv")
for(target in cpds){
  df <- as.data.frame(matrix(0,length(ids), length(interv_ids)))
  names(df) <- interv_ids
  rownames(df) <- ids
  target_list[[target]] <- df
}
#i=1
for(i in (1:length(dirs))){ ## for each sample
  dir <- dirs[i]
  patient_id <- ids[i]
  print(paste0(i,"_",patient_id))
#  file <- list.files(dir)[100]
  for(file in list.files(dir)){
      interv_id <- str_remove(file, ".csv")
      file <- paste0(dir, file)
      intervent_table <- read.csv(file)
      intervent_table$rxn <- str_remove_all(intervent_table$rxn, "EX_|_e0")
      flux_table <- intervent_table[,c("rxn", "flux")]
      fluxo_table <- intervent_table[,c("rxn", "o.flux")]
      names(fluxo_table) <- c("rxn", "flux")
      fluxo_table$rxn <- paste0(fluxo_table$rxn, "_o")
      intervent_table2 <- rbind(flux_table, fluxo_table)
#      target <- cpds[360]
      for (target in cpds){
        ind <- which(intervent_table2$rxn==target)
        if(length(ind)==0){
          flux=0
        }
        else{
          flux <- intervent_table2$flux[ind]
        }
        target_list[[target]][patient_id, interv_id] <- flux
      }
    }
}
for(target in names(target_list)){
  target_table <- target_list[[target]]
  write.csv(target_table, paste0("/home/samer/projects/intervention4/targets/", target, ".csv"),
            row.names = TRUE)
}
