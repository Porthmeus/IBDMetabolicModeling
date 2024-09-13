rm(list=ls())
setwd("/home/samer/projects/intervention4/rscripts/")
library(stringr)
targets <- read.csv("/home/samer/projects/intervention4/significantMetabolites.csv")
targets <- targets[targets$dataset=="MicrobiomeGS",]
target_ids <- c(paste0(targets$cpd, "_o"))
patient_ids <- list.files("/home/samer/projects/intervention4/interv_per_comm")
dirs <- list.dirs("/home/samer/projects/intervention4/interv_per_comm")
dirs <- paste0(dirs[2:length(dirs)], "/")
interv_ids <- c()
for (dir in dirs){
  interv_ids <- union(interv_ids, list.files(dir))
}
interv_list <- list()
for (interv in interv_ids){
  interv_table <- as.data.frame(matrix(0, nrow = length(patient_ids), ncol = length(target_ids)))
  names(interv_table) <- c(target_ids)
  rownames(interv_table) <- patient_ids
  interv_list[[interv]] <- interv_table
}
target <- target_ids[1]
for (target in target_ids){
  target_table  <- read.csv(paste0("/home/samer/projects/intervention4/targets/", target, ".csv"))
  for (interv in names(interv_list)){
    if (interv %in% names(target_table)){
      ind <- match(row.names(target_table), rownames(interv_list[[interv]]))
      interv_list[[interv]][ind, target] <- target_table[, interv]
    }
  }
}

for (interv in names(interv_list)){
  current_table <- interv_list[[interv]]
  write.csv(current_table, paste0("/home/samer/projects/intervention4/intervention_tables/", interv),
              row.names = T)
}

