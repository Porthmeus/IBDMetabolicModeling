rm(list=ls())
setwd("/home/samer/projects/intervention4/rscripts/")
library(stringr)
files <- list.files("/home/samer/projects/intervention4/comm_intervention_tables/")
file_pathes <- paste0("/home/samer/projects/intervention4/comm_intervention_tables/", files)
i=1
for (i in 1:length(file_pathes)){
  print(i)
  comm_table <- read.csv(file_pathes[i], row.names = 1)
  comm_table <- comm_table-comm_table[,"orig"]
  comm_table[abs(comm_table)<10^-6] <- 0
  write.csv(comm_table, paste0("/home/samer/projects/intervention4/comm_diff_tables/", files[i]))
}


