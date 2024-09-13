rm(list=ls())
#Rscript /home/sukem124/Rscripts/models_growth.R path path/growth_table.csv
library("stringr")
library("sybil") 
library("cplexAPI") 
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
sybil::SYBIL_SETTINGS("SOLVER_CTRL_PARM",data.frame("CPX_PARAM_THREADS"=1L)); ok <- 1
args <- commandArgs(trailingOnly=TRUE)
#args <- c("/zfshome/sukem124/humanizedi/", "/zfshome/sukem124/humanized/growth_table.csv")
files <- list.files(args[[1]], pattern=".RDS")
files_full <- paste0(args[[1]], files)
growth <- as.data.frame(matrix(NA,length(files),2))
names(growth) <- c("Model", "Growth")
growth$Model <- str_remove(files, ".RDS")
for (i in 1:length(files)){
	print(paste0(i," ", files[[i]]))
	model <- readRDS(files_full[[i]])
	model <- changeObjFunc (model, react = "bio1")
	opt <- optimizeProb(model)
	growth[i,2] <- opt@lp_obj
}

write.csv(growth, args[[2]])
print("done!")
