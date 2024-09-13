rm(list=ls())
#usage:
args <- commandArgs(trailingOnly = TRUE)
args <- c("/home/sukem124/IBD2/emed/HRGM_models_emed.RDS", "/home/sukem124/IBD2/emed/HRGM_emed.growth", "/home/sukem124/intervention2/emed/mic_table.csv", "/home/sukem124/gapseq_models/HRGM.EX_masses", "1", "5", "/home/sukem124/intervention2/test/1_community_models.RDS")

library("MicrobiomeGS2")
source("/home/sukem124/Rscripts/adjustNutrition.R")
source("/home/sukem124/Rscripts/construct_community_models_single.R")
source("/home/sukem124/Rscripts/runCommunityOptimizations_Par.R")
source("/home/sukem124/Rscripts/runMetaboliteInterventions_Par.R")
library("sybil")
library("cplexAPI")

sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
sybil::SYBIL_SETTINGS("SOLVER_CTRL_PARM",data.frame("CPX_PARAM_THREADS"=1L)); ok <- 1
models <- readRDS (args[[1]])
mic.table <- read.table(args[[3]], sep="\t", header=TRUE, stringsAsFactors=FALSE,row.names=1)

###################################

    ##Kick out strains with no growth in the diet
growth <- read.table(args[[2]], row.names = 1, header = T, sep = ",")
rownames(growth) <- growth$Model
growth <- growth[names(models),]
growing_ind <- growth$Growth > 1e-2
growing_models <- growth$Model[growing_ind]
growing_models_in_mic <- intersect(growing_models, rownames(mic.table))
mic.table <- mic.table[growing_models_in_mic,]
models=models[rownames(mic.table)]
masses <- read.csv(args[[4]], header=TRUE, stringsAsFactors=FALSE,row.names=1)
masses$intervention <- 2*(1000/masses[metaName,"Mass"])
rownames(masses) <- paste("EX_",rownames(masses),"_e0", sep = "")
pFBAcoeff=1e-6
nr.cores = 24
nutrition = NULL
communities <- mic.table[,args[[5]]:args[[6]]]
out <- runMetaboliteInterventions_Par(communities, models, pFBAcoeff, nr.cores, masses, nutrition)

saveRDS(out, args[[7]])
