rm(list=ls())
#usage:
#Rscript runSimulations_range.R path/models.RDS path/models_growth.csv path/OTU_abundance.table.csv first_col_num last_col_num path/result_comm_model.RDS
args <- commandArgs(trailingOnly = TRUE)


#args <- c("/zfshome/sukem124/aging/mice_baseline/metamouse-20220304.RDS", "/zfshome/sukem124/aging/mice_baseline/mice_baseline_growth.csv", "/zfshome/sukem124/aging/mice_baseline/df_MAGabundances_20210907.txt", "0.005", "11", "15", "/zfshome/sukem124/aging/mice_baseline/batch_files/11_community_models.RDS")

if (length(args)!=7){
	stop("Error: not the correct number of input arguments (should be 7)")
}

library("MicrobiomeGS2")
source("/zfshome/sukem124/Rscripts/adjustNutrition.R")
source("/zfshome/sukem124/Rscripts/construct_community_models_single.R")
source("/zfshome/sukem124/Rscripts/runCommunityOptimizations_Par.R")
library("sybil")
library("cplexAPI")

sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1 #set the solver into cplexAPI
sybil::SYBIL_SETTINGS("SOLVER_CTRL_PARM",data.frame("CPX_PARAM_THREADS"=1L)); ok <- 1 #use only one thread
models <- readRDS (args[[1]]) #load the single gapseq models
mic.table <- read.table(args[[3]], sep="\t", header=TRUE, stringsAsFactors=FALSE,row.names=1) #lo0ad the abundance (OTU) table

###################################

    ##Kick out strains with no growth in the diet
growth <- read.table(args[[2]], row.names = 1, header = T, sep = ",") #load the individual bacteria growth
rownames(growth) <- growth$Model
growth <- growth[names(models),]
growth_thr <- args[[4]] #threashold for kicking out low growing bcteria
growing_ind <- growth$Growth >= growth_thr
growing_models <- growth$Model[growing_ind]
growing_models_in_mic <- intersect(growing_models, rownames(mic.table))
mic.table <- mic.table[growing_models_in_mic,]
models=models[rownames(mic.table)]

mic.table <- mic.table[,args[[5]]:args[[6]]] #take subset of the communities

#for testing:
pFBAcoeff=1e-6 #an FBA parameter
nr.cores = 24
nutrition=NULL #if we have sample specific nutrition, we use the nutrition file
i=1
abunCutoff = 0.001
origDietFactor=0.01

performCommunityOptimizations <- function(mic.table, outFile, nutritionFile, useAbsorption=FALSE){
    # creat a new function for community optimization
    strt<-Sys.time()
    out <- list()
        nutrition <- NULL
        out <- runCommunityOptimizations_Par(mic.table, models, pFBAcoeff=1e-6,
					     nr.cores = 24, nutrition=NULL)
    saveRDS(out, file = outFile)
    cat(paste0("Elapsed time: ",round(Sys.time()-strt, digits = 2),"\n"))
}
return
performCommunityOptimizations(mic.table,
    # the function that we created earlier here
    args[[7]],  ## the results file
    NULL  ## The diet is already encoded in the models, no need to define anything here
)

