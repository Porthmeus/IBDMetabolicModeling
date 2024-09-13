##Run community optimization 
runCommunityOptimizations_Par <- function(mic.table, models, pFBAcoeff ,nr.cores=24, nutrition=NULL){        
# set up the paralelization
    require(foreach)
    require(doParallel)
    require(parallel) # to detect the number of cores
    file.create("comm_mods.log") # creating a blank file
    if(is.na(nr.cores)) nr.cores <- detectCores()-1
    cl<-makeCluster(nr.cores)
    registerDoParallel(cl) ##from doParallel, To register doParallel to be used with foreach,
    # you must call the registerDoParallel function
    out <- list()
    index <- 0

    ##Go through all communities in the subsample
    out <- foreach(i = 1:dim(mic.table)[2], .packages=c('sybil','MicrobiomeGS2','cplexAPI'),.export=c("construct_community_models_single")) %dopar% {
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
    sybil::SYBIL_SETTINGS("SOLVER_CTRL_PARM",data.frame("CPXPARAM_Threads"=1L,"CPXPARAM_LPMethod"=1L)); ok <- 1
       sybil::SYBIL_SETTINGS("METHOD", "primopt")        
        res <- list() # result object
        cat(paste0("\r Optimizing community ",i," of ",dim(mic.table)[2],"\n"), file = "comm_mods2.log", append = T)        
        modelName <- colnames(mic.table)[i] # we take each sample (community)
    model <-construct_community_models_single(mic.table, i, models, abunCutoff= 1e-3, nutrition = nutrition, modelName)

        ##Get initial solver object
        modj_warm <- sysBiolAlg(model$modj,
                    algorithm = "mtfEasyConstraint2",
                    easyConstraint=model$coupling,
                             ##easyConstraint=list(),
                    pFBAcoeff = pFBAcoeff)

        ##Get solution for basic problem
        solj1 <- optimizeProb(modj_warm)
	
	if (solj1$stat==1){
            ##Now reconstrain the model and determine biomass production part of the flux
            fluxes <- solj1$fluxes[1:model$modj@react_num]

            ## Get community growth
	    out.gr <- solj1$fluxes[grep("EX_bio1", model$modj@react_id)]

            # Get metabolic interchange
library("MicrobiomeAGORA")
            met.interchange <- MicrobiomeAGORA::get_metabolic_interchange(model$modj, solj1)
            names(solj1$fluxes) <- model$modj@react_id

            ##Get reduced costs
            redCosts <- getDjCPLEX(modj_warm@problem@oobj@env,modj_warm@problem@oobj@lp,0,model$modj@react_num-1)
            redCosts <- redCosts[grepl("^EX_",model$modj@react_id)]
            names(redCosts) <- model$modj@react_id[grepl("^EX_",model$modj@react_id)]
            redCosts <- sort(redCosts,decreasing=TRUE)

    	    outObj <- list(solj = solj1, community.growth = out.gr, met.interchange = met.interchange, models=model$model.IDs,redCosts=redCosts)
            print(paste("Growth:",out.gr))
            cat(paste0("\r Optimizing community ",i," of ",dim(mic.table)[2]," Growth:",out.gr,"\n"), file = "comm_mods2.log", append = T)
            
	    return(outObj)
	}else {
	    NULL
	}
    }

    names(out) <- colnames(mic.table)
    return(out)
}
