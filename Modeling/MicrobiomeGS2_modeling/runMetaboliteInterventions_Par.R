##Run metabolite interventions
runMetaboliteInterventions_Par <- function(communities, models, pFBAcoeff, nr.cores = NA, masses, nutrition = NULL){
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
    ##Go through all community communities
    out <- foreach(i = 1:dim(communities)[2],.packages=c('sybil','MicrobiomeGS2','cplexAPI'),.export=c("construct_community_models_single")) %dopar% {
  sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1
    sybil::SYBIL_SETTINGS("SOLVER_CTRL_PARM",data.frame("CPXPARAM_Threads"=1L,"CPXPARAM_LPMethod"=1L)); ok <- 1
       sybil::SYBIL_SETTINGS("METHOD", "primopt")        
       #sybil::SYBIL_SETTINGS("SOLVER_CTRL_PARM",data.frame("CPX_PARAM_THREADS"=1L)); ok <- 1   
        res <- list()
        cat(paste0("\r Optimizing community ",i," of ",dim(communities)[2],"\n"), file = "comm_mods2.log", append = T) 
    modelName <- colnames(communities)[i]
    model <- construct_community_models_single(communities,i, models, abunCutoff = 0.001, nutrition)
    index <- index + 1
    ##Get indices of exchange reactions
    exIndices <- which(grepl("^EX_",model$modj@react_id))

    ##Get indices of internal exchange reactions

    ##Get original lower bounds
    exLBs <- model$modj@lowbnd[exIndices]
    exUBs <- model$modj@uppbnd[exIndices]

    ##Get initial solver object
    modj_warm <- sysBiolAlg(model$modj, algorithm = "mtfEasyConstraint2", easyConstraint=model$coupling, pFBAcoeff = pFBAcoeff)

    ##Get solution for basic problem
        solBase <- optimizeProb(modj_warm)
        solj <- solBase
        names(solj$fluxes) <- model$modj@react_id 

    ##Copy optimization object to use as replacement for cases
    ##with low objective function value that require much more
    ##iterations in the next step
    copy <- cloneProbCPLEX(modj_warm@problem@oobj@env,modj_warm@problem@oobj@lp)

    out.gr <- solBase$fluxes[grep("EX_bio", model$modj@react_id)]
   
 ##Get reduced costs for exchanges -> don't need to add them since
    ##they have no influence on the objective function value
    redCosts <- getDjCPLEX(modj_warm@problem@oobj@env,modj_warm@problem@oobj@lp,0,model$modj@react_num-1)
        
library("MicrobiomeAGORA")
    interchange <- MicrobiomeAGORA::get_metabolic_interchange(model$modj, solj)
    na_ind <- is.na(names(solj[["fluxes"]]))
    res[[paste("orig",sep="-")]] <- list(solj = solj, community.growth=out.gr, met.interchange=interchange, models=model$model.IDs, fluxes = solj[["fluxes"]][!ind],redCosts=redCosts)

    ##Get base for quicker runs if objective has changed
    base <- getBaseCPLEX(modj_warm@problem@oobj@env,modj_warm@problem@oobj@lp)

    ##Do metabolite interventions
    start_time <- Sys.time()
for (l in 1:length(exIndices)) if(model$modj@react_id[exIndices[l]] %in% rownames(masses)){
      metaName <- model$modj@react_id[exIndices[l]]

      ##First removal of metabolite
      newLBs <- exLBs
      newLBs[l] <- exLBs[l]/100

      ## Reset old bounds
      if (l>1){
        chgColsBndsCPLEX(modj_warm@problem@oobj@env,
                         modj_warm@problem@oobj@lp,
                         ##exIndices[i-1]-1, ##Since CPLEX starts numbering with 0
                         ##exLBs[i-1],exUBs[i-1])
                         exIndices-1,
                         exLBs,exUBs)
      }

      chgColsBndsCPLEX(modj_warm@problem@oobj@env,
                       modj_warm@problem@oobj@lp,
                       exIndices[l]-1, ##Since CPLEX starts numbering with 0
                       newLBs[l],exUBs[l])

      ##Set previous basis
      copyBaseCPLEX(modj_warm@problem@oobj@env,
                    modj_warm@problem@oobj@lp,
                    base$cstat,base$rstat)

      ##Solve again
      solRem <- NULL
      if ((exLBs[l]!=0) && (solBase$fluxes[exIndices[l]]!=0)){
        ##if (exLBs[l]!=0){
        solRem <- optimizeProb(modj_warm)
        ##    print("Had to solve")
      }else{
        solRem <- solBase
      }
      if (solRem$stat==1){
        names(solRem$fluxes) <- model$modj@react_id
        out.gr <- solRem$fluxes[grep("EX_bio", model$modj@react_id)]
        interchange <- MicrobiomeAGORA::get_metabolic_interchange(model$modj, solRem)
        na_ind <- is.na(names(solRem[["fluxes"]]))
        res[[paste(metaName,"rem",sep="-")]] <- list(solj = solRem,community.growth=out.gr, met.interchange=interchange, models=model$model.IDs, fluxes = solRem[["fluxes"]][!na_ind])
}
      ##Second adding 2g of the metabolite

      ##Replace optimization object in case the previous optimization resulted in almost no growth
      if (solRem$obj/solBase$obj < 0.1){
        ##Free memory of previous problem
        delProbCPLEX(modj_warm@problem@oobj@env, modj_warm@problem@oobj@lp)

        modj_warm@problem@oobj@lp <- copy

        ##Create new copy
        copy <- cloneProbCPLEX(modj_warm@problem@oobj@env,modj_warm@problem@oobj@lp)
      }

      newLBs <- exLBs
      newLBs[l] <- newLBs[l] - masses[metaName, "intervention"]
      chgColsBndsCPLEX(modj_warm@problem@oobj@env,
                       modj_warm@problem@oobj@lp,
                       exIndices[l]-1, ##Since CPLEX starts numbering with 0
                       newLBs[l],exUBs[l])

      ##Set previous basis
      copyBaseCPLEX(modj_warm@problem@oobj@env,
                    modj_warm@problem@oobj@lp,
                    base$cstat,base$rstat)

      solAdd <- NULL
      if (redCosts[exIndices[l]]!=0)
        solAdd <- optimizeProb(modj_warm)
      else
        solAdd <- solBase

      if (solAdd$stat==1){
        names(solAdd$fluxes) <- model$modj@react_id
        out.gr <- solAdd$fluxes[grep("EX_bio", model$modj@react_id)]
        interchange <- MicrobiomeAGORA::get_metabolic_interchange(model$modj, solAdd)
        na_ind <- is.na(names(solAdd[["fluxes"]]))
        res[[paste(metaName,"add",sep="-")]] <- list(solj = solAdd,community.growth=out.gr, met.interchange=interchange, models=model$model.IDs, fluxes = solAdd[["fluxes"]][!na_ind])
      }
    }
   ##Free memory of CPLEX problem
    delProbCPLEX(modj_warm@problem@oobj@env, modj_warm@problem@oobj@lp)
  return(res)
 }
  names(out) <- colnames(communities)
  return(out)
}