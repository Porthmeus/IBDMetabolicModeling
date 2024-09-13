construct_community_models_single <- function(mic.table, i, models, abunCutoff = 0.001, nutrition=NULL, modelName=NULL) {
    library("MicrobiomeGS2")
    source("/zfshome/sukem124/Rscripts/adjustNutrition.R")
  mod.rel <- t(t(mic.table[-1,])/colSums(mic.table[-1,])) #make the relative abundancies
  mod.rel <- data.table(as.table(mod.rel))
  colnames(mod.rel) <- c("spec","sample","ratio") # spec = species
  mod.rel <- mod.rel[ratio >= abunCutoff]

  spec.ratio <- list() #creat a list of the the species abundance in each sample (but it is not needed to do it on all samples becuase this runs only on one sample (this should be changed in the future by other intelegent life forms))
  for(l in 1:ncol(mic.table)) {
    spec.ratio[[l]] <- mod.rel[sample==colnames(mic.table)[i],.(spec,ratio)]
  }
  ind <- which(names(models) %in% spec.ratio[[i]]$spec) #get the indecies of the bacterial models that we want to include in the community
  model.joined <-  MicrobiomeGS2::join_mult_models(models[ind]) #a microbiomeGS function to creat the compartments and put the models in them

  # introduce new reaction (community biomass) with the biomasses of the bacteria as constraints (stochiometries)

  model.joined$model.IDs <- merge(model.joined$model.IDs,spec.ratio[[i]],by.x="model.name",by.y="spec")

  bm.mets <- paste0(model.joined$model.IDs$model.id,"_cpd11416[c0]") #the id of the biomass reaction


  model.joined$modj <- addReact(model = model.joined$modj, id = "EX_bio1", met = bm.mets, Scoef = -model.joined$model.IDs$ratio, reversible = F, reactName = "joined Biomass with fixed ratios") #add the new reaction that we created

  model.joined$modj@obj_coef <- rep(0, length(model.joined$modj@react_id))
  model.joined$modj@obj_coef[grep("EX_bio1", model.joined$modj@react_id)] <- 1 #make the new reaction as the objectiove function

  model.joined$modj@uppbnd[grep("M[0-9]+_EX_cpd11416", model.joined$modj@react_id)] <- 0
  # get coupling constraints for later optimizations
  model.joined$coupling <- MicrobiomeGS2::get_coupling_constraints_mult(model.joined$modj)

    if(!is.null(nutrition)){
      nut_ind <- which(rownames(nutrition) == modelName)
      model.joined <- adjustNutrition(model.joined, nutrition[nut_ind,])
    }
    return(model.joined)
}
