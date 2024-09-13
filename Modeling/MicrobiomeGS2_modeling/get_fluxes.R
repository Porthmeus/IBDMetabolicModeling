# extract flux and flux.o from community models
# Rscript get_fluxes.R path_to_comm_models.RDS path_to_output

rm(list=ls())
fncols <- function(data, cname) { 
  add <-cname[!cname%in%names(data)] #what are the columns that do not exist in the column names of data
  if(length(add)!=0) data[add] <- NA #if there are some then, give them a value of NA in the dataframe
  data
}
# a function to check if columns exist in a dataframe
# and to add them with NA value if they don't exist
###########
args <- commandArgs(trailingOnly = TRUE)
comm_files <- paste0(args[[1]], list.files(args[[1]]))
all_flux <-as.data.frame(t(matrix(0,1,1)))
all_flux.o <- as.data.frame(t(matrix(0,1,1)))
all_react <- as.data.frame(t(matrix(0,1,1)))
growth <- c()
com_count <- 0
#i=2
for (i in 1:length(comm_files)){
  current_comms <- readRDS(comm_files[[i]])
  current_comms <- current_comms[!sapply(current_comms,is.null)] #remove null from community models
  print(paste0(i, "_", length(current_comms)))
#  l=1
  if (length(current_comms)==0){
    next
    }
    for (l in 1:length(current_comms)){
        com_count = com_count+1
        growth <- c(growth, current_comms[[l]][["solj"]][["obj"]])
        samp_metabolites <- current_comms[[l]][["met.interchange"]][["rxn"]]
        samp_reacts <- current_comms[[l]][["solj"]][["fluxes"]]
        reacts_df <- as.data.frame(samp_reacts)
        names(reacts_df) <- "flux"
        reacts_df$react <- names(samp_reacts)
        reacts_df <- reacts_df[!is.na(reacts_df$react),]
        reacts_df <- reacts_df[grepl("rxn", reacts_df$react),]
        reacts_df$react <- gsub("M\\d+_", "", reacts_df$react)
        reacts_df <- aggregate(reacts_df$flux,
                               list(reacts_df$react), FUN = sum)
        names(reacts_df) <- c("react", "flux")
        all_react <- fncols(all_react, reacts_df$react)
        ind <- match(reacts_df$react, colnames (all_react))
        all_react[com_count,ind] = reacts_df$flux
        rownames(all_react)[[com_count]] <- names(current_comms)[[l]]
        all_flux <- fncols(all_flux, samp_metabolites)
        ind2 <- match(samp_metabolites, colnames (all_flux))
        all_flux[com_count,ind2] = current_comms[[l]][["met.interchange"]][["flux"]]
        rownames(all_flux)[[com_count]] <- names(current_comms)[[l]]
        all_flux.o <- fncols(all_flux.o, samp_metabolites)
        ind3 <- match(samp_metabolites, colnames (all_flux.o))
        all_flux.o[com_count,ind3] = current_comms[[l]][["met.interchange"]][["o.flux"]]
        rownames(all_flux.o)[[com_count]] <- names(current_comms)[[l]]
    }
#    all_flux[abs(all_flux)<10^-6] =0
#    all_flux.o[abs(all_flux.o)<10^-6] =0
#    all_react[abs(all_react)<10^-6] =0
    all_flux$growth <- growth
}
all_flux = subset(all_flux, select = -c(V1))
all_flux.o = subset(all_flux.o, select = -c(V1))
all_react = subset(all_react, select = -c(V1))
growth <- as.data.frame(growth)
rownames(growth)<-rownames(all_flux)
write.csv(all_flux, paste0(args[2], "flux_among_bacteria.csv"))
write.csv(all_flux.o, paste0(args[2], "flux_with_host.csv"))
write.csv(all_react, paste0(args[2], "internal_reacts.csv"))
write.csv(growth, paste0(args[2], "community_growth.csv"))