# Porthmeus
# 25.01.21

# count the reactions and the reactions in the subsystems for each extracted model
sys.stdout = sys.stderr = open(snakemake.log[0], 'w')

import pandas as pd
import cobra as cb
import os
import numpy as np
import re

extModels = snakemake.input["cntxtMod"] # list of reactions of the context specific model
model = snakemake.input["model"] # original metabolic model to start with
samples = snakemake.params["samples"] # samples which have been used to extract the context specific models



# read the model to get the rxnIDs and the subsystems
mod_ori= cb.io.read_sbml_model(model)
rxns = [x.id for x in mod_ori.reactions]
groups = [x.id for x in mod_ori.groups]


# create data frame with zeros for each reaction and sample to fill it with ids afterwards
rxnIDMat = np.zeros((len(rxns),len(samples)), dtype = bool)
rxnIDMat = pd.DataFrame(rxnIDMat, columns = samples, index = rxns)


# read the extracted models and check the reactions therein
for sample in samples:
    
    extModel = [x for x in extModels if re.search(sample, x)][0]
    extMod = pd.read_csv(extModel)
    rxns = extMod.rxn
    rxnIDMat.loc[rxns, sample] = True

# count the reactions in the subsystems
subSysMat = np.zeros((len(groups), len(samples)),dtype = int)
subSysMat = pd.DataFrame(subSysMat, columns = samples, index = groups)

for group in groups:
    rxns = [x.id for x in mod_ori.groups.get_by_id(group).members]
    subSysMat.loc[group,:] = np.sum(rxnIDMat.loc[rxns,:],0)


# write the matrices to disk
rxnIDMat.to_csv(snakemake.output["rxnIdMat"])
subSysMat.to_csv(snakemake.output["subSysMat"])

