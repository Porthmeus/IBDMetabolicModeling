# Porthmeus
# 07.10.21

sys.stdout = sys.stderr = open(snakemake.log[0],"w")

import cobra as cb
import pandas as pd
from corpse import simpleFastcore
import numpy as np

dbg = False

if dbg:
    model = "resources/models/colormore22.xml"
    diet = "resources/diets/colormore22_Matjes.csv"
    out = "results/data/consistentModels/CnstMod.Matjes.colormore22.csv"
else:
    model = snakemake.input["model"]
    diet = snakemake.input["diet"]
    out = snakemake.output["cnst_rxns"]

mod = cb.io.read_sbml_model(model)
diet = pd.read_csv(diet, index_col= 0)

# adjust the diet
diet.loc[np.isnan(diet.iloc[:,0]),:] = 0
for rxn in mod.exchanges:
    if rxn.id in diet.index:
        rxn.lower_bound = -1*diet.loc[rxn.id,:][0]
    else:
        rxn.lower_bound = 0

# run fastcc
fastcore = simpleFastcore(model = mod)
fastcore.FVA_consistency()
cnst_mod = fastcore.get_model()

# export the tissue specific model
cnst_df = pd.DataFrame({"rxn" : [x.id for x in cnst_mod.reactions]})
cnst_df.to_csv(out)

