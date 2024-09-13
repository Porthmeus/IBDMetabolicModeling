# Porthmeus
# 09.02.21

# calculate a FVA for each model from the extraction
sys.stdout = sys.stderr = open(snakemake.log[0], 'w')


import cobra
import pandas as pd
import numpy as np
import warnings

# load the data
model = cobra.io.read_sbml_model(snakemake.input["model"])
extract = pd.read_csv(snakemake.input["extract"])
diet = pd.read_csv(snakemake.input["diet"])
output = snakemake.output["out"]

# set the options
threads = snakemake.threads
frac_opt = snakemake.params["frac_opt"]
pfba_fac = snakemake.params["pfba_fac"]
setMaxBounds = snakemake.params["setMaxBounds"]

# reduce model
modRed = model.copy()
rxns = [x.id for x in modRed.reactions]

rm_rxns = [x for x in rxns if x not in extract.loc[:,"rxn"].to_list()]
# security check
if (extract.shape[0] + len(rm_rxns) != len(rxns)) or len(rm_rxns) == 0 or extract.shape[0] == 0:
    raise Exception("Something went wrong with the model extraction, please check the extracted reactions file and the model, compare namespaces!")

rm_rxns = [modRed.reactions.get_by_id(x) for x in rm_rxns]
modRed.remove_reactions(rm_rxns)

diet = diet.set_index("ex_rxn",drop = False)
# adjust diet
for ex_rxn in modRed.exchanges:
    ex_rxn.lower_bound = 0
    if ex_rxn.id in list(diet.index):
        if not pd.isna(diet.loc[ex_rxn.id][1]):
            ex_rxn.lower_bound = -1*diet.loc[ex_rxn.id][1]

# if toggled, set the boundaries to max 1000/-1000
inf = float("Inf")
if setMaxBounds:
    for rxn in modRed.reactions:
        if rxn.upper_bound == inf:
            rxn.upper_bound = 1000
        if rxn.lower_bound == -inf:
            rxn.lower_bound = -1000

# run the FVA
try:
    FVA = cobra.flux_analysis.flux_variability_analysis(modRed,
            processes = threads,
            fraction_of_optimum= frac_opt,
            pfba_factor= pfba_fac)
except Exception as optErr:
    err = optErr
    FVA = {"minimum" : [np.nan]*len(modRed.reactions), "maximum": [np.nan] * len(modRed.reactions)}
    FVA = pd.DataFrame(FVA, index = [x.id for x in modRed.reactions])
    warnings.warn(str(optErr) + "\n Will return a table with only NaN")

# write the output
FVA.to_csv(output)


