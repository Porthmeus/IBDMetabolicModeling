# Porthmeus
# 07.10.21

sys.stdout = sys.stderr = open(snakemake.log[0],"w")

import cobra as cb
import pandas as pd
import numpy as np
import warnings
from corpse import simpleFastcore

dbg = False

if dbg:
    model = "resources/models/colormore3D.xml"
    diet = "resources/diets/colormore3D_MatjesAbsorption.csv"
    #core = "results/data/SPLIT_CoreRxns/SPLITGL25|L50|GU75_colormore22-F02243_L1_S10_L002.csv"
    core = "results/data/coreRxns/coreRxns.GL10-L50-GU90.colormore3D.csv"
    out = "temp/test.csv"
    sample = "H17446_L1_S13_L003"
    cnst = "results/data/consistentModels/CnstMod.MatjesAbsorption.colormore3D.csv"
else:
    model = snakemake.input["model"]
    diet = snakemake.input["diet"]
    core = snakemake.input["core"]
    cnst = snakemake.input["cnst_rxns"]
    invalRxn = snakemake.input["inval_rxns"]
    sample = snakemake.params["sample"]
    out = snakemake.output["cntxt"]

mod = cb.io.read_sbml_model(model)
diet = pd.read_csv(diet, index_col= 0)
core_rxn = pd.read_csv(core, index_col=0)
cnst_rxn = pd.read_csv(cnst, index_col=0)
if dbg:
    invalRxn = pd.DataFrame({"rxn":["RE1587R"]})
else:
    invalRxn = pd.read_csv(invalRxn, index_col= 0)

# adjust the diet
diet.loc[np.isnan(diet.iloc[:,0]),:] = 0
for rxn in mod.exchanges:
    if rxn.id in diet.index:
        rxn.lower_bound = -1*diet.loc[rxn.id,:][0]
    else:
        rxn.lower_bound = 0

# get the core set
rxn = core_rxn[core_rxn.loc[:,sample]==1].index.to_list()
# add the objective function(s) if not present
c_rxn = [x.id for x in mod.reactions if x.objective_coefficient == 1]
for c_r in c_rxn:
    if c_r not in rxn:
        rxn.append(c_r)

# remove one reactions which screws the whole thing

r = invalRxn["rxn"].to_list()
rxn2 = [x for x in rxn if x not in r]


# get the consistent model
cnst_rxns = cnst_rxn.loc[:,"rxn"].to_list()
noncon_rxns = [x.id for x in mod.reactions if x.id not in cnst_rxns]
mod2 = mod.copy()
mod2.remove_reactions(noncon_rxns)


# run fastcore 
fastcore = simpleFastcore(model = mod2, core_set = rxn2)
try:
    fastcore.fastcore()
    tis_mod = fastcore.get_model()
except Exception as e:
    print(e)
    warnings.warn("WARNING: fastcore is complaining about inconsistent model - will try to make it consistent with fastcc!")
    fastcore = simpleFastcore(model = mod, core_set= rxn2)
    fastcore.fastcc()
    fastcore.fastcore()
    tis_mod = fastcore.get_model()

# export the tissue specific model
tis_df = pd.DataFrame({"rxn" : [x.id for x in tis_mod.reactions]})
tis_df.to_csv(out)



#try:
#    # try to run the shortcut - for some reason that does not work very often, if so take the long run
#    fastcore = simpleFastcore(model = mod2, core_set = rxn)
#    fastcore.fastcore()
#    tis_mod = fastcore.get_model()
#except Exception as exc:
#    print(exc)
#    try:
#        zero_cutoff = mod2.tolerance/10
#        fastcore.zero_cutoff = zero_cutoff
#        fastcore.fastcore()
#        tis_mod = fastcore.get_model()
#    except Exception as exc:
#        print(exc)
#        try:
#            zero_cutoff = zero_cutoff/10
#            warnings.warn("WARNING:Initial consistent model is not working, will reduce the zero_cutoff to {cutoff}".format(cutoff = zero_cutoff))
#            fastcore.zero_cutoff = zero_cutoff
#            fastcore.fastcore()
#            tis_mod = fastcore.get_model()
#        except Exception as exc:
#            print(exc)
#            try:
#                zero_cutoff = zero_cutoff/10
#                warnings.warn("WARNING:Initial consistent model is not working, will reduce the zero_cutoff to {cutoff}".format(cutoff = zero_cutoff))
#                fastcore.zero_cutoff = zero_cutoff
#                fastcore.fastcore()
#                tis_mod = fastcore.get_model()
#            except Exception as exc:
#                print(exc)
#                zero_cutoff = 0.0
#                warnings.warn("WARNING:Initial consistent model is not working, will reduce the zero_cutoff to {cutoff}".format(cutoff = zero_cutoff))
#                fastcore.zero_cutoff = zero_cutoff
#                fastcore.fastcore()
#                tis_mod = fastcore.get_model()
#

