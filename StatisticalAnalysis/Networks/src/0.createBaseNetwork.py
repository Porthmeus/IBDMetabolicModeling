# Porthmeus
# 03.04.23

import cobra as cb
import pandas as pd
import numpy as np
from itertools import product

col3d = cb.io.read_sbml_model("../resources/Host/colormore3D_git.xml")

#rxn = col3d.reactions[2233]
#[x.id for x in rx.metabolites.keys()]

dat_all = pd.DataFrame()
for rxn in col3d.reactions:
    left = [x.id for x in rxn.metabolites if rxn.metabolites[x] <0]
    right = [x.id for x in rxn.metabolites if rxn.metabolites[x] >0]
    dat = pd.DataFrame(product(left,right), columns = ["from.id","to.id"])
    dat["reaction"] = rxn.id
    dat["subsystem"] = rxn.subsystem
    if rxn.reversibility:
        dat["arrows"] = "to;from"
    else: 
        dat["arrows"] = "to"
    dat_all = pd.concat([dat_all, dat])

dat_all.to_csv("../results/baseEdges.csv",index = False)

