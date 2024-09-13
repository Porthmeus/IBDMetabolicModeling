# Porthmeus
# 05.10.21

import cobra as cb
import pandas as pd
import numpy as np
import corpse
import re


dbg = False

if dbg:
    TPM = "results/data/TPM_emed_future.csv"
    conversionTable = "resources/modelGenes2ensembl/colormore3DGenes2ensembl.csv"
    threshold_string = "GL25-L50-GU75"
    out = "results/data/coreRxns/coreRxns." + threshold_string + ".colormore3D.csv"
    rxnExpr = "results/RxnExpression/rxnExpr."+ threshold_string + ".colormore3D.csv"
    sbml = "resources/models/colormore3D.xml"
    threads = 3
else:
    # save stuff to the log file
    sys.stdout = sys.stderr = open(snakemake.log[0], 'w')
    TPM = snakemake.input["TPM"]
    conversionTable = snakemake.input["conversionTable"]
    threshold_string = snakemake.params["thresholds"]
    out = snakemake.output["CoreRxnMatrix"]
    rxnExpr = snakemake.output["TPMRxnMatrix"]
    sbml = snakemake.input["sbml"]
    threads = snakemake.threads


# extract the thresholds
thresholds = [float(re.sub("\D","",x)) if re.sub("\D","",x) != '' else None for x in threshold_string.split("-")]
while len(thresholds) < 3:
    thresholds.append(None)
gl,local,gu = thresholds


# load data
array = pd.read_csv(TPM, index_col = 0)
convTab = pd.read_csv(conversionTable, dtype = str)

# get the conversion table right - there might be genes missing in the TPM matrix, or there might be ambiguity between the two columns in the conversion table
subset = [i for i in range(convTab.shape[0]) if convTab.iloc[i,1] in array.index]
convTab = convTab.iloc[subset]
subset = list(convTab.iloc[:,1])
array = array.loc[subset,:]
array.index = convTab.iloc[:,0]
array = array.drop_duplicates()

# load the model
mod = cb.io.read_sbml_model(sbml)

# map the expression to the metabolic network
rxnArray = corpse.omicsMapper().mapExpressionToReaction(model = mod,
        dataframe = array, num_cores=threads)

#rxnArray = omicsMapper().mapExpressionToReaction(model = mod,
#        dataframe = array)

# get those reactions which fulfill the thresholds
coreRxn,thrlds_string = corpse.coreSetFinder().getCoreSet(array = rxnArray,
        global_lower = gl,
        global_upper = gu,
        local = local)

# write the binary matrix for the reactions
coreRxn.to_csv(out)

# write the mapping of the expression values to the model
rxnArray.to_csv(rxnExpr)

