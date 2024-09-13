# Analysis of the metabolic network reconstructions for IBD patients

The main idea of the project is to reconstruct metabolic networks for tissue samples of the gut of IBD patients. These will be used to run an FVA in order to aim for three things:

1. Understand the etiology of the disease on the metabolic level
2. Find biomarkers which can predict the Response/Remission of a patient for a treatment
3. Find targets for a possible intervention therapy to directly treat, or improve therapy

The analysis are based on the following data sources:

1. The expression values mapped to the metabolic network via the GPR rules - termed **rxnExpr**
2. The presence/absence (**PA**) of the reaction after reconstruction in the network
3. The results of the **FVA** converted to ranges and centers for each reaction

Following analysis are being done:

1. A (G)LMM for each reaction individually to predict the given outcome
2. A subsystem enrichment as one of the tow:
    1. GSEA - gene set enrichment - based on the estimates from the (G)LMMs
    2. HGT - hyper-geometric test (fisher's exact) - based on the reactions significantly associated by the (G)LMMs
3. An analysis of the involved metabolites

## Etiology

Is the simplest of the analysis - basically we correlate the different data sources to the change of the clinical scores (HB/Mayo). Currently we consider only the intra-patient dependence of the data as a random factor in the models.

As an additional layer of information, I added an analysis with an interaction term with the disease, in order to get specific reactions which change for each of the diseases. This might become important later on.


## Biomarkers

Is more complicated, because we need to exclude certain time points. For Remission, we exclude the last time point for the inference of the reactions important for the association with HB/Mayo. We are still considering which data to use for the Response of the patients. I tend to include only data up until day 14 to predict Response.

## Intervention

This analysis is with Response or Remission as outcome variable, but the reaction prediction is corrected by the inflammatory state of the patients - HB/Mayo is a confounding effect in the model analysis.

