# miTarget_MM_pipeline

The pipeline is using the transcriptomic data to map these to the Recon3D model and extract context specific models from it.

# Usage

The pipeline expects a `resource` directory where all required files as input are stored. These files are the following:

| File                                             | Description                                                                                                                                                                                                                                                                                                                                                                                                             |
|--------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| TPM_matrix.csv                                   | A matrix file containing with genes/transcripts in rows and samples in columns. This should be some kind of normalized readcounts and should be comparable between genes **and** between samples. Usually expressed as TPM values.                                                                                                                                                                                          |
| META_data.csv                                    | Containing all the meta data relevant for the project. This table should contain as a minimum a column `SeqID` containing the sample identifier (column names) that have been used in the `TPM_matrix.csv`.                                                                                                                                                                                                                 |
| models/{model}.xml                               | The pipeline can handle multiple models at once. Every model where the data should be applied to, should be placed in the `resources/models/` subdirectory with the model name as filename. The model name will be used in the pipeline, thus be careful how you call it.                                                                                                                                                 |
| modelGene2TPM/{model}Genes2TPM.csv               | Different models use different namespaces, also in the GPRs. To match the gene/transcript identifiers. The table is expected to have a header (column names do not matter) with the first column corresponding to the gene.ids of the model GPRs and the second column being the gene.ids from the TPM_matrix.csv                                                                                                       |
| diet/{model}\_{dietName}.csv                      | Due to the same namespace issues, there needs to be a dedicated diet for each model. This should be placed into this subdirectory and with the corresponding model name and a diet name in the file name.                                                                                                                                                                                                               |
| invalidFastcoreRxns/invalRxns.{diet}.{model}.csv | Is a file you can add to exclude certain reactions from the core set. This is sometimes necessary, as in a very unfortunate constellation of core reactions. In that case fastcore fails, but excluding some reactions can help. If this is not needed add a file which just contains: `,"rxn"` and save it with the correct name. This features is deprecated and will be removed in the next iteration of the pipeline. |

After all files have been placed correctly, simply run the snakemake as ususal.

# Core reactions

Currently the pipeline runs a thresholding apporoach to define the core reaction which will be fed to fastcore. By default the pipeline applies 9 different thresholds - a mixture of global and local thesholds with different percentiles of the data.
\+ GL<percentile> - Global lower - all genes lower than that percentile are considered not active (not in core), percentile in comparison to all genes of all samples.
\- GU<percentile> - Global upper - all genes higher than that percentile are considered active (in core), percentile in comparison to all genes and all samples
\- L<percentile> - Local threshold - all genes, where this gene is higher than the percentile are considered active, percentile in comparison to only the gene in question

To change that behavior or the thresholds change the values in the snakemake file.

# CPLEX support

The Fastcore implementation was designed to be used with CPLEX, which requires manual installation as it is proprietary software. However, there is an academic license which is free to use if you are working for an academic institution. The registration and the download of CPLEX can be found [here](https://www.ibm.com/academic/topic/data-science). To configure python to work with CPLEX you need to tell python where to find your installation. In the environment requirements file `workflow/envs/python.yaml` adjust the `PYTHONPATH` variable to the directory of your installation of CPLEX.

If this is not properly configured/installed, fastcore will fall back onto GLPK as the linear solver. This can also work, but it might just fail as well.