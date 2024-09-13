# Porthmeus
#
# 27.10.20

# load the read count data, get the gene length information from the ensemble identifier and transform the data into TPM values

rule createTPMMatrix:
    input:
        RC = "resources/RC_emed_future.csv",
    output:
        TPM =  "results/data/TPM_emed_future.csv",
        FPKM = "results/data/FPKM_emed_future.csv"
    params:
        organism = "homo sapiens",
    conda:
        "../envs/R.yaml"
    log:"logs/createTPMMatrix.log"
    resources:
        mem_mb = 16000,
        time = "01:00:00",
        cpus_per_task = 1,
        task_per_node = 1
    script:"../scripts/createTPMMatrix.R"
