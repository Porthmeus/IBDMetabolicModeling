# Porthmeus
# 18.08.21

# get the core genes dependent on the expression in each sample

rule getCoreRxns:
    input:
        TPM = "results/data/TPM_emed_future.csv",
        conversionTable = "resources/modelGenes2ensembl/{model}Genes2ensembl.csv",
        sbml = "resources/models/{model}.xml"
    output:
        CoreRxnMatrix = "results/data/coreRxns/coreRxns.{thrld}.{model}.csv",
        RxnExprMatrix = "results/RxnExpression/rxnExpr.{thrld}.{model}.csv"
    log: "logs/getCoreRxns.{thrld}_{model}.csv"
    conda: "../envs/python.yaml"
    threads: 8
    params:
        thresholds = "{thrld}"
    script:
        "../scripts/getCoreRxns.py"
