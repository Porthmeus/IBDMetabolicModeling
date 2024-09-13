# Porthmeus
# 11.10.21

rule getContextModel:
    input:
        model = "resources/models/{model}.xml",
        diet = "resources/diets/{model}_{diet}.csv",
        core = "results/data/coreRxns/coreRxns.{thrld}.{model}.csv",
        cnst_rxns = "results/data/consistentModels/CnstMod.{diet}.{model}.csv",
        inval_rxns = "resources/invalidFastcoreRxns/invalRxns.{diet}.{model}.csv"
        #inval_rxns = "results/data/invalidFastcoreRxns/invalRxns.{diet}.{model}.csv"
    params:
        sample = "{sample}"
    output:
        cntxt = "results/data/ContextSpecificModels/CntxtRxns.{thrld}.{sample}.{diet}.{model}.csv"
    log: "logs/getContextModel_{thrld}.{sample}.{diet}.{model}.log"
    conda: "../envs/python.yaml"
    threads: 2
    resources:
        mem_mb = 8000
    script: "../scripts/getContextModel.py"
