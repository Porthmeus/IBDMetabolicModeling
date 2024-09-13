# Porthmeus
# 11.02.21

# calculate a the FVA for each extracted model
rule calcFVA_zeroBiomass:
    input:
        model = "resources/models/{model}.xml",
        diet = "resources/diets/{model}_{diet}.csv",
        extract = "results/data/ContextSpecificModels/CntxtRxns.{thrld}.{sample}.{diet}.{model}.csv"
    output:
        out = temp("results/data/FVA/FVAzeroBiomass.{thrld}.{sample}.{diet}.{model}.csv")
    log:
        "logs/FVAzeroBiomass.{thrld}.{sample}.{diet}.{model}.log"
    conda:
        "../envs/python.yaml"
    params:
        frac_opt = 0,
        pfba_fac = None,
        setMaxBounds = True
    resources:
        mem_mb = "16000",
        time = "08:00:00",
        task_per_node = 1,
    threads: 1
    conda: "../envs/python.yaml"
    script: "../scripts/calcFVA.py"

        
