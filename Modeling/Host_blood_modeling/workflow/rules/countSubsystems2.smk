# Porthmeus
# 19.01.21

# load all the extraced models and count the reactions in the subsystem as well as their presence in general

rule countSubsystems:
    input:
        cntxtMod = expand("results/data/ContextSpecificModels/CntxtRxns.{{thrld}}.{sample}.{{diet}}.{{model}}.csv",
                sample = samples),
        model = "resources/models/{model}.xml"
    output:
        rxnIdMat = "results/ModelAnalysis/rxnIdMat.{thrld}.{diet}.{model}.csv",
        subSysMat = "results/ModelAnalysis/SubSysMat.{thrld}.{diet}.{model}.csv"
    log: "logs/countSubsystems.{thrld}.{diet}.{model}.log"
    params:
        samples = expand("{sample}", sample = samples)
    resources:
        mem_mb = 16000,
        time = "05:00:00",
        cpus_per_task = 8,
        task_per_node = 1
    conda: "../envs/python.yaml"
    script: "../scripts/countSubsystems2.py"

