# Porthmeus
# 01.03.21

# combine the results of the FVA into a single matrix of ranges

rule summarizeFVA:
    input:
        rxns ="results/ModelAnalysis/rxnIdMat.{thrld}.{diet}.{model}.csv",
        samples = expand("results/data/FVA/FVAzeroBiomass.{{thrld}}.{sample}.{{diet}}.{{model}}.csv",
                sample = samples)
    params:
        sampleIDs = expand("{sample}",sample=samples)

    output:
        range = "results/FVA/rangeFVA.{thrld}.{diet}.{model}.csv",
        min = "results/FVA/minFVA.{thrld}.{diet}.{model}.csv",
        max = "results/FVA/maxFVA.{thrld}.{diet}.{model}.csv",
        center = "results/FVA/centerFVA.{thrld}.{diet}.{model}.csv",
        RDS = "results/FVA/FVA.{thrld}.{diet}.{model}.RDS"
    log : "logs/summaryFVA_{thrld}.{diet}.{model}.log"
    conda : "../envs/R.yaml"
    script : "../scripts/summarizeFVA2.R"
