# Porthmeus
# 21.01.22

# take the best explaining data set and prepare result tables for the presence/absence of reactions to HB/Mayo scores


# libraries
## @knitr loadLibraries
require(data.table)
require(caret)


# load data
## @knitr loadData_stat
FVA.center <- read.csv("../../../results/FVA/centerFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)

clinic <- fread("../../../resources/ClinicalData.csv")
meta <- fread("../../../resources/META_emed_future.csv")
blcklst <- fread("../../../resources/sampleBlacklist.csv")

# merge meta and clincal data
meta[,PatientID_Time := paste(PatientID, Time_seq, sep ="_")]
meta <- meta[Time_seq != -1,]
meta <- merge(meta,clinic, by="PatientID_Time", all.x=TRUE, suffix = c(".meta",""))


# remove degraded samples
blckSmpls <- blcklst[relevance < 2, Sample]
meta <- meta[!(SeqID %in% blckSmpls),]
FVA.center <- FVA.center[,!(colnames(FVA.center) %in% blckSmpls)]
meta <- meta[SeqID %in% colnames(FVA.center),]

# remove near zero variance reactions
near0var_idx <- apply(FVA.center,1, var) == 0 # nearZeroVar(t(FVA.center))
FVA.center <- FVA.center[!near0var_idx,]

# scale values
FVA.center <- t(scale(t(FVA.center)))

# get the exchange reactions
ex_rxns <- grep("^EX_", rownames(FVA.center), value =TRUE)

#sigs.FVA.HBMayo <- fread("../results/FVA_LMMsRandomPatient_fullCoefTable.csv")
#sigs.FVA.HBMayo[,padj := p.adjust(`Pr(>|t|)` , method = "BH")]
#sigs.PA.HBMayo <- fread("../results/PA_LMMsRandomPatient_fullCoefTable.csv")
#sigs.rxnExpr.HBMayo <- fread("../results/rxnExpr_LMMsRandomPatient_fullCoefTable.csv")
#sigs.FVA.Response <- fread("../results/FVA_GLMMsResponseAll_coefs.csv")
#sigs.PA.Response <- fread("../results/PA_GLMMsResponseAll_coefs.csv")
#sigs.rxnExpr.Response<- fread("../results/rxnExpr_GLMMsResponseAll_coefs.csv")
#sigs.FVA.Remission <- fread("../results/FVA_GLMMsRemissionAll_coefs.csv")
#sigs.PA.Remission <- fread("../results/PA_GLMMsRemissionAll_coefs.csv")
#sigs.rxnExpr.Remission<- fread("../results/rxnExpr_GLMMsRemissionAll_coefs.csv")
#sigs.FVA.Remission2 <- fread("../results/FVA_GLMMsRemissionWoHBMayoCorrection_coefs.csv")
#sigs.PA.Remission2 <- fread("../results/PA_GLMMsRemissionWoHBMayoCorrection_coefs.csv")
#sigs.rxnExpr.Remission2 <- fread("../results/rxnExpr_GLMMsRemissionWoHBMayoCorrection_coefs.csv")
#
#
#rxns.exsig <- unique(c(ex_rxns,
#                       sigs.rxnExpr.HBMayo[padj <0.05, rxn],
#                       sigs.PA.HBMayo[padj < 0.05, rn],
#                       sigs.FVA.HBMayo[padj < 0.05,rxn],
#                       sigs.rxnExpr.Response[padj < 0.05, rxn],
#                       sigs.PA.Response[padj < 0.05, rxn],
#                       sigs.FVA.Response[padj < 0.05,rxn],
#                       sigs.PA.Remission [padj < 0.05, rxn],
#                       sigs.rxnExpr.Remission[padj < 0.05, rxn],
#                       sigs.FVA.Remission [padj < 0.05, rxn],
#                       sigs.PA.Remission2 [padj < 0.05, rxn],
#                       sigs.rxnExpr.Remission2 [padj < 0.05, rxn],
#                       sigs.FVA.Remission2 [padj < 0.05, rxn]))

# get only the significant reactions
dr <- file.path("..","results")
sig.files <- list.files(dr, pattern = "_coefs.csv")
lst.rxns.exsig <- sapply(sig.files, function(x) fread(file.path(dr,x))[padj <= 0.05, rxn])

# combine exchange reactions and the significant ones to reduce dimensionality
rxns.exsig <- unique(c(ex_rxns,do.call(c,lst.rxns.exsig)))
rxns.exsig.bool <- rownames(FVA.center) %in% rxns.exsig 

# calculate correlation coefficients
corre <- cor(t(FVA.center[rxns.exsig.bool,]))


corre <- corre[ex_rxns, setdiff(colnames(corre),ex_rxns)]
corre.melt <- melt(data.table(corre, keep.rownames = TRUE),
                   id.vars = "rn",
                   variable.name = "rxn",
                   value.name = "r")
# add a absolute r
corre.melt[,r_abs:= sqrt(r^2)]

# write the table to a file
gz <- gzfile(file.path(dr, "FVA.center_ExCorrelations.csv.gz"), "w")
write.csv(corre.melt, row.names= FALSE, file = gz)
close(gz)
