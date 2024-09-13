#Porthmeus
# 23.03.22


require(data.table)
require(stringr)

# simple script to get a table which holds the most relevant reaction-to-exchange relations
FVA.center.excor <- fread("../results/FVA.center_ExCorrelations.csv.gz")
FVA.range.excor <- fread("../results/FVA.range_ExCorrelations.csv.gz")
PA.excor <- fread("../results/PA_ExCorrelations.csv.gz")


# read the statistical results
res.dr <- file.path("..","results")
coef.files <- list.files(res.dr, pattern = "_coefs.csv")
lst.coefs <- lapply(coef.files, function(x) fread(file.path(res.dr,x))[padj <=0.05,] )
names(lst.coefs) <- coef.files
# split the information of the file name into seperate columns
for(f in coef.files){
    name.tab <- t(
                  as.data.frame(
                                unlist(
                                       strsplit(
                                                gsub("_coefs.csv","",f),
                                                split = "\\."
                                                )
                                       )
                                )
                  )
    colnames(name.tab) <- c("data_src","stat_model","est_src","aim")
    lst.coefs[[f]] <- cbind(lst.coefs[[f]], name.tab)
    # add some more information to the FVA
    if(name.tab[1,"data_src"] == "FVA"){
        lst.coefs[[f]][,coef:=paste(name.tab[1,"data_src"], coef, sep = ".")]
    } else {
        lst.coefs[[f]][,coef:=name.tab[1,"data_src"]]
    }

}
# combine to a single table
clmn.cmn <- unique(unlist(sapply(lst.coefs, colnames)))
clmn.cmn.diff <- unique(unlist(sapply(lapply(lst.coefs,colnames), setdiff, x = clmn.cmn)))
clmn.cmn <- setdiff(clmn.cmn,clmn.cmn.diff)
sigs.all <- do.call(rbind,lapply(lst.coefs, function(x) x[,..clmn.cmn]))

#sigs.FVA.HBMayo <- fread("../results/FVA_LMMsRandomPatient_fullCoefTable.csv")
#sigs.FVA.HBMayo[,padj := p.adjust(`Pr(>|t|)` , method = "BH")]
#sigs.PA.HBMayo <- fread("../results/PA_LMMsRandomPatient_fullCoefTable.csv")
#sigs.PA.HBMayo[,rxn := rn]
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
## correct the correlation coefficient for the FVA centers
#FVA.center.excor[,r := -1*r]
#
## correct the LMM coefficients for the Exchange reactions for explaining the HB/Mayo coefficients:
#sigs.FVA.HBMayo[coef == "center" & grepl("^EX_", rxn), Estimate := Estimate * -1]
#sigs.FVA.Response[coef == "center" & grepl("^EX_", rxn), Estimate := Estimate * -1]
#sigs.FVA.Remission[coef == "center" & grepl("^EX_", rxn), Estimate := Estimate * -1]
#sigs.FVA.Remission2[coef == "center" & grepl("^EX_", rxn), Estimate := Estimate * -1]

# combine the correlation tables
FVA.center.excor[, cor_src := "FVA.center"]
FVA.range.excor[, cor_src := "FVA.range"]
PA.excor[, cor_src := "PA"]
All.excor <- rbind(FVA.center.excor,
                   FVA.range.excor,
                   PA.excor)

# calculate summary statistics across the different data set correlations
All.excor.sum <- All.excor[, .(
                               bst_cor = cor_src[which.max(r_abs)],
                               r_mean = mean(r),
                               r_abs_mean = mean(r_abs),
                               N_sets = .N,
                               r_maxDev = diff(range(r)),
                               r_abs_maxDev = diff(range(r_abs)),
                               direction = -1*sum(r<0) +sum(r>0)
                               ), by = .(rn,rxn)]






### HB/Mayo
#sigs.FVA.HBMayo <- sigs.FVA.HBMayo[padj < 0.05, .(rxn, coef = paste("FVA",coef,sep="."), Estimate, est_src = "HBMayo")]
#sigs.PA.HBMayo <- sigs.PA.HBMayo[padj < 0.05, .(rxn, coef = "PA", Estimate,  est_src = "HBMayo")]
#sigs.rxnExpr.HBMayo <- sigs.rxnExpr.HBMayo[padj < 0.05, .(rxn, coef = "rxnExpr", Estimate,  est_src = "HBMayo")]
#
### Response
#sigs.FVA.Response <- sigs.FVA.Response[padj < 0.05, .(rxn, coef = paste("FVA",coef,sep="."), Estimate, est_src = "Response")]
#sigs.PA.Response <- sigs.PA.Response[padj < 0.05, .(rxn, coef = "PA", Estimate,  est_src = "Response")]
#sigs.rxnExpr.Response <- sigs.rxnExpr.Response[padj < 0.05, .(rxn, coef = "rxnExpr", Estimate,  est_src = "Response")]
#
### Remission
#sigs.FVA.Remission <- sigs.FVA.Remission[padj < 0.05, .(rxn, coef = paste("FVA",coef,sep="."), Estimate, est_src = "Remission")]
#sigs.PA.Remission <- sigs.PA.Remission[padj < 0.05, .(rxn, coef = "PA", Estimate,  est_src = "Remission")]
#sigs.rxnExpr.Remission <- sigs.rxnExpr.Remission[padj < 0.05, .(rxn, coef = "rxnExpr", Estimate,  est_src = "Remission")]
#
## Remission2
#sigs.FVA.Remission2 <- sigs.FVA.Remission2[padj < 0.05, .(rxn, coef = paste("FVA",coef,sep="."), Estimate, est_src = "Remission2")]
#sigs.PA.Remission2 <- sigs.PA.Remission2[padj < 0.05, .(rxn, coef = "PA", Estimate,  est_src = "Remission2")]
#sigs.rxnExpr.Remission2 <- sigs.rxnExpr.Remission2[padj < 0.05, .(rxn, coef = "rxnExpr", Estimate,  est_src = "Remission2")]
#
#sigs.all <- rbind(sigs.FVA.HBMayo,
#                  sigs.PA.HBMayo,
#                  sigs.rxnExpr.HBMayo,
#                  sigs.FVA.Response,
#                  sigs.PA.Response,
#                  sigs.rxnExpr.Response,
#                  sigs.FVA.Remission,
#                  sigs.PA.Remission,
#                  sigs.rxnExpr.Remission,
#                  sigs.FVA.Remission2,
#                  sigs.PA.Remission2,
#                  sigs.rxnExpr.Remission2)


sigs.all.sum <- sigs.all[,.(
                                bst_est = coef[which.max(abs(Estimate))],
                                mean_est = mean(Estimate),
                                mean_est_abs = mean(sqrt(Estimate^2)),
                                max_est = Estimate[which.max(abs(Estimate))],
                                N_sets = .N,
                                direction = -1*sum(Estimate < 0) + sum(Estimate>0)
                                ), by = .(rxn, est_src,aim,stat_model)]
# leave only those correlations which have at least |rho| > 0.5 to filter week correlations
All.excor.sum <- All.excor.sum[r_abs_mean > 0.5,]

# calculate a score for the correlation coefficient
# the score is the product of the absolute correlation coefficient and the mean estimate for the respective reactions to its factor (HB/Mayo or Response) in question: r_mean * mean_est

All.excor.sum.extended <- merge(All.excor.sum, sigs.all.sum, by = "rxn", all.x=TRUE, allow.cartesian=TRUE)
All.excor.sum.sum <- All.excor.sum.extended[,.(bst_cor = bst_cor,
                                               r_mean = r_mean,
                                               r_abs_mean = r_abs_mean,
                                               r_sets = N_sets.x,
                                               r_abs_maxDev = r_abs_maxDev,
                                               r_direction = direction.x,
                                               bst_est = bst_est,
                                               est_mean = mean_est,
                                               est_mean_abs = mean_est_abs,
                                               est_max = max_est,
                                               est_sets = N_sets.y,
                                               est_direction = direction.y,
                                               score_mean = mean(r_mean*mean_est),
                                               score_absabs_mean = mean(r_abs_mean * mean_est_abs),
                                               score_abs_mean = mean(r_abs_mean * mean_est)
                                               ),
                                            by = .(rxn, rn, est_src,aim,stat_model)]


#setkey(sigs.all.sum, "rxn")
#reac <- All.excor.sum[,rxn]
#All.excor.sum[, HBMayo_score := sigs.all.sum[est_src == "HBMayo",][reac,mean_est,allow.cartesian=TRUE] * r_mean]
#All.excor.sum[, Response_score := sigs.all.sum[est_src == "Response",][reac,mean_est,allow.cartesian=TRUE] * r_mean]
#All.excor.sum[, Remission_score := sigs.all.sum[est_src == "Remission",][reac,mean_est,allow.cartesian=TRUE] * r_mean]
#All.excor.sum[, Remission2_score := sigs.all.sum[est_src == "Remission2",][reac,mean_est,allow.cartesian=TRUE] * r_mean]
#

## leave only the top1 percent of exchange reactions per reaction in the table
#top1 <- All.excor.sum[,.(HB = abs(HBMayo_score) > quantile(HBMayo_score,0.99,na.rm =TRUE), 
#                         Resp = abs(Response_score) > quantile(Response_score,0.99,na.rm =TRUE)
#                         ), by = "rxn"]
#
#top1[,keep:=FALSE]
#top1[HB==TRUE | Resp == TRUE,keep :=TRUE]
#All.excor.sum <- All.excor.sum[top1[,keep],]


# add the scores of the exchange reactions which have been directly detected in the initial analysis
ex.rxns <- grep("^EX_", sigs.all.sum[,rxn], value=TRUE)
sigs.ex.sum <- sigs.all.sum[rxn %in% ex.rxns,]

#sigs.ex.sum <- dcast(sigs.all.sum[rxn %in% ex.rxns,], rxn~est_src, value.var = "mean_est")
#ex.excor.sum <- sigs.ex.sum[,.(
#                                                  rn = rxn,
#                                                  rxn = rxn,
#                                                  bst_cor = NA,
#                                                  r_mean = 1,
#                                                  r_abs_mean = 1,
#                                                  N_sets = NA,
#                                                  r_maxDev = 0,
#                                                  r_abs_maxDev = 0,
#                                                  direction = 1,
#                                                  HBMayo_score = HBMayo,
#                                                  Response_score = Response,
#                                                  Remission_score = Remission,
#                                                  Remission2_score = Remission2)]

ex.excor.sum <- sigs.ex.sum[,.(
                               rn = rxn,
                               rxn = rxn,
                               est_src = est_src,
                               aim = aim,
                               stat_model = stat_model,
                               bst_cor = NA,
                               r_mean = 1,
                               r_abs_mean = 1,
                               r_sets = 3,
                               r_abs_maxDev = 0,
                               r_direction = 3,
                               bst_est = bst_est, 
                               est_mean = mean_est,
                               est_mean_abs= mean_est_abs,
                               est_max = max_est,
                               est_sets = N_sets,
                               est_direction = direction,
                               score_mean = mean_est,
                               score_absabs_mean = mean_est_abs,
                               score_abs_mean = mean_est)
                            ]
clmns.excor <- setdiff(colnames(All.excor.sum.sum),colnames(ex.excor.sum))
extra.clmns<-  data.table(matrix(NA, ncol = length(clmns.excor), nrow = 1, dimnames = list(x=1,y=clmns.excor)))
ex.excor.sum <- cbind(ex.excor.sum,extra.clmns)
clmns.excor <- colnames(All.excor.sum.sum)
ex.excor.sum <-  ex.excor.sum[,..clmns.excor]
All.excor.sum <- rbind(All.excor.sum.sum, ex.excor.sum)

# finally add the absolute scores for filtering
#All.excor.sum[,HBMayo_score_abs := sqrt(HBMayo_score^2)]
#All.excor.sum[,Response_score_abs := sqrt(Response_score^2)]
#All.excor.sum[,Remission_score_abs := sqrt(Remission_score^2)]
#All.excor.sum[,Remission2_score_abs := sqrt(Remission2_score^2)]

# to get an idea how to interact with the metabolites in terms of therapy, add the information if the metabolite is usually taken up or exported
fva.max <- read.csv("../../../results/FVA/maxFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names= 1)
fva.min <- read.csv("../../../results/FVA/minFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names = 1)

# get only those reactions still in the table
ex.rxns <- unique(All.excor.sum[,rn])
fva.max <- fva.max[ex.rxns,]
fva.min <- fva.min[ex.rxns,]

# remove numerical instabilities
fva.max[abs(fva.max) <1E-6] <- 0
fva.min[abs(fva.min) <1E-6] <- 0

# check for each reaction and sample whether it is used for export and/or for import
#import.export <- (apply(fva.max, 1, function(x) x>0) - apply(fva.min, 1, function(x) x<0) +2)

# 0 = no import nor export, 1 = import, 2 = both, 3 = export
import.export <- ( (fva.max > 0) - (fva.min < 0) +2 )  +(-2*(fva.max == 0 & fva.min ==0)) 

# get the majority vote for each reaction
i.e <- apply(import.export, 1, table)
i.e <- sapply(i.e, function(x) as.integer(names(sort(x, decreasing = TRUE))))
i.e <- sapply(i.e, function(x) c("import","both","export")[x][1])

# add the information to the table
All.excor.sum[,main.direction := i.e[rn]]

# add some human readible information
subsystems <- fread("../dat/subsystems.csv")
rxnAnno <- fread("../../../../resources/recon-store-reactions-1.tsv")
rxnAnno[, Metabolite := gsub(" exchange", "", gsub("Exchange of ", "", description))]
rxnAnno[, rxnBase := gsub("\\[.\\]","",abbreviation)]

# add the subsystems
All.excor.sum <- merge(All.excor.sum, subsystems, all.x = TRUE, by.x = "rxn", by.y = "reaction")

# add the metabolites which are associated with the exchange
All.excor.sum[,rxnBase := gsub("\\[.\\]","",rn)]
All.excor.sum <- merge(All.excor.sum, unique(rxnAnno[,.(rxnBase, Metabolite)]), by = "rxnBase", all.x = TRUE)

# add the reaction descriptions
All.excor.sum <- merge(All.excor.sum, unique(rxnAnno[,.(abbreviation, description)]), by.x = "rxn", by.y = "abbreviation", all.x = TRUE)

# add the information which compartment the exchange is directed to
All.excor.sum[,compartment:= "Colon"]
All.excor.sum[grep("EX_.*\\[e\\]", rn),compartment:= "Blood"]
 

write.csv(file="../results/ExCorrelation_summary2.csv", All.excor.sum, row.names = FALSE)

