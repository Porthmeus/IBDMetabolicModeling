# Porthmeus
# 10.03.22

# create some diversity measures for the different analysis values we have

require(data.table)
require(vegan)
require(lme4)
require(lmerTest)
require(ggplot2)
require(cowplot)

# load data
rxnCount <- read.csv("../../../results/ModelAnalysis/rxnIdMat.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
rxnExpr <- read.csv("../../../results/RxnExpression/rxnExpr.GL10-L50-GU90.colormore3D.csv", row.names = 1)
# load FVA results
FVA.min <- read.csv("../../../results/FVA/minFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
FVA.max <- read.csv("../../../results/FVA/maxFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
# convert the values to ranges and centers
FVA.range <- FVA.max-FVA.min
FVA.range[abs(FVA.range) < 1E-6] <- 0
FVA.center <- (FVA.max+FVA.min)/2
FVA.center[abs(FVA.center) < 1E-6] <- 0


clinic <- fread("../../../resources/ClinicalData.csv")
meta <- fread("../../../resources/META_emed_future.csv")
blcklst <- fread("../../../resources/sampleBlacklist.csv")

# make rxnCounts a integer matrix
rxnCount2 <- rxnCount
rxnCount2[,] <- 0
rxnCount2[rxnCount == "True"] <- 1
rxnCount <- rxnCount2

# merge meta and clincal data
meta[,PatientID_Time := paste(PatientID, Time_seq, sep ="_")]
meta <- meta[Time_seq != -1,]
meta <- merge(meta,clinic, by="PatientID_Time", all.x=TRUE, suffix = c(".meta",""))


# remove degraded samples
blckSmpls <- blcklst[relevance < 2, Sample]
meta <- meta[!(SeqID %in% blckSmpls),]
meta <- meta[SeqID %in% colnames(rxnCount),]

rxnExpr <- rxnExpr[,!(colnames(rxnExpr) %in% blckSmpls)]
rxnCount <- rxnCount[,!(colnames(rxnCount) %in% blckSmpls)]
sel<-!(colnames(FVA.min) %in% blckSmpls)
FVA.range <- FVA.range[,sel]
FVA.center <- FVA.center[,sel]


# simply count the reactions retrieved per sample
shannon.rxnExpr <- diversity(t(rxnExpr), index = "shannon")
richness.PA <- colSums(rxnCount)/nrow(rxnCount)
shannon.FVArange <-diversity(t(FVA.range), index = "shannon")
# transform to positive values for the center
shannon.FVAcenter <-diversity(t(abs(FVA.center)), index = "shannon")
 
meta[, richness.PA:= richness.PA[SeqID]]
meta[, shannon.rxnExpr:= shannon.rxnExpr[SeqID]]
meta[, shannon.FVArange:= shannon.FVArange[SeqID]]
meta[, shannon.FVAcenter:= shannon.FVAcenter[SeqID]]
meta[, modelSize := colSums(rxnCount)[SeqID]]

# plot diversity vs. HBMayo
pa.mod <- lmer(data = meta,
               formula = HB_Mayo_impu ~ richness.PA + (1|PatientID))
pa.plt <- ggplot(meta, aes(y = richness.PA, x = HB_Mayo_impu))+
    geom_point() +
    annotate(geom = "text", x = min(meta[,HB_Mayo_impu]),y=max(meta[,richness.PA]) , label = paste("p =", round(coef(summary(pa.mod))[2,5], digits = 3)), hjust = 0, vjust = 0)+
    ggtitle("PA richness")+
    geom_smooth(method = "lm")+
    theme_bw()

rxnExpr.mod <- lmer(data = meta,
               formula = HB_Mayo_impu ~ shannon.rxnExpr + (1|PatientID))
rxnExpr.plt <- ggplot(meta, aes(y = shannon.rxnExpr, x = HB_Mayo_impu))+
    geom_point() +
    annotate(geom = "text", x = min(meta[,HB_Mayo_impu]),y=max(meta[,shannon.rxnExpr]) , label = paste("p =", round(coef(summary(rxnExpr.mod))[2,5], digits = 3)), hjust = 0, vjust = 0)+
    ggtitle("RxnExpr Shannon")+
    geom_smooth(method = "lm")+
    theme_bw()

FVArange.mod <- lmer(data = meta,
               formula = HB_Mayo_impu ~ shannon.FVArange + (1|PatientID))
FVArange.plt <- ggplot(meta, aes(y = shannon.FVArange, x = HB_Mayo_impu))+
    geom_point() +
    annotate(geom = "text", x = min(meta[,HB_Mayo_impu]),y=max(meta[,shannon.FVArange]) , label = paste("p =", round(coef(summary(FVArange.mod))[2,5], digits = 3)), hjust = 0, vjust = 0)+
    stat_smooth(method = "lm") +
    ggtitle("FVA range Shannon")+
    theme_bw()

FVAcenter.mod <- lmer(data = meta,
               formula = HB_Mayo_impu ~ shannon.FVAcenter + (1|PatientID))
FVAcenter.plt <- ggplot(meta, aes(y = shannon.FVAcenter, x = HB_Mayo_impu))+
    geom_point() +
    annotate(geom = "text", x = min(meta[,HB_Mayo_impu]),y=max(meta[,shannon.FVAcenter]) , label = paste("p =", round(coef(summary(FVAcenter.mod))[2,5], digits = 3)), hjust = 0, vjust = 0)+
    stat_smooth(method = "lm") +
    ggtitle("FVA center Shannon")+
    theme_bw()

div.plt <- plot_grid(rxnExpr.plt, pa.plt, FVArange.plt, FVAcenter.plt)
ggsave(div.plt, file = "../results/diversityPlots.HBMayo.pdf", width = 8, height = 8)



# save the meta data with the diverstiy measures
write.csv(meta, file = "../results/diversityMeasures.csv", row.names = FALSE)

# Remission
rem.dat <- meta[Remission != "C",]
rem.dat[, Remission := factor(Remission)]
pa.mod <- glmer(data = rem.dat, 
                family = "binomial",
                formula = Remission ~ richness.PA + HB_Mayo_impu + (1|PatientID))
pa.plt <- ggplot(rem.dat, aes(y = richness.PA, x = Remission))+
    geom_boxplot() +
    annotate(geom = "text", x = 0.5,y=max(rem.dat[,richness.PA]) , label = paste("p =", round(coef(summary(pa.mod))[2,4], digits = 3)), hjust = 0, vjust = 0)+
    ggtitle("PA richness")+
    theme_bw()

rxnExpr.mod <- glmer(data = rem.dat,
                family = "binomial",
               formula = Remission ~ shannon.rxnExpr + HB_Mayo_impu+(1|PatientID))
rxnExpr.plt <- ggplot(rem.dat, aes(y = shannon.rxnExpr, x = Remission))+
    geom_boxplot()+
    annotate(geom = "text", x = 0.5,y=max(rem.dat[,shannon.rxnExpr]) , label = paste("p =", round(coef(summary(rxnExpr.mod))[2,4], digits = 3)), hjust = 0, vjust = 0)+
    ggtitle("RxnExpr Shannon")+
    theme_bw()

FVArange.mod <- glmer(data = rem.dat,
                   family = "binomial",
                   formula = Remission ~ HB_Mayo_impu +shannon.FVArange + (1|PatientID))
FVArange.plt <- ggplot(rem.dat, aes(y = shannon.FVArange, x = Remission))+
    geom_boxplot()+
    annotate(geom = "text", x = 0.5,y=max(rem.dat[,shannon.FVArange]) , label = paste("p =", round(coef(summary(FVArange.mod))[3,4], digits = 3)), hjust = 0, vjust = 0)+
    ggtitle("FVA range Shannon")+
    theme_bw()

FVAcenter.mod <- glmer(data = rem.dat,
                   family = "binomial",
                   formula = Remission ~ shannon.FVAcenter + HB_Mayo_impu + (1|PatientID))
FVAcenter.plt <- ggplot(rem.dat, aes(y = shannon.FVAcenter, x = Remission))+
    geom_boxplot()+
    annotate(geom = "text", x = 0.5,y=max(rem.dat[,shannon.FVAcenter]) , label = paste("p =", round(coef(summary(FVAcenter.mod))[2,4], digits = 3)), hjust = 0, vjust = 0)+
    ggtitle("FVA center Shannon")+
    theme_bw()

div.plt <- plot_grid(rxnExpr.plt, pa.plt, FVArange.plt, FVAcenter.plt)
ggsave(div.plt, file = "../results/diversityPlots.Remission.pdf", width = 8, height = 8)

# Response
rem.dat <- meta[Response != "C",]
rem.dat[, Response := factor(Response)]
pa.mod <- glmer(data = rem.dat, 
                family = "binomial",
                formula = Response ~ richness.PA + HB_Mayo_impu + (1|PatientID))
pa.plt <- ggplot(rem.dat, aes(y = richness.PA, x = Response))+
    geom_boxplot() +
    annotate(geom = "text", x = 0.5,y=max(rem.dat[,richness.PA]) , label = paste("p =", round(coef(summary(pa.mod))[2,4], digits = 3)), hjust = 0, vjust = 0)+
    ggtitle("PA richness")+
    theme_bw()

rxnExpr.mod <- glmer(data = rem.dat,
                family = "binomial",
               formula = Response ~ shannon.rxnExpr + HB_Mayo_impu+(1|PatientID))
rxnExpr.plt <- ggplot(rem.dat, aes(y = shannon.rxnExpr, x = Response))+
    geom_boxplot()+
    annotate(geom = "text", x = 0.5,y=max(rem.dat[,shannon.rxnExpr]) , label = paste("p =", round(coef(summary(rxnExpr.mod))[2,4], digits = 3)), hjust = 0, vjust = 0)+
    ggtitle("RxnExpr Shannon")+
    theme_bw()

FVArange.mod <- glmer(data = rem.dat,
                   family = "binomial",
                   formula = Response ~ HB_Mayo_impu + shannon.FVArange + (1|PatientID))
FVArange.plt <- ggplot(rem.dat, aes(y = shannon.FVArange, x = Response))+
    geom_boxplot()+
    annotate(geom = "text", x = 0.5,y=max(rem.dat[,shannon.FVArange]) , label = paste("p =", round(coef(summary(FVArange.mod))[3,4], digits = 3)), hjust = 0, vjust = 0)+
    ggtitle("FVA range Shannon")+
    theme_bw()

FVAcenter.mod <- glmer(data = rem.dat,
                   family = "binomial",
                   formula = Response ~ HB_Mayo_impu + shannon.FVAcenter + (1|PatientID))
FVAcenter.plt <- ggplot(rem.dat, aes(y = shannon.FVAcenter, x = Response))+
    geom_boxplot()+
    annotate(geom = "text", x = 0.5,y=max(rem.dat[,shannon.FVAcenter]) , label = paste("p =", round(coef(summary(FVAcenter.mod))[2,4], digits = 3)), hjust = 0, vjust = 0)+
    ggtitle("FVA center Shannon")+
    theme_bw()

div.plt <- plot_grid(rxnExpr.plt, pa.plt, FVArange.plt, FVAcenter.plt)
ggsave(div.plt, file = "../results/diversityPlots.Response.pdf", width = 8, height = 8)


# plot the size of the models
HB.mod <- lmer(data = meta,
               formula = HB_Mayo_impu ~ modelSize + (1|PatientID))
modsize.HB <- ggplot(meta, aes(x = HB_Mayo_impu, y = modelSize)) +
    geom_point() +
    geom_smooth(method = "lm") +
    annotate(geom = "text", x = min(meta[,HB_Mayo_impu]),y=max(meta[,modelSize]) , label = paste("p =", round(coef(summary(HB.mod))[2,5], digits = 3)), hjust = 0, vjust = 0)+
    theme_bw()

rem.dat <- meta[Remission != "C",]
rem.dat[, Remission := factor(Remission)]
Rem.mod <- glmer(data = rem.dat,
                 family = "binomial",
                 formula = Remission ~ modelSize + HB_Mayo_impu + (1|PatientID))
modsize.Rem <- ggplot(rem.dat, aes(x = Remission, y = modelSize)) +
    geom_boxplot() +
    annotate(geom = "text", x = 0.5,y=max(rem.dat[,modelSize]) , label = paste("p =", round(coef(summary(Rem.mod))[2,4], digits = 3)), hjust = 0, vjust = 0)+
    theme_bw()
                 
rem.dat <- meta[Response != "C",]
rem.dat[, Response := factor(Response)]
Res.mod <- glmer(data = rem.dat,
                 family = "binomial",
                 formula = Response ~ modelSize + HB_Mayo_impu + (1|PatientID))

modsize.Res <- ggplot(rem.dat, aes(x = Response, y = modelSize)) +
    geom_boxplot() +
    annotate(geom = "text", x = 0.5,y=max(rem.dat[,modelSize]) , label = paste("p =", round(coef(summary(Res.mod))[2,4], digits = 3)), hjust = 0, vjust = 0)+
    theme_bw()
modSize.plt <- plot_grid(modsize.HB, modsize.Res, modsize.Rem, nrow = 1)
ggsave(modSize.plt, file = "../results/modsizePlots.pdf", width = 8, height = 4)
