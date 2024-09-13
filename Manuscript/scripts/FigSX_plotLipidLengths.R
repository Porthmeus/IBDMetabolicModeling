# Porthmeus
# 13.11.23

# plot the length of the differnt fatty acids vs the different phenotypes of interest.

require(data.table)
require(lme4)
require(lmerTest)
require(ggplot2)
require(cowplot)

rsrc.dir <- file.path("..","data","Metabolomics")

# load the data
met.melt <- fread(file.path(rsrc.dir, "metaAndMetabolomics.csv"))
met.melt[Diagnosis == "CD", norm_HBMayo := HB_Mayo_impu/16]
met.melt[Diagnosis == "UC", norm_HBMayo := HB_Mayo_impu/12]


# lysoPC
subs <- list("lysoPC" = "lysoPC\\.a\\.C([0-9]+)\\.[0-9]*",
             "^PC"  = "^PC\\.a[a|e]\\.C([0-9]+)\\.[0-9]*",
             "SM" = "SM.*C([0-9]+).*",
             "^Cer" = "^Cer\\.d1[6|8]\\.[0-9]\\.([0-9]+).*",
             "^HexCer" = "^HexCer\\.d1[6|8]\\.[0-9]\\.([0-9]+).*",
             "^Hex2Cer" = "^Hex2Cer\\.d1[6|8]\\.[0-9]\\.([0-9]+).*",
             "^Hex3Cer" = "^Hex3Cer\\.d1.\\.[0-9].([0-9]+).*",
             "^CE" = "^CE\\.([0-9]+).*",
             "^DG" = "^DG.*_([0-9]+).*",
             "^TG" = "^TG.*_([0-9]+).*")
dat <- list()
for(i in 1:length(subs)){
    print(names(subs)[i])
    dat[[i]] <- cbind(met.melt[grep(names(subs)[i], metabolite),], group = gsub("\\^","", names(subs)[i]))
    dat[[i]][, C:= as.integer(gsub(subs[[i]], "\\1", metabolite))]
}
dat <- do.call(rbind, dat)
# add choline
chln <- met.melt[metabolite == "Choline",]
chln[,c("group","C") := list("Choline",2)]
dat <- rbind(dat,chln)
# calculate some extra variables
dat[,C.norm := C*Concentration]
dat[,C.z := C*z.score]
dat[,HighC := C>24]



grps <- c("PC","lysoPC","SM")

# Remission 
models <- list()
for(grp in grps){
    dat.mod <- dat[Remision != "C" & group == grp,.(z.score = median(z.score,na.rm = TRUE)), by = .(PatientID, Time_seq, Remision)]
    dat.wil <- dcast(data = dat.mod,
                     Remision+PatientID ~ paste0("t",Time_seq),
                     value.var = "z.score")
    tmps <- paste0("t",dat.mod[,unique(Time_seq)])
    p.vals <- c()
    for(tm in tmps){
        w.t <- wilcox.test(data = dat.wil[!is.na(get(tm)),], get(tm)~Remision, paired = FALSE)
        p.vals <- c(p.vals, w.t$p.value)
    }
    models[[grp]] <- data.table(group = grp,
                        Time_seq = as.integer(gsub("^t","",tmps)),
                        p.value = p.vals)
}
coefs <- do.call(rbind, models)
coefs[,padj :=p.adjust(p.value, method = "BH")]

dat.plt <- dat[Remision != "C" & group %in% grps, .(z.score = median(z.score,na.rm = TRUE)), by = .(PatientID, Time_seq, Remision, group)]
rem.boxplot <-
    ggplot(dat.plt,aes(fill = Remision, x = factor(Time_seq), y = z.score)) +
    geom_boxplot(alpha = 0.7)+
    geom_text(data = coefs,inherit.aes = FALSE, aes(label = paste0("p=",round(padj, digits = 4)), x = factor(Time_seq)), y = quantile(dat.plt[,z.score],1))+
    labs(title = "Remission", x = "Time point [d]",y = "Z-score concentration", fill = "Remission")+
    facet_wrap(~group)+
    scale_fill_brewer(palette="Dark2")+
    theme_bw() 


## Response
models <- list()
for(grp in grps){
    dat.mod <- dat[Responder != "C" & group == grp,.(z.score = median(z.score,na.rm = TRUE)), by = .(PatientID, Time_seq, Responder)]
    dat.wil <- dcast(data = dat.mod,
                     Responder+PatientID ~ paste0("t",Time_seq),
                     value.var = "z.score")
    tmps <- paste0("t",dat.mod[,unique(Time_seq)])
    p.vals <- c()
    for(tm in tmps){
        w.t <- wilcox.test(data = dat.wil[!is.na(get(tm)),], get(tm)~Responder, paired = FALSE)
        p.vals <- c(p.vals, w.t$p.value)
    }
    models[[grp]] <- data.table(group = grp,
                        Time_seq = as.integer(gsub("^t","",tmps)),
                        p.value = p.vals)
}
coefs <- do.call(rbind, models)
coefs[,padj :=p.adjust(p.value, method = "BH")]


dat.plt <- dat[Responder != "C" & group %in% grps, .(z.score = median(z.score,na.rm = TRUE)), by = .(PatientID, Time_seq, Responder, group)]
resp.boxplot <-
    ggplot(dat.plt,aes(fill = Responder, x = factor(Time_seq), y = z.score)) +
    geom_boxplot(alpha = 0.7)+
    geom_text(data = coefs,inherit.aes = FALSE, aes(label = paste0("p=",round(padj, digits = 4)), x = factor(Time_seq)), y = quantile(dat.plt[,z.score],1))+
    labs(title = "Response", x = "Time point [d]",y = "Z-score concentration", fill = "Remission")+
    facet_wrap(~group)+
    scale_fill_brewer(palette="Dark2")+
    theme_bw() 


## HBMayo
models <- list()
for(grp in grps){
    dat.mod <- dat[group == grp,.(norm_HBMayo = unique(norm_HBMayo), z.score = median(z.score,na.rm = TRUE)), by = .(PatientID, Time_seq, Remision)]
    mod <- lmer(data = dat.mod,
        formula = z.score~  norm_HBMayo + (1|PatientID))
    coefs <- data.table(coef(summary(mod)), keep.rownames = TRUE, group = grp)[grep("HBMayo",rn),]
    models[[grp]] <- coefs
}
coefs <- do.call(rbind,models)
coefs[,padj := p.adjust(`Pr(>|t|)`, method = "BH")]
coefs

dat.plt <- dat[group %in% grps, .(norm_HBMayo = unique(norm_HBMayo), z.score = median(z.score,na.rm = TRUE)), by = .(PatientID, Time_seq, group, Remision)]
dat.plt[,group := factor(group, levels = grps)]
hb.plot<-ggplot(dat.plt,aes( x = norm_HBMayo, y = z.score)) + 
    #geom_boxplot(alpha = 0.7)+
    geom_point()+
    geom_text(data = coefs,
              inherit.aes = FALSE,
              hjust = 0,
              aes(label = paste0("b=", round(Estimate, digits = 4), "; p=",round(padj, digits = 4))), x = 0, y = quantile(dat.plt[,z.score],1))+
    geom_smooth(method = "lm")+
    facet_wrap(~group)+
    labs(title = "HBI/Mayo", x = "Disease activity",y = "Z-score concentration")+
    theme_bw() 


#pairwise.wilcox.test(dat[Remision != "C" & group == "PC", z.score], dat[Remision != "C" & group == "PC", paste(Remision, Time_seq)])


#a <- dcast(data = dat[group %in% c("SM","PC","lysoPC") & !is.na(z.score),], PatientID_Time+norm_HBMayo~group, value.var = "z.score", fun.aggregate = median)
#
#ggplot(a, aes(x=norm_HBMayo, y = lysoPC/SM))+
#    geom_point() +
#    geom_smooth(method = "lm")
plt.all <- plot_grid(hb.plot, resp.boxplot, rem.boxplot, ncol = 1)

ggsave(plt.all, file = "../temp/MembraneLipids.pdf", 
       width = 9,
       height = 10)

# plot choline
dat.plt <- met.melt[metabolite == "Choline",]

ggplot(dat.plt[Responder != "C",], aes(x = factor(Time_seq), y = z.score, fill = Remision))+
    geom_boxplot()

