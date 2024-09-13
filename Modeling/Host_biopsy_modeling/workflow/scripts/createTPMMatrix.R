# Porthmeus
# 27.10.20

# convert read count matrix with ensemble identifier as row names into a tpm matrix

# log the errors and warnings
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

require(GenomicFeatures)

# load data and remove NAs, correct column/sample names
df_RC <-read.csv(snakemake@input[["RC"]], 
           header=TRUE, 
           row.names = 1
)
df_RC[is.na(df_RC)] <- 0
colnames(df_RC) <- gsub("^X", "", colnames(df_RC))


# load gff data from ensembl
print("loading data from the ENSEMBL, this might take some minutes...")
txdb_GRCh38 <- makeTxDbFromEnsembl(organism = snakemake@params[["organism"]], 
                                    server = "ensembldb.ensembl.org")

# calculate effective length of each gene
print("calculating effective lengths")
exons <- exonsBy(txdb_GRCh38, by = "gene")
effLen <- sum(width(reduce(exons)))
effLen <- effLen/1000

# calculate TPM and RPKM

# remove all genes present in the RC matrix but not present in the human genome
# (wherever these come from) and removeall genes from the effLen vector which
# are not present in the read count matrix (wherever these are gone to)
effLen <- effLen[names(effLen) %in% rownames(df_RC)]
df_RC <- df_RC[names(effLen),]


print("calculate FPKM")
df_FPKM <- t(t(df_RC)/(apply(df_RC, 2, function(x){
                            sum(x, na.rm = TRUE)})))
df_FPKM <- df_FPKM/effLen
df_FPKM <- df_FPKM * 1e6
df_FPKM <- as.data.frame(df_FPKM)

print("calculate_TPM")
df_TPM <- df_RC/(effLen)
sizeFactor <- apply(df_TPM,2, sum, na.rm =TRUE)
df_TPM <- t(t(df_TPM)/sizeFactor)
df_TPM <- df_TPM*1e6
df_TPM <- as.data.frame(df_TPM)


# write the data frames

write.csv(df_TPM, file = snakemake@output[["TPM"]])
write.csv(df_FPKM, file = snakemake@output[["FPKM"]])
