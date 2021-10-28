library(tximport)

samples <- c("CMV1", "CMV2", "CMV4", "CMV5", "MYO1", "MYO2", "MYO4", "MYO5")
files <- file.path('results/tables/salmon/', samples, 'quant.sf.gz')
names(files) <- samples
t2g <- read.table('intermediate/transcript2gene.csv', header = F, sep = "\t")
read_counts <- tximport(files, type="salmon", txOut=TRUE, countsFromAbundance="dtuScaledTPM", tx2gene=t2g)
countData <- read_counts$counts
head(countData)

dim(countData)[1]
keep <- which(rowSums(countData)>0)
length(keep)
countData <- countData[keep,]

countDataCMV <- countData[,c("CMV1", "CMV2", "CMV4", "CMV5")]
countDataMYO <- countData[,c("MYO1", "MYO2", "MYO4", "MYO5")]

write.table(countData, 'results/tables/salmon/salmon_TPM.csv', row.names = T, sep = '\t', quote = F)
write.table(countDataCMV, 'results/tables/salmon/salmon_TPM_CMV.csv', row.names = T, sep = '\t', quote = F)
write.table(countDataMYO, 'results/tables/salmon/salmon_TPM_MYO.csv', row.names = T, sep = '\t', quote = F)
