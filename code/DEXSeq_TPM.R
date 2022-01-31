d <- read.table('results/tables/DEXSeq/DEXSeq_results.csv', header = T, sep = "\t")

read_counts <- d[,grep('countData.', colnames(d))]
head(read_counts)

gene_length <- d$genomicData.width
head(gene_length)

condition <- c(rep('CMV', 4), rep('MYO', 4))
  
read_counts_length <- read_counts/gene_length
tpm <- t(t(read_counts_length)*1e6/colSums(read_counts_length))
colnames(tpm) <- paste(colnames(tpm), 'tpm', sep = '_')
# calculate means
mean_tpm <- cbind(rowMeans(tpm[,c(1:4)]), rowMeans(tpm[,c(5:8)]))
colnames(mean_tpm) <- paste(c('CMV', 'MYO'), 'tpm', sep = '_')
# merge
tpm_d <- cbind(d, mean_tpm, tpm)
# move columns
lastColumns <- c('exon_id', 'transcripts')
tpm_d <- tpm_d[,c(setdiff(colnames(tpm_d), lastColumns), lastColumns)]

write.table(tpm_d, 'results/tables/DEXSeq/DEXSeq_results_TPM.csv',
            row.names = F, sep = '\t', quote = F)
