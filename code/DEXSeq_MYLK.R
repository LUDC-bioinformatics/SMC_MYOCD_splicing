library(reshape2)

d <- read.table('results/tables/DEXSeq/DEXSeq_results_TPM.csv', header = T, sep = "\t")

g <- d[d$groupID=="ENSG00000065534.19_MYLK",]

#exon32 <- 123612298-123614349

mean_gene_exp_CMV_MYO <- c(14.2430957267229, 58.7654231920079)

g_tpm <- g[,c('featureID', 'countData.CMV1_tpm', 'countData.CMV2_tpm', 'countData.CMV4_tpm', 'countData.CMV5_tpm',
              'countData.MYO1_tpm', 'countData.MYO2_tpm', 'countData.MYO4_tpm', 'countData.MYO5_tpm')]
rownames(g_tpm) <- g$featureID
countsNorm <- log2(as.matrix(g_tpm)+1)
colnames(countsNorm) <- gsub('_tpm|countData.', '', colnames(countsNorm))
# melt
norm_melt <- melt(countsNorm)
colnames(norm_melt) <- c('bin', 'samplename', 'normalized_counts')
# add annonation
meta_exp_levels <- data.frame(condition=c(rep('CMV', 4), rep('MYO', 4)), samplename=c('CMV1', 'CMV2', 'CMV4', 'CMV5',
                                                                                      'MYO1', 'MYO2', 'MYO4', 'MYO5'))
top_genes_norm <- merge(norm_melt, meta_exp_levels, by='samplename')
# reverse order
top_genes_norm <- top_genes_norm[order(top_genes_norm$bin, decreasing = T),]
# plot
ggplot(top_genes_norm) +
  geom_point(aes(x = bin,
                 y = normalized_counts,
                 color = condition,
                 shape=condition),
             position=position_jitter(w=0.2,h=0)) +
  xlab("Genes") +
  ylab("log2(TPM counts)") +
  ggtitle('Exon bins of MYLK') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))


write.table(g, 'results/tables/DEXSeq/DEXSeq_results_MYLK.csv',
            row.names = F, sep = '\t', quote = F)
