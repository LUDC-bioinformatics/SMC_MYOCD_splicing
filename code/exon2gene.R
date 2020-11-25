# Rscript exon2gene.R {input_d} {input_exon_gene} {output}

ff <- commandArgs(trailingOnly = T)

d <- read.table(ff[1], header = T)
# head(d)

exon_gene <- read.table(ff[2])
colnames(exon_gene) <- c('gene', 'exon')
exon_gene$gene <- gsub("\\..*","", exon_gene$gene)
exon_geneUniq <- unique(exon_gene)
# head(exon_geneUniq)

dgene <- merge(d, exon_geneUniq, by = "exon")
# dim(d)
# dim(dgene)

dgene <- dgene[,c(which(colnames(dgene)=="gene"),which(colnames(dgene)!="gene"))]
names(dgene) <- gsub("X","", names(dgene))
dgene$exon <- gsub("\\..*","", dgene$exon)
# sum(duplicated(dgene$exon))
# head(dgene, 2)

write.table(dgene, ff[3],
            row.names = F, sep = "\t", quote = F)
