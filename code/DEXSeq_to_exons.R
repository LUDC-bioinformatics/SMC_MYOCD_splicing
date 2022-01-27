library(DEXSeq)
library(GenomicFeatures)
library(GenomicRanges)

setwd('/home/dmytro/Science/myocardin_splicing_Ola_Karl/')

res <- read.table('results/tables/DEXSeq/DEXSeq_results.csv', header = T, sep = '\t')
transcriptDb <- makeTxDbFromGFF("data/reference/gencode.v35.primary_assembly.annotation.gtf", format="gtf")

exonBins <- GRanges(
  seqnames=res$genomicData.seqnames, 
  ranges=IRanges( 
    res$genomicData.start, res$genomicData.end),    
  strand=res$genomicData.strand,
  mcols=res[,c('groupID', 'featureID')])

exonsByTranscript <- unlist(exonsBy(transcriptDb, "tx", use.names=TRUE))

bins_exons_overlap <- as.data.frame(findOverlaps(exonsByTranscript, exonBins))

# re-order by hits order
exons_ordered <- as.data.frame(exonsByTranscript[bins_exons_overlap$queryHits,], row.names = NULL)
bins_ordered <- as.data.frame(exonBins[bins_exons_overlap$subjectHits,])

# combine bins and exons with some name corrections
exons_ordered$exon_id <- NULL
colnames(exons_ordered) <- c("exon_chr", "exon_start", "exon_end",
                             "exon_width","exon_strand",
                             "exon_id","exon_number")
bins_exons <- cbind(bins_ordered[,c("mcols.groupID", "mcols.featureID")], exons_ordered)
colnames(bins_exons)[c(1,2)] <- c('DEXSeq_groupID', 'DEXSeq_featureID')
bins_exons <- bins_exons[order(bins_exons$exon_chr, bins_exons$exon_start, bins_exons$exon_end),]
write.table(bins_exons,'results/tables/DEXSeq/DEXSeq_results_bins2exons.csv',
            row.names = F, sep = '\t')

# add exon_id to the results
bins_exons_aggregate <- aggregate(bins_exons$exon_id, list(bins_exons$mcols.groupID, bins_exons$mcols.featureID),
                                  paste, collapse=",")
colnames(bins_exons_aggregate)[3] <- "exon_id"
res_bins_exons <- merge(res, bins_exons_aggregate, by=c(1,2))
# remove `_exon_id` from file name. it is added for over-write protection:
write.table(res_bins_exons,'results/tables/DEXSeq/DEXSeq_results_exon_id.csv', 
            row.names = F, sep = '\t')
