---
title: "Myocardin splicing"
author:
   name: "Dmytro Kryvokhyzha"
   email: dmytro.kryvokhyzha@med.lu.se
   affiliation: LUDC Bioinformatics Unit
date: "14 januari, 2021"
output:
  html_document:
    keep_md: true
---



## PI

Name: [Ola Hansson](https://www.ludc.lu.se/ola-hansson-pi) & 
      [Karl Swärd](https://www.lunduniversity.lu.se/lucat/user/d6a67258cd7d2448a83f44c29cbda4e8)

Email: [ola.hansson@med.lu.se](mailto:ola.hansson@med.lu.se), 
       [karl.sward@med.lu.se](mailto:karl.sward@med.lu.se)


## Project

Smooth muscle cells from the human coronary artery were treated either with *Ad-CMV-null* or *Ad-CMV-MYOCD*, i.e. adenovirsuses that express nothing or myocardin under control of CMV promoter. There were 4 virus induced and 4 control samples. Timepoint: ~ 96h or 120h or 140h. They were RNA-Seq sequenced using Stranded Total RNA Library Prep kit.

The raw sequencing data is stored in `201022_NB501805_0127_AHY7KCBGXG`.

We want to **test whether myocarding is a splicing regulator**. To that end, we look for
splicing variation in several genes that are over-expressed by viral transduction of
myocarding.

- CMV - control
- MYO - virus induced.

There is no need to control for confounding factors because the parental cells are
split into wells that get either null or *MYOCD*. 

Expression of myocardin (*MYOCD*) can be used as control.

As positive control genes verifying that *MYOCD* over expression worked we can look at *MYH11*, *MCAM* and *SYNM* for example. In addition there are hundreds of genes that respond well, but not as much as those above. It would be important to know if much less responsive target genes have also moved in the experiment. Examples include *RRAS*, *RBPMS*, and *RBPMS2*. Two of these (*RBPMS* and *RBPMS2*) are the splicing factors that we believe are regulated by *MYOCD* and that effectuate downstream changes in splicing. Important to know consequently if they change.

We also have qPCR information in *MYH11, MCAM, KCNMB1*:

![](data/MYH11-MCAM-KCNMB1.png)

Differentially splicing of the following genes is of interest:

*ACTN1, SMTN, NFIC, BICD2, TEAD1, TPM1, TPM2, PDLIM7, FHL1, CSR2,
TGFB11, PKD1, MYL6, PDLIM5, FLNA, MBNL1, ZFR, FNBP1, SORBS1, PIP5K1C,
TNRC18, SVIL, KLHL42, LARGE1, MICAL3*

This is provided that *RBPMS* and *RBPMS2* are differentially expressed.

## Data

The raw data and all files of this project are stored on the Indigo cluster at LUDC:

`indigo.crc.med.lu.se:/ludc/Active_Projects/myocardin_splicing/`

## Prerequisites

You need to install [Conda](https://conda.io) and load the pre-configured conda
environment. It should also install all the required programs.


```bash
conda env create -f conf/conda.yml
conda activate myocardin
```

## Quality Control (QC)

Quality Control (QC) of the fastq files and other data are summarized in
`results/reports/multiQC`

## Map and count reads at exon and gene level

1. Map reads to the genome reference with [STAR](https://github.com/alexdobin/STAR) 
in 2-pass mode.

2. Count reads mapped to genes and exons with [featureCounts](http://bioinf.wehi.edu.au/featureCounts/).

Rules are added to the same `code/Snakefile`:


```bash
snakemake -s code/Snakefile \
   -j 100 \
   -p --use-conda \
   --cluster-config conf/cluster.yml \
   --cluster "condor_qsub -o logs/{rule}.out -e logs/{rule}.err -l procs={cluster.cores},mem={cluster.ram} -V"
```

Results:

- `results/tables/featureCounts/featureCounts_counts_gene.csv.gz` - gene counts.
- `results/tables/featureCounts/featureCounts_counts_exon.csv.gz` - exon counts data.
- `results/tables/featureCounts/featureCounts_counts_exon.csv.jcounts` - number of reads supporting each exon-exon junction.

## Differential gene expression

Performed with [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html):


```bash
R -e 'rmarkdown::render("code/DESeq.Rmd", output_dir="results/reports/")'
```

Results:

- `results/reports/DESeq.html` - report describing the analysis.

- `results/tables/DESeq/DESeq_all.xlsx` - differential expression results.

- `results/tables/DESeq/DESeq_interesting.xlsx` - differential expression results
  for provided list of potential candidate genes.

Column names in the tables:

- *baseMean* - mean of normalized counts for all samples. 
- *log2FoldChange* - log2 fold change (MLE): condition treatment vs control.
- *lfcSE* - standard error.
- *stat* - Wald statistic.
- *pvalue* - Wald test p-value.
- *padj* - BH adjusted p-values.
- *\_tpm* - transcript per million (TPM) normalized count data.

**TPM normalization**.

Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
Count up all the RPK values in a sample and divide this number by 1,000,000. This is your "per million" scaling factor.
Divide the RPK values by the "per million" scaling factor. This gives you TPM.

When you use TPM, the sum of all TPMs in each sample are the same. This 
makes it easier to compare the proportion of reads that mapped to a gene
in each sample. In contrast, with RPKM and FPKM, the sum of the 
normalized reads in each sample may be different, and this makes it 
harder to compare samples directly.

Here’s an example. If the TPM for gene A in Sample 1 is 3.33 and the TPM
in sample B is 3.33, then I know that the exact same proportion of 
total reads mapped to gene A in both samples. This is because the sum of
the TPMs in both samples always add up to the same number (so the 
denominator required to calculate the proportions is the same, 
regardless of what sample you are looking at.)


## Differential exon usage

Performed with [DEXSeq](https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html):


```bash
mkdir -p results/reports/DEXSeqReport/
R -e 'rmarkdown::render("code/DEXSeq.Rmd", output_dir="results/reports/DEXSeqReport/")'
```

Results:

- `results/reports/DEXSeqReport/DEXSeq.html` - report describing the analysis.

- `results/reports/DEXSeqReport/testForDEU.html` - DEXSeq results.

- `results/reports/DEXSeqReport/files/*.svg` - figures that can be used in a publication.

- `results/tables/DEXSeq/DEXSeq_results.csv` - results of the differential exon usage analysis.

- `results/tables/DEXSeq/genelevel_FDR.csv` - gene-level FDR values.


