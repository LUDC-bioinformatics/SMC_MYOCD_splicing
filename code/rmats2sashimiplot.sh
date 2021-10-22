#!/bin/bash

mkdir intermediate/rmats2sashimiplot
mkdir results/figures/rmats2sashimiplot

for i in SE A5SS A3SS MXE RI;
  do
    rmats2sashimiplot \
      --b1 intermediate/STAR/CMV1_pass2/CMV1.Aligned.out.sorted.bam,intermediate/STAR/CMV2_pass2/CMV2.Aligned.out.sorted.bam,intermediate/STAR/CMV4_pass2/CMV4.Aligned.out.sorted.bam,intermediate/STAR/CMV5_pass2/CMV5.Aligned.out.sorted.bam \
      --b2 intermediate/STAR/MYO1_pass2/MYO1.Aligned.out.sorted.bam,intermediate/STAR/MYO2_pass2/MYO2.Aligned.out.sorted.bam,intermediate/STAR/MYO4_pass2/MYO4.Aligned.out.sorted.bam,intermediate/STAR/MYO5_pass2/MYO5.Aligned.out.sorted.bam \
      -t $i \
      -e results/tables/rMATS_bam/${i}.MATS.JC.sign.txt \
      --l1 CMV \
      --l2 MYO \
      -o intermediate/rmats2sashimiplot/rmats2sashimiplot_${i} \
      --group-info conf/grouping.gf
      
      mv intermediate/rmats2sashimiplot/${i}/Sashimi_plot results/figures/rmats2sashimiplot/${i};
  done
