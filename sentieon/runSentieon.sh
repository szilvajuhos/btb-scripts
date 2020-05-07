#!/bin/bash -x
# simple nextflow run to use sentieon
# $1 -- input TSV
# $2 -- additional things (i.e. -resume)

#nextflow run nf-core/sarek -r 2.5.2 -profile munin --input ../tsv/P2233_101T_P2233_120N.tsv --sentieon --tools DNAseq,DNAscope,TNscope,ASCAT,ControlFREEC,Mutect2,Manta,Strelka,HaplotypeCaller --skip_qc --pon /data1/PON/vcfs/BTB.PON.vcf.gz --pon_index /data1/PON/vcfs/BTB.PON.vcf.gz.tbi $1

nextflow run nf-core/sarek -r dev -profile munin --step mapping --input $INPUT --sentieon --skip_qc all $1
