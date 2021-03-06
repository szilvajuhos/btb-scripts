#!/bin/bash
set -e
set -o pipefail
set -u

#nextflow -log ${1%.tsv}.sentieon.log run nf-core/sarek -r dev -profile munin --tools TNscope,merge --sentieon --input $1 
#
#INPUT=results/Preprocessing/TSV/sentieon_recalibrated.tsv
#
#for vc in haplotypecaller manta mutect2 strelka; do  
#  nextflow -log ${1%.tsv}.${vc}.log run nf-core/sarek -r dev -profile munin --step variantcalling --tools ${vc} --input $INPUT --no_gvcf
#done
#
#nextflow -log ${1%.tsv}.CNV.log run nf-core/sarek -r dev -profile munin --step variantcalling --tools controlfreec,ascat --input $INPUT
#
#
## Annotation
## Manta diploid and somatic calls
#nextflow -log manta.ann.${1%.tsv}.log run nf-core/sarek -r dev -profile munin --step annotate --tools snpEff --input "results/VariantCalling/*/Manta/Manta_*{diploid,somatic}SV.vcf.gz"
#
## Strelka somatic
#nextflow -log strelka.ann.${1%.tsv}.log run nf-core/sarek -r dev -profile munin --step annotate --tools merge --input "results/VariantCalling/*/Strelka/*somatic_{indels,snvs}.vcf.gz" 
#
## Mutect2 annotation
#nextflow -log mutect2.ann.${1%.tsv}.log run nf-core/sarek -r dev -profile munin --step annotate --tools merge --input "results/VariantCalling/*/Mutect2/*vcf.gz"
#
## HaplotypeCaller
#nextflow -log haplotypecaller.ann.${1%.tsv}.log run nf-core/sarek -r dev -profile munin --step annotate --tools merge --input "results/VariantCalling/*/HaplotypeCaller/*vcf.gz"

# now do annotation:
# add CAF, TOPMED, dbSNP and SweGEn allele frequencies
module load vcfanno
PFDIR=results/PF
mkdir -p results/PF
date
for f in `find results/Annotation/ -name "*VEP*vcf.gz"`; do
  BASE=`basename -s vcf.gz $f`
  echo "Annotating $f"
  vcfanno ~/dev/btb-scripts/gms/3AFs.toml $f > ${PFDIR}/${BASE}AF.vcf &
done
wait
date

