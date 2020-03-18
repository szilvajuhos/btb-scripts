#!/bin/bash -x
GRCH38=/data1/references/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
REFDIR=/data1/references/annotations/GATK_bundle
HAPMAP=${REFDIR}/hapmap_3.3.hg38.vcf.gz                             # training with high-confidence TPs
OMNI1000G=${REFDIR}/1000G_omni2.5.hg38.vcf.gz                       # training with some FPs
PHASE11000G=${REFDIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz
DBSNP146=${REFDIR}/dbsnp_146.hg38.vcf.gz                            # validation
LOCALBASE=`basename $1`
VQSROUT=${LOCALBASE%.vcf*}

echo "######################################################################"
echo "                        SNP recalibration"
echo "######################################################################"

singularity exec /data1/containers/sarek-latest.simg gatk VariantRecalibrator \
  -R ${GRCH38} \
  -V $1 \
  --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${HAPMAP} \
  --resource:omni,known=false,training=true,truth=false,prior=12.0 ${OMNI1000G} \
  --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${PHASE11000G} \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP146} \
  -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
  -mode SNP \
  -O ${VQSROUT}.VQSR.SNP.recal.vcf \
  --tranches-file ${VQSROUT}.VQSR.SNP.tranches \
  --rscript-file ${VQSROUT}.SNP.plots.R      # needs r-ggplot2 conda package

singularity exec /data1/containers/sarek-latest.simg gatk ApplyVQSR \
  -R ${GRCH38} \
  -V $1 \
  -O ${VQSROUT}.SNP.recalibrated.vcf.gz \
  --truth-sensitivity-filter-level 99.0 \
  --tranches-file  ${VQSROUT}.VQSR.SNP.tranches \
  --recal-file ${VQSROUT}.VQSR.SNP.recal.vcf \
  -mode SNP


echo "######################################################################"
echo "                        indel recalibration"
echo "######################################################################"


MILLSINDELS=${REFDIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz

singularity exec /data1/containers/sarek-latest.simg gatk VariantRecalibrator \
  -R ${GRCH38} \
  -V ${VQSROUT}.SNP.recalibrated.vcf.gz \
  --resource:mills,known=false,training=true,truth=true,prior=12.0 ${MILLSINDELS} \
  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${DBSNP146} \
  -an DP -an FS -an QD -an SOR -an MQRankSum -an ReadPosRankSum \
  -mode INDEL \
  -O ${VQSROUT}.VQSR.INDEL.recal.vcf \
  --tranches-file ${VQSROUT}.VQSR.INDEL.tranches \
  --rscript-file ${VQSROUT}.INDEL.plots.R

singularity exec /data1/containers/sarek-latest.simg gatk ApplyVQSR \
  -R ${GRCH38} \
  -V ${VQSROUT}.SNP.recalibrated.vcf.gz \
  -O ${VQSROUT}.SNP.INDEL.recalibrated.vcf.gz \
  --truth-sensitivity-filter-level 99.0 \
  --tranches-file  ${VQSROUT}.VQSR.INDEL.tranches \
  --recal-file ${VQSROUT}.VQSR.INDEL.recal.vcf \
  -mode INDEL


