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
   -O ${VQSROUT}.VQSR.recal.vcf \
   --tranches-file ${VQSROUT}.VQSR.tranches \
   --rscript-file ${VQSROUT}.plots.R      # needs r-ggplot2 conda package


singularity exec /data1/containers/sarek-latest.simg gatk ApplyVQSR \
  -R ${GRCH38} \
  -V $1 \
  -O ${VQSROUT}.recalibrated.vcf.gz \
  --truth-sensitivity-filter-level 99.0 \
  --tranches-file  ${VQSROUT}.VQSR.tranches \
  --recal-file ${VQSROUT}.VQSR.recal.vcf \
  -mode SNP


echo "######################################################################"
echo "                        SNP recalibration"
echo "######################################################################"



