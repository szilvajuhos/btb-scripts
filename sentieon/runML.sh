#!/bin/bash -x
# skeleton script to run TNScope with machine-learning 
# $1 -- TUMOR BAM
# $2 -- NORMAL BAM
#
# We are expecting the recalibration tables being at the same directory 
# as the corresponding BAM, having a name like ${file%.bam}recal.table
# i.e. like tumor.bam , tumor_anything_recal.table

module load sentieon
module load samtools
module load bcftools

REFERENCE=/data1/references/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
ML_MODEL=/data1/software/sentieon/SentieonTNscopeModel_GiAB_HighAF_LowFP-201711.05.model
ML_THRESHOL=0.81

TBASE=`basename $1`
TDIR=`dirname $1`
TUMOR_DEDUPED_BAM=$1
TUMOR_RECAL_DATA_TABLE=${1%*bam}*recal.table
# sample name
TUMOR=`samtools view -H ${TUMOR_DEDUPED_BAM}| awk '/^@RG/{print $4}'| uniq| cut -b 4-`

NBASE=`basename $2`
NDIR=`dirname $2`
NORMAL_DEDUPED_BAM=$2
NORMAL_RECAL_DATA_TABLE=${2%.bam}*recal.table
# sample name
NORMAL=`samtools view -H ${NORMAL_DEDUPED_BAM}| awk '/^@RG/{print $4}'| uniq| cut -b 4-`

TMP_VARIANT_VCF=${TUMOR}_TMP_VARIANT_VCF
VARIANT_VCF=TNScope_${TUMOR}_vs_${NORMAL}_raw.vcf
FILTER_VARIANT_VCF=${VARIANT_VCF%_raw*}_filtered.vcf

sentieon driver -t 46 -r $REFERENCE \
  -i $TUMOR_DEDUPED_BAM -q $TUMOR_RECAL_DATA_TABLE \
  -i $NORMAL_DEDUPED_BAM -q $NORMAL_RECAL_DATA_TABLE \
  --algo TNscope --tumor_sample $TUMOR --normal_sample $NORMAL \
  --clip_by_minbq 1 --max_error_per_read 3 \
  --min_init_tumor_lod 2.0 --min_base_qual 10 --min_base_qual_asm 10 \
  --min_tumor_allele_frac 0.00005 $TMP_VARIANT_VCF

sentieon driver -t 46 -r $REFERENCE --algo TNModelApply \
  --model $ML_MODEL -v $TMP_VARIANT_VCF $VARIANT_VCF

bcftools filter -s "ML_FAIL" -i "INFO/ML_PROB > $ML_THRESHOLD" $VARIANT_VCF \
  -O z -m x -o $FILTER_VARIANT_VCF

