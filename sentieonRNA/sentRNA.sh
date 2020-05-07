#!/bin/bash -ex
# *******************************************
# Script to perform RNA variant calling
# using a single sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# *******************************************

# Update with the fullpath location of your sample fastq
fastq_1=$1
fastq_2=$2
group=`basename $1 _R1_001.fastq.gz`
sample=`echo $group| awk -F_ '{print $1"_"$2}'`
platform="ILLUMINA"

# Update with the location of the reference data files
REFBASE=/data2/fusiontest/GRCh38
#fasta=${REFBASE}/Homo_sapiens.GRCh38.dna.toplevel.fa
fasta=/data1/references/igenomes/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
#dbsnp=${REFBASE}/GCF_000001405.38.noCHR.vcf.gz
dbsnp=/data1/references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz
known_Mills_indels=/data1/references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
known_1000G_indels=/data1/references/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz
#comment if STAR should generate a new genomeDir
star_fasta=${REFBASE} #/genomeDir
#uncomment if you would like to use a bed file
#interval_file="/home/regression/references/b37/TruSeq_exome_targeted_regions.b37.bed"

# Update with the location of the Sentieon software package and license file
export SENTIEON_INSTALL_DIR=/data1/software/sentieon/sentieon-genomics-201911/
export SENTIEON_LICENSE=/data1/software/sentieon/sentieon.lic
star_binary="singularity exec /data1/containers/nfcore-rnafusion-dev.img /opt/conda/envs/nf-core-rnafusion-dev/bin/STAR"

# Other settings
nt=48 #number of threads to use in computation, set to number of cores in the server
workdir="$PWD/${sample}" #Determine where the output files will be stored

echo "******************************************"
echo "0. Setup"
echo "******************************************"
mkdir -p $workdir
#logfile=$workdir/run.log
#exec >$logfile 2>&1
cd $workdir
if [ ! -z "$interval_file" ]; then
  driver_interval_option="--interval $interval_file"
  realign_interval_option="--interval_list $interval_file"
fi

# ******************************************
# 1. Mapping reads with STAR
# ******************************************
if [ -z "$star_fasta" ]; then
  star_fasta="genomeDir"
  # The genomeDir generation could be reused
  mkdir $star_fasta
  $star_binary --runMode genomeGenerate --genomeDir $star_fasta --genomeFastaFiles $fasta --runThreadN $nt
fi
#perform the actual alignment and sorting
$star_binary --twopassMode Basic --genomeDir $star_fasta --runThreadN $nt --outSAMtype BAM SortedByCoordinate --twopass1readsN -1 --sjdbOverhang 150 --readFilesIn $fastq_folder/$fastq_1 $fastq_folder/$fastq_2 --readFilesCommand zcat --outSAMattrRGline ID:$group SM:$sample PL:$platform
mv Aligned.sortedByCoord.out.bam sorted.bam
$SENTIEON_INSTALL_DIR/bin/sentieon util index sorted.bam

# ******************************************
# 2. Metrics
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver $driver_interval_option -r $fasta -t $nt -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o gc-report.pdf gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qd-report.pdf qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o mq-report.pdf mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o is-report.pdf is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads. It is possible
# to mark instead of remove duplicates
# by ommiting the --rmdup option in Dedup
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i sorted.bam --algo LocusCollector --fun score_info score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt deduped.bam 

# ******************************************
# 2a. Coverage metrics
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --algo CoverageMetrics coverage_metrics

# ******************************************
# 4. Split reads at Junction
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --algo RNASplitReadsAtJunction --reassign_mapq 255:60 splitted.bam

# ******************************************
# 6. Base recalibration
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver $driver_interval_option -r $fasta -t $nt -i splitted.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table
$SENTIEON_INSTALL_DIR/bin/sentieon driver $driver_interval_option -r $fasta -t $nt -i splitted.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o recal_plots.pdf recal.csv

# ******************************************
# 7. HC Variant caller for RNA
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver $driver_interval_option -r $fasta -t $nt -i splitted.bam -q recal_data.table --algo Haplotyper -d $dbsnp --trim_soft_clip --emit_conf=20 --call_conf=20 output-hc-rna.vcf.gz

