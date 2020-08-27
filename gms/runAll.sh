#nextflow -log ${1%.tsv}.sentieon.log run nf-core/sarek -r dev -profile munin --tools TNscope,merge --sentieon --input $1 

INPUT=results/Preprocessing/TSV/sentieon_recalibrated.tsv

for vc in haplotypecaller manta mutect2 strelka; do  
  nextflow -log ${1%.tsv}.${vc}.log run nf-core/sarek -r dev -profile munin --step variantcalling --tools ${vc} --input $INPUT --no_gvcf
done

nextflow -log ${1%.tsv}.CNV.log run nf-core/sarek -r dev -profile munin --step variantcalling --tools controlfreec,ascat --input $INPUT
