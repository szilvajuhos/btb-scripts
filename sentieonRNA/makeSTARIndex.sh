#!/bin/bash -x
CONTAINER=/data1/containers/nfcore-rnafusion-dev.img
REFERENCE=ENSEMBL/Homo_sapiens.GRCh38.dna.toplevel.fa
GTF=ENSEMBL/Homo_sapiens.GRCh38.99.chr_patch_hapl_scaff.gtf
# calculate read length
READLENGTH=`zcat $1 | head -2| awk '{getline;print length($1)-1}'`

singularity exec ${CONTAINER} \
  STAR --runThreadN 48 \
  --runMode genomeGenerate \
  --genomeDir ./ENSEMBL \
  --genomeFastaFiles ${REFERENCE} \
  --sjdbGTFfile ${GTF} \
  --limitGenomeGenerateRAM 256000000000 \
  --sjdbOverhang ${READLENGTH}
