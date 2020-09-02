#!/bin/bash -x
# First start NGB server on 8080 like
# java -jar catgenome-latest.jar


# the CLI interface
NGB=/home/szilva/sr/NGB/ngb-cli/bin/ngb

# add FASTA reference
${NGB} reg_ref ${HOME}/genome/Homo_sapiens_assembly38.fasta -n GRCh38 -t

# add corresponding annotation
#${NGB} reg_file GRCh38 ${HOME}/genome/gencode.v28.chr_patch_hapl_scaff.annotation.sorted.gtf -n GRCh38_genes -t
${NGB} reg_file GRCh38 ${HOME}/genome/gencode.v29.basic.annotation.sorted.gtf -n GRCh38_genes -t

# associate the two 
${NGB} add_genes GRCh38 GRCh38_genes

# Add datasets is like:

# - register dataset and associate with the reference and annotation
#${NGB} reg_dataset GRCh38 TEST_SET GRCh38_genes -t
#${NGB} add	TEST_SET	sample.bam sample.vcf

