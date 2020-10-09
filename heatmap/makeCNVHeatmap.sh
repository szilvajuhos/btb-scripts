#!/bin/bash
set -euo pipefail 

# removing old data
rm -rf *csv heatmap.png

#for f in P7708_105T P4551_218T; do
#for f in `cat current.samples|head `; do 
for f in `cat current.samples`; do 
	echo $f; 
	for c in `seq 1 22` X Y;do 
		CHROM=chr${c}; 
		echo -n $f",">>${CHROM}.csv; 
		python it.py -t cnvs -f CNVs/${f}*CNVs -c $CHROM -i ~/genome/Homo_sapiens_assembly38.fasta.fai -s 250000 >> ${CHROM}.csv;
	done; 
done

python tight.py
display heatmap.png
