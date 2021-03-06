#!/bin/bash
set -euo pipefail 

# removing old data
rm *csv
for f in `cat current.samples`; do 
	echo $f; 
	for c in `seq 1 22` X Y;do 
		CHROM=chr${c}; 
		echo -n $f",">>${CHROM}.csv; 
		# making new data
		python it.py -t ratio -f ratios/${f}*ratio.txt -c $CHROM -i ~/genome/Homo_sapiens_assembly38.fasta.fai >> ${CHROM}.csv;
	done; 
done

python tight.py
