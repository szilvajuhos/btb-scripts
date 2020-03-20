# btb-scripts
Collection of ad-hoc scripts for BTB

## addStrelkaAFs.py  
Calculates *somatic* SNV and indel allele frequencies from Strelka calls. To get a histogram with bcftools and R use like

```
python3 addStrelkaAFs.py -v testdata/Strelka_TUMOR_vs_NORMAL_somatic_snvs_snpEff_VEP.ann.vcf | \
bcftools query -f '%SAF\n'  | \
Rscript -e "png(\"snvs_SAF_histogram.png\"); vcf<-read.csv(\"/dev/stdin\"); hist(vcf[,1])"
```



## gatkVQSR.sh
Applies Variant Quality Score Recalibration for HaplotypeCaller VCF files using singularity container from nf-core/sarek

