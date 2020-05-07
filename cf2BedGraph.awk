
# massage input file like
# awk '!/Chrom/{print}' P13714_154.pileup_BAF.txt| sort -V -k1,1 -k2,2n 

BEGIN{
  getline;  # have to start with getline to get the FILENAME
  print "track type=bedGraph name=\""FILENAME" BAF\" description=\""FILENAME" BAF\" color=100,100,100"
  chr=$1;
  start=$2;
  baf=$3;
}
# it is like
# 1       10181   0.0495495       -1      -1      -1      -1      -1
# 1       10231   0.075   -1      -1      -1      -1      -1
# 1       10257   0.130178        -1      -1      -1      -1      -1
#
# and we want to have
#
# track type=bedGraph name="P13714_154.pileup_BAF.txt BAF" description="P13714_154.pileup_BAF.txt BAF" color=100,100,100
# 1   10181   10230 0.0495495
# 1   10231   10256 0.075
#
{
  if(chr==$1) {
    printf("%s\t%d\t%d\t%f\n",chr,start,$2-1,baf);
  } 
  # at change of chromosomes we are only updating
  chr=$1;
  start=$2;
  baf=$3
}
