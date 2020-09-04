import json
from intervaltree import Interval, IntervalTree

# read file
with open('test.json', 'r') as myfile:
    data=myfile.read()

# parse file
obj = json.loads(data)

# show values
#print("ENSEMBL id: " + str(obj['id']))
#print("Gene: " + str(obj['display_name']))
#print("Description: " + str(obj['description']))
#print("Assembly: " + str(obj['assembly_name']))
#print("Species: " + str(obj['species']))
#print(json.dumps(obj, indent=4, sort_keys=True))

# now go through the Transcript list
transcripts = obj['Transcript']
chromosome = "chr" + str(obj['seq_region_name'])
exonTree = IntervalTree()
for trs in transcripts:
    # go through each transcript, and store coordinate intervals
    for exon in trs['Exon']:
        start = exon['start']
        end = exon['end']
        exonTree[start:end] = (start,end)

exonTree.merge_overlaps()
for exon in sorted(exonTree.items()):
    print(chromosome + "\t" + str(exon.begin) + "\t" + str(exon.end))
#        print("chr" + str(exon['seq_region_name']) + "\t" + str(exon['start']) + "\t" + str(exon['end']))
