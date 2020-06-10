#from intervaltree import Interval, IntervalTree
#t = IntervalTree()
#t.addi(0,13,1)
#t.addi(14,26,2)
#t.addi(27,38,0.8)
## show all values
#print(sorted(t))
## query
#
#print( sorted(t[10])[0].data )
#
#segment = IntervalTree(sorted(t[10:30]))
#data = []
#for interval_obj in segment:
#    data.append(interval_obj.data)
#
#print(max(data))


from intervaltree import Interval, IntervalTree
import csv

t = IntervalTree()
with open('chr3.bed', 'r') as f:
    reader = csv.reader(f, dialect='excel', delimiter=' ')
    for row in reader:
        start = int(row[1])
        end = int(row[2])
        data = float(row[3]) 
        if start == end:
            end = start + 1
        t.addi( start, end, data)

line = ""
for i in range(0,198295559,1000000):    # size of chr3
    value = 0       # default value is zero
    if len(t[i]) != 0:
        data = []
        for interval_obj in t[i]:
            data.append(interval_obj.data)
        value = max(data)
    line = line + str(value) + ","
line = line + "0"
print(line)
