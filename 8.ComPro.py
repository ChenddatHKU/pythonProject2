import glob

outfile = open('/Users/chendd/LowFreMut/screening/diversity', 'w')
header = "\t".join(['amplicon', 'wt', 'single', 'double', 'multiple'])
outfile.write(header + "\n")

summaryComposition = {}
amps = ['flu3amp1', 'flu3amp2', 'flu3amp3', 'flu3amp4', 'flu3amp5', 'flu3amp6', 'flu3amp7', 'flu3amp8', 'flu3amp9']
mutsC = ['wt', 'single', 'double', 'multiple', 'depth']
for amp in amps:
    summaryComposition[amp] = {}
    for muts in mutsC:
        summaryComposition[amp][muts] = 0

input = open('/Users/chendd/LowFreMut/BAT1-OP/lib1-9-wt.m', 'r')
for line in input.readlines():
    amp = line.rstrip().rsplit('\t')[0]
    summaryComposition[amp]['depth'] += 1
    if line.rstrip().rsplit('\t')[2] == "WT":
        summaryComposition[amp]['wt'] += 1
    elif line.rstrip().rsplit('\t')[2] != "WT" and line.rstrip().rsplit('\t')[2].count('-') == 0:
        summaryComposition[amp]['single'] += 1
    elif line.rstrip().rsplit('\t')[2] != "WT" and line.rstrip().rsplit('\t')[2].count('-') == 1:
        summaryComposition[amp]['double'] += 1
    elif line.rstrip().rsplit('\t')[2] != "WT" and line.rstrip().rsplit('\t')[2].count('-') > 1:
        summaryComposition[amp]['multiple'] += 1
print(summaryComposition)

sumP = {}
for amp in amps:
    sumP[amp] = {}
    for mut in mutsC:
        sumP[amp][mut] = 0
        sumP[amp][mut] = summaryComposition[amp][mut] / summaryComposition[amp]['depth']
    out = "\t".join(
        [amp, str(sumP[amp]['wt']), str(sumP[amp]['single']), str(sumP[amp]['double']), str(sumP[amp]['multiple']), ])
    outfile.write(out + '\n')

print(sumP)
##line 10: change input file names.