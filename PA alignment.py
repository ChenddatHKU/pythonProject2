from Bio import AlignIO
import os
import glob
from Bio import SeqIO

def statistics (infile):
    n = 0
    list = []
    list_nonL = []
    align = AlignIO.read(infile, 'fasta')
    for record in align:
        n += 1
        if len(record) < 716 or record[0] != 'M': continue
        if record.seq[106] != 'X' and '-':
            list.append(record.seq[106])
            if record.seq[106] != 'L':
                list_nonL.append(record.seq[106])
                print(record, record.seq)
    # print(len(list), n)
    year = infile.rsplit('/')[6].rsplit('.')[0]
    totalsequence = n
    nonLlist = list_nonL
    nonLPercentage  = len(list_nonL) / len(list) * 100

    return [year, totalsequence, nonLlist, nonLPercentage]

    # print("total sequence# = " + str(n))
    # print("nonLlist =" + str(list_nonL))
    # print("non L%" + " in" + ' ' + str(n) + "=" + str(len(list_nonL)/len(list) *100))

outfile = open('/Users/chendd/Documents/HKU_Experiment/B Yamagata PA alignment/result', 'w')
header = "\t".join(['year', 'totalsequence', 'nonLlist', 'nonLPercentage'])
outfile.write(header + "\n")
filenames = sorted(glob.glob('/Users/chendd/Documents/HKU_Experiment/B Yamagata PA alignment/*.fas')) #change the location of NGS raw data accordingly
for filename in filenames:
  print(filename)
  fileID = filename.rsplit('/')[6].rsplit('.')[0]       # change the method to extract fileID according to the file name of your own.
  result = statistics(filename)
  print(result)
  outfile.write(str(result[0]) + "\t" + str(result[1]) + "\t" + str(result[2]) + "\t" + str(result[3]) + "\n")


# ##    if sth wrong, check the raw fas file, mainly check the location of L. sometimes "-" will be added when do multiple alignment, thus the location of L maybe changed accordingly depending on how many "-" added.
# align = AlignIO.read('/Users/chendd/Documents/HKU_Experiment/PA alignment/190120-200720PA.fas', 'fasta')
# n = 0
# list = []
# list_nonL = []
# for record in align:
#     print(record.seq, record.seq[106]) # change the number accordingly.
#     n += 1
#     if len(record) < 716 or record[0] != 'M': continue
#     if record.seq[106] != 'X' and '-':
#         list.append(record.seq[106])
#         if record.seq[106] != 'L':
#             list_nonL.append(record.seq[106])
#             print(record, record.seq)
#
# totalsequence = n
# nonLlist = list_nonL
# nonLPercentage  = len(list_nonL) / len(list) * 100
# print(totalsequence, nonLlist, nonLPercentage)

