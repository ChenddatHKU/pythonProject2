import os
import sys
import glob


def readframe(aa, pos):
  wtaa  = aa[0]
  mutaa = aa[1]
  pos   = (pos-11)/3+1 # change with different gene  11 = len(nucleotide number before ATG) + 1
  return wtaa+str(pos)+mutaa


def adjustmutpos(muts, offset):
  if muts == 'WT': return muts   ### what does this mean? why return muts?
  else:
    muts  = muts.rsplit('-') #type(muts) =  list?
    mlist = []
    for m in muts:
      pos = str(int(m[1:-1])+offset-1)
      mlist.append(m[0]+pos+m[-1])
    return '-'.join(mlist)  #try this function it seems I can join each content in the lits vie a "-"


def list2string(l):
  newl = []
  for i in l:
    newl.append(str(i))
  return newl

# read in offset values for each amplicon
infile  = open('/Users/chendd/Documents/HKU_Experiment/ErrorCorrect/Source code/PAlib240-master/Fasta/flu3offset','r')
offsetH = {}
for line in infile.readlines():  #change xreadline() to readline() duo to upgrade of python version.
  line = line.rstrip().rsplit("\t")
  offsetH[line[0]] = [int(line[1]),int(line[2])]
infile.close()

#READ IN BARCODE FILE
infile   = open('/Users/chendd/Documents/HKU_Experiment/ErrorCorrect/Source code/PAlib240-master/Fasta/BarCode','r')
barcodes = {}
pops     = []
for line in infile.readlines():  #change xreadline() to readline() duo to upgrade of python version.
  line = line.rstrip().rsplit("\t")
  barcodes[line[0]] = line[1]
  pops.append(line[1])
infile.close()

#GENOTYPE AND DEPTH HASH INITIATE
GenotypeH = {}
depthH    = {}
WTcount   = {}
#assign key values to GenotypeH, depthH and WTcount
for pop in pops:
  GenotypeH[pop] = {}   #why do not induce dic in dic of GenotypeH? while dpethH and WTcount was added with subdic?
  depthH[pop]    = {}
  WTcount[pop]   = {}
  for amp in offsetH.keys():
    depthH[pop][amp] = 0
    WTcount[pop][amp] = 0

#############################MAIN##################################
infile = open('/Users/tmp_flu3MC/AllM','r')
for line in infile.readlines():
  line = line.rstrip().rsplit("\t")
  if line[1] not in barcodes.keys(): continue
  amp    = line[0]
  offset = offsetH[amp][0]
  pop    = barcodes[line[1]]
  muts   = adjustmutpos(line[2],offset)
 # if muts == 'G689T' : print line
  depthH[pop][amp] += 1
  if muts == 'WT':
    WTcount[pop][amp] += 1
  else:

    if muts in GenotypeH[pop].keys():   #change according to the update of python version
      GenotypeH[pop][muts] += 1
    else:
      GenotypeH[pop][muts] = 1   # these four lines assign key/value pair for GenotypeH[pop]
infile.close()

#SUMMARIZE GENOTYPE
Genotypes = []
for pop in pops:
  Genotypes.extend(GenotypeH[pop].keys())
Genotypes = list(set(Genotypes))  # set remove duplicates, then the out most list turn the set into a list again.
#Genotypes is a full list of all mutations, without duplicates.

Genotypes_out = open('/Users/tmp_flu3MC/Genotypes','w')
for G in Genotypes:
  Genotypes_out.write(G + "\n")
Genotypes_out.close()

#OUTPUT
#WTCOUNT AND DEPTH OUTPUT
outfile = open('/Users/tmp_flu3MC/depth','w')
outfile.write('pop'+"\t"+'amp'+"\t"+'WTcount'+"\t"+'Depth'+"\n")
for pop in pops:
  for amp in offsetH.keys():
    outfile.write("\t".join([pop, amp, str(WTcount[pop][amp]), str(depthH[pop][amp])])+"\n")
outfile.close()

#GENOTYPES
outfileA = open('/Users/tmp_flu3MC/AGenotypes','w')
outfileS = open('/Users/tmp_flu3MC/SGenotypes','w')
header   = 'Genotype'+"\t"+"\t".join(pops)
outfileA.write(header+"\n")
outfileS.write(header+"\n")

for G in Genotypes:
  counts = []
  for pop in pops:
    if G in GenotypeH[pop].keys():
      counts.append(GenotypeH[pop][G])
    else:
      counts.append(0)
  out = G+"\t"+"\t".join(list2string(counts))
  outfileA.write(out+"\n")
  if '-' not in G:   # this line remove double/multiple mutations, thus Sgenotypes stores only single mutants.
    outfileS.write(out+"\n")
outfileA.close()
outfileS.close()

#SINGLE MUTS INFO FILE
outfileC = open('/Users/tmp_flu3MC/flu3MC_T4','w')   #change the location of outfile, but do not change the file name.
header   = ['Pos','Genotype','Frame1','Frame2']
for pop in pops:
  header.extend([pop, pop+'_WT',pop+'_Dep'])
outfileC.write("\t".join(header)+"\n")

#READ IN PROTEIN LEVEL DATA
infile = open('/Users/chendd/Documents/HKU_Experiment/ErrorCorrect/Source code/PAlib240-master/Fasta/flu3info','r')   #read in array file. How to creat the array file?
Rframe = {}
for line in infile.readlines():
  array = line.rstrip().rsplit("\t")
  ID = ''.join(array[0:3])
  WTaa1 = array[3]
  Mutaa1 = array[4]
  WTaa2 = array[5]
  Mutaa2 = array[6]
  Rframe[ID] = [WTaa1+Mutaa1, WTaa2+Mutaa2]
infile.close()

# single mutation
for G in Genotypes:
  if '-' not in G and 'N' not in G:
    pos = int(G[1:-1])
#    if not Rframe.has_key(G): print G; continue
    F1  = readframe(Rframe[G][0],pos)
    F2  = readframe(Rframe[G][1],8)
    out = [str(pos),G,F1,F2]
    for pop in pops:
      wtc = 0
      dep = 0
      for amp in offsetH.keys():
        if pos >= offsetH[amp][0] and pos <= offsetH[amp][1]:
          wtc += WTcount[pop][amp]
          dep += depthH[pop][amp]
      if G in GenotypeH[pop].keys(): out.append(str(GenotypeH[pop][G]))
      else: out.append(str(0))
      out.append(str(wtc))
      out.append(str(dep))
    if len(out) < 5: continue
    else:
      outfileC.write("\t".join(out)+"\n")
outfileC.close()
#READ IN PROTEIN LEVEL DATA
infile = open('Fasta/flu3info','r')
Rframe = {}
for line in infile.readlines():
  array = line.rstrip().rsplit("\t")
  ID = ''.join(array[0:3])
  WTaa1 = array[3]
  Mutaa1 = array[4]
  WTaa2 = array[5]
  Mutaa2 = array[6]
  Rframe[ID] = [WTaa1+Mutaa1, WTaa2+Mutaa2]
infile.close()
for G in Genotypes:
  if '-' not in G and 'N' not in G:
    pos = int(G[1:-1])
    F1  = readframe(Rframe[G][0],pos)
    F2  = readframe(Rframe[G][1],8)   #why there is a F2, what's the biological meaning of F2?
    out = [str(pos),G,F1,F2]
    for pop in pops:
      wtc = 0
      dep = 0
      for amp in offsetH.keys():
        if pos >= offsetH[amp][0] and pos <= offsetH[amp][1]:   #append wtc and dep of amp.
          wtc += WTcount[pop][amp]
          dep += depthH[pop][amp]
      if G in GenotypeH[pop].keys(): out.append(str(GenotypeH[pop][G]))   # append reads count of G
      else: out.append(str(0))
      out.append(str(wtc))
      out.append(str(dep))
    outfileC.write("\t".join(out)+"\n")
outfileC.close()

