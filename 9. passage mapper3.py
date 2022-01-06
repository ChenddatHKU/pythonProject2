#!/usr/bin/python
import os
import sys
import glob
# from string import atof   #delete duo to update of python version
def floor(n):
  if float(n) < 0.0001: return 0.0001   #replace atof by float   #太不要脸了，<=0.001, regard as 0.001
  else: return n


infile  = open('/Users/chendd/LowFreMut/passage/changed pos','r')
outfile = open('/Users/chendd/LowFreMut/passage/HCMut','w')
silfile = open('/Users/chendd/LowFreMut/passage/sil','w')
misfile = open('/Users/chendd/LowFreMut/passage/mis','w')
nonfile = open('/Users/chendd/LowFreMut/passage/non','w')
header  = "\t".join(['Pos','Mut','AA','P4+BXA-1-abFre','P4+BXA-2-abFre','P4-1-abFre','P4-2-abFre','DS-1', 'DS-2', 'avDS'])
outfile.write(header+"\n")
countline  = 0
mincov     = 15000
count      = 0   #what's its function?
counthc    = 0
for line in infile.readlines():   #replace xreadlines by readlines
  if countline == 0: countline += 1; continue   #skip header
  if countline == 1: countline +=1; continue   #I skip two lines duo to file(changed pos), when I creat it, I didnt't add \n after header.
  line       = line.rstrip().rsplit("\t")
  Pos        = int(line[0])
  G          = line[1]
  F1         = line[2]
  F2         = line[3]

  withBXA_P4_1         = float(line[4])
  withBXA_P4_1_WT      = float(line[5])
  withBXA_P4_1_Dep     = float(line[6])

  withBXA_P4_2      = float(line[7])
  withBXA_P4_2_WT   = float(line[8])
  withBXA_P4_2_Dep  = float(line[9])

  passage_p4_1      = float(line[10])
  passage_p4_1_WT   = float(line[11])
  passage_p4_1_Dep  = float(line[12])

  passage_p4_2      = float(line[13])
  passage_p4_2_WT   = float(line[14])
  passage_p4_2_Dep  = float(line[15])

  withBXA_P4_1_fre = withBXA_P4_1/ withBXA_P4_1_Dep
  withBXA_P4_2_fre = withBXA_P4_2/ withBXA_P4_2_Dep
  passage_p4_1_fre = passage_p4_1/ passage_p4_1_Dep
  passage_p4_2_fre = passage_p4_2/ passage_p4_2_Dep

  count += 1
  if withBXA_P4_1_fre < 0.0004 and withBXA_P4_2_fre < 0.0004: continue   #minimum frequency?
  counthc += 1
  DS_1 = withBXA_P4_1_fre/floor(passage_p4_1_fre)
  DS_2 = withBXA_P4_2_fre/floor(passage_p4_2_fre)
  avDS = (DS_1+DS_2) /2

  out = "\t".join([str(Pos),G,F1, str(withBXA_P4_1_fre), str(withBXA_P4_2_fre), str(passage_p4_1_fre), str(passage_p4_2_fre), str(DS_1), str(DS_2), str(avDS)])
  outfile.write(out+"\n")
  if Pos > 24 and Pos < 750:
    if F1[-1] == '_': nonfile.write(out+"\n")
    elif F1[0] != F1[-1]: misfile.write(out+"\n")
    elif F1[0] == F1[-1]: silfile.write(out+"\n")
infile.close()
silfile.close()
misfile.close()
nonfile.close()
print(count, counthc)
#line 69, why DNA_RF    = (DNA_A+DNA_B)/(DNA_A_WT+DNA_B_WT), while there is a problem. Determine which one is better.