# !/usr/bin/python
import os
import sys
import glob
# from string import atof   #delete duo to update of python version

def floor(n):
  if float(n) < 0.001:
    return str(0.001)  # replace atof by float   #太不要脸了，<=0.001, regard as 0.001
  else:
    return n

infile = open('/Users/chendd/LowFreMut/screening/changed pos', 'r')
outfile = open('/Users/chendd/LowFreMut/screening/Normalized to Dep/HCMut', 'w')
silfile = open('/Users/chendd/LowFreMut/screening/Normalized to Dep/sil', 'w')
misfile = open('/Users/chendd/LowFreMut/screening/Normalized to Dep/mis', 'w')
nonfile = open('/Users/chendd/LowFreMut/screening/Normalized to Dep/non', 'w')
outfile2 = open('/Users/chendd/LowFreMut/screening/Normalized to Dep/HCMu2(mis)', 'w')
header = "\t".join(['Pos', 'Mut', 'AA', 'HBXA-1-RF', 'HBXA-2-RF', 'LBXA-1-RF', 'LBXA-2-RF', 'passage-1-RF','passage-2-RF', 'HBXA-C-1-RF', 'HBXA-C-2-RF', 'LBXA-C-1-RF', 'LBXA-C-2-RF','passage-C-1-RF','passage-C-2-RF','Fit', "FitC", "DS-HBXA", "DS-LBXA", "DS-HBXA-C", "DS-LBXA-C", "passage_1_fit", "passage-2-fit"])
outfile.write(header + "\n")
silfile.write(header + "\n")
misfile.write(header + "\n")
nonfile.write(header + "\n")
outfile2.write(header + "\n")


countline = 0
mincov = 15000
count = 0  # what's its function?
counthc = 0
for line in infile.readlines():  # replace xreadlines by readlines
  if countline == 0: countline += 1; continue  # skip header
  if countline == 1: countline +=1; continue
  line = line.rstrip().rsplit("\t")

  Pos = int(line[0])
  G = line[1]
  F1 = line[2]
  # F2         = line[3]
  DNA = float(line[4])
  DNA_WT = float(line[5])
  DNA_Dep = float(line[6])

  WT = float(line[7])
  WT_WT = float(line[8])
  WT_Dep = float(line[9])

  HBXA_1 = float(line[10])
  HBXA_1_WT = float(line[11])
  HBXA_1_Dep = float(line[12])

  LBXA_1 = float(line[13])
  LBXA_1_WT = float(line[14])
  LBXA_1_Dep = float(line[15])

  passage_1 = float(line[16])
  passage_1_WT = float(line[17])
  passage_1_Dep = float(line[18])

  HBXA_C_1 = float(line[19])
  HBXA_C_1_WT = float(line[20])
  HBXA_C_1_Dep = float(line[21])

  LBXA_C_1 = float(line[22])
  LBXA_C_1_WT = float(line[23])
  LBXA_C_1_Dep = float(line[24])

  passage_C_1 = float(line[25])
  passage_C_1_WT = float(line[26])
  passage_C_1_Dep = float(line[27])

  HBXA_2 = float(line[28])
  HBXA_2_WT = float(line[29])
  HBXA_2_Dep = float(line[30])

  LBXA_2 = float(line[31])
  LBXA_2_WT = float(line[32])
  LBXA_2_Dep = float(line[33])

  passage_2 = float(line[34])
  passage_2_WT = float(line[35])
  passage_2_Dep = float(line[36])

  HBXA_C_2 = float(line[37])
  HBXA_C_2_WT = float(line[38])
  HBXA_C_2_Dep = float(line[39])

  LBXA_C_2 = float(line[40])
  LBXA_C_2_WT = float(line[41])
  LBXA_C_2_Dep = float(line[42])

  passage_C_2 = float(line[43])
  passage_C_2_WT = float(line[44])
  passage_C_2_Dep = float(line[45])

  WTfreq = WT / WT_Dep  # absolute frequency
  mutfreq = DNA / DNA_Dep  # absolute frequency
  count += 1
  if mutfreq < WTfreq * 4: continue  # filter out low confidence mutations
  if mutfreq < 6 / float(mincov): continue  # minimum reads count
  counthc += 1

  if Pos > 24 and Pos < 745:
    DNA_RF = str(DNA / DNA_Dep)   ## mut DAN baseline relative frequency to wt
    HBXA_1_RF = str(HBXA_1 / HBXA_1_Dep)
    HBXA_2_RF = str(HBXA_2 / HBXA_2_Dep)

    LBXA_1_RF = str(LBXA_1 / LBXA_1_Dep)
    LBXA_2_RF = str(LBXA_2 / LBXA_2_Dep)

    passage_1_RF = str(passage_1 / passage_1_Dep)
    passage_2_RF = str(passage_2 / passage_2_Dep)


    HBXA_C_1_RF = str(HBXA_C_1 / HBXA_C_1_Dep)
    HBXA_C_2_RF = str(HBXA_C_2 / HBXA_C_2_Dep)

    LBXA_C_1_RF = str(LBXA_C_1 / LBXA_C_1_Dep)
    LBXA_C_2_RF = str(LBXA_C_2 / LBXA_C_2_Dep)

    passage_C_1_RF = str(passage_C_1 / passage_C_1_Dep)
    passage_C_2_RF = str(passage_C_2 / passage_C_2_Dep)

# Relative fitness score calculation
    RF_1 = str(float(passage_1_RF)  / float(DNA_RF))
    RF_2 = str(float(passage_2_RF)  / float(DNA_RF))
    # RF_C_1 = str(float(passage_C_1_RF)  / float(DNA_RF))
    # RF_C_2 = str(float(passage_C_2_RF)  / float(DNA_RF))

    fit         = str( (float(passage_1_RF) + float(passage_2_RF)) / float(DNA_RF) / 2)  # fitness = average fre of infection / fre of DNA
    fitC        = str( (float(passage_C_1_RF) + float(passage_C_2_RF)) / float(DNA_RF) / 2)

# differential selection calculation
    passage_RF = str((float(passage_1_RF) + float(passage_2_RF)) / 2)
    passage_C_RF = str((float(passage_C_1_RF) + float(passage_C_2_RF)) / 2)



    dRF_HBXA    = str( (float(HBXA_1_RF) + float(HBXA_2_RF))/ float(passage_RF)/2 )
    dRF_LBXA    = str( (float(LBXA_1_RF) + float(LBXA_2_RF))/ float(passage_RF)/2 )
    dRF_HBXA_C  = str((float(HBXA_C_1_RF) + float(HBXA_C_2_RF))/ float(passage_C_RF))
    dRF_LBXA_C  = str((float(LBXA_C_1_RF) + float(LBXA_C_2_RF))/ float(passage_C_RF))








    out = "\t".join([str(Pos), G, F1, HBXA_1_RF, HBXA_2_RF, LBXA_1_RF, LBXA_2_RF, passage_1_RF, passage_2_RF, HBXA_C_1_RF, HBXA_C_2_RF, LBXA_C_1_RF, LBXA_C_2_RF, passage_C_1_RF, passage_C_2_RF,fit, fitC, dRF_HBXA, dRF_LBXA, dRF_HBXA_C, dRF_LBXA_C, RF_1, RF_2])
    outfile.write(out + "\n")
    if F1[-1] != F1[0] and F1[-1] != '_':
      outfile2.write(out + "\n")

  if Pos > 24 and Pos < 745:
    if F1[-1] == '_':
      nonfile.write(out + "\n")
    elif F1[0] != F1[-1]:
      misfile.write(out + "\n")
    elif F1[0] == F1[-1]:
      silfile.write(out + "\n")

infile.close()
silfile.close()
misfile.close()
nonfile.close()
print(count, counthc)
# line 69, why DNA_RF    = (DNA_A+DNA_B)/(DNA_A_WT+DNA_B_WT), while there is a problem. Determine which one is better.