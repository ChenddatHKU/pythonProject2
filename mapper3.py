#!/usr/bin/python
import os
import sys
import glob
# from string import atof   #delete duo to update of python version

def floor(n):
  if float(n) < 0.001: return str(0.001)   #replace atof by float   #太不要脸了，<=0.001, regard as 0.001
  else: return n

infile  = open('/Users/tmp_flu3MC/changed pos','r')
outfile = open('/Users/tmp_flu3MC/HCMut','w')
silfile = open('/Users/tmp_flu3MC/sil','w')
misfile = open('/Users/tmp_flu3MC/mis','w')
nonfile = open('/Users/tmp_flu3MC/non','w')
header  = "\t".join(['Pos','Mut','AA','DNA_A','DNA_B','Tra_A','Tra_B','Inf_A1','Inf_A2','Fit'])
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
  WT         = float(line[4])
  WT_WT      = float(line[5])
  WT_Dep     = float(line[6])
  DNA_A      = float(line[7])
  DNA_A_WT   = float(line[8])
  DNA_A_Dep  = float(line[9])
  DNA_B      = float(line[10])
  DNA_B_WT   = float(line[11])
  DNA_B_Dep  = float(line[12])
  Tra_A      = float(line[13])
  Tra_A_WT   = float(line[14])
  Tra_A_Dep  = float(line[15])
  Tra_B      = float(line[16])
  Tra_B_WT   = float(line[17])
  Tra_B_Dep  = float(line[18])
  Inf_A1     = float(line[19])
  Inf_A1_WT  = float(line[20])
  Inf_A1_Dep = float(line[21])
  Inf_A2     = float(line[22])
  Inf_A2_WT  = float(line[23])
  Inf_A2_Dep = float(line[24])
  if DNA_A_Dep < mincov or DNA_B_Dep < mincov or Tra_A_Dep < mincov or Tra_B_Dep < mincov or Inf_A1_Dep < mincov or Inf_A2_Dep < mincov: continue   #quality control, filter out mutations depth < 15000
  WTfreq  = WT/WT_Dep   #absolute frequency
  DNAfreq = (DNA_A+DNA_B)/(DNA_A_Dep+DNA_B_Dep)   #average absolute frequency
  Trafreq = Tra_A/Tra_A_Dep   #absolute frequency
  count += 1
  if DNAfreq < WTfreq*4: continue   #filter out low confidence mutations
  if DNAfreq < 5/float(mincov): continue   #minimum frequency?
  counthc += 1
  DNA_A_RF  = str(DNA_A/DNA_A_WT)
  DNA_B_RF  = str(DNA_B/DNA_B_WT)
  Tra_A_RF  = str(Tra_A/Tra_A_WT)
  Tra_B_RF  = str(Tra_B/Tra_B_WT)
  Inf_A1_RF = str(Inf_A1/Inf_A1_WT)
  Inf_A2_RF = str(Inf_A2/Inf_A2_WT)

  DNA_RF    = (DNA_A+DNA_B)/(DNA_A_WT+DNA_B_WT)   #average fre of DNA
  Tra_RF    = Tra_A/Tra_A_WT   #average fre of Tra
  Inf_RF    = (Inf_A1+Inf_A2)/(Inf_A1_WT+Inf_A2_WT)   #average fre of inf
  fit = str((float(Inf_A1_RF)+float(Inf_A2_RF))/DNA_RF/2)   # fitness = average fre of infection / fre of DNA
  out = "\t".join([str(Pos),G,F1,DNA_A_RF,DNA_B_RF,Tra_A_RF,Tra_B_RF,Inf_A1_RF,Inf_A2_RF,fit])
  outfile.write(out+"\n")
  if Pos > 400 and Pos < 1800:
    if F1[-1] == '_': nonfile.write(str(floor(fit))+"\n")
    elif F1[0] != F1[-1]: misfile.write(str(floor(fit))+"\n")
    elif F1[0] == F1[-1]: silfile.write(str(floor(fit))+"\n")
infile.close()
silfile.close()
misfile.close()
nonfile.close()
print(count, counthc)
#line 69, why DNA_RF    = (DNA_A+DNA_B)/(DNA_A_WT+DNA_B_WT), while there is a problem. Determine which one is better.