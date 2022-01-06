import pandas as pd

# this script is used to transform the HCMut file into an array, followed by plotting the array into sequence logos

# summary of mutated positions
infile = open('/Users/tmp_flu3MC/HCMut','r')
header = infile.readline()
poss = []
AA_poss = []
for line in infile.readlines():
    pos = line.rstrip().rsplit('\t')[0]
    AA_pos = line.rstrip().rsplit('\t')[2][1:-1]
    poss.append(pos)
    AA_poss.append(AA_pos)
infile.close()

Newposs = list(set(poss))
NewAA_poss = list(set(AA_poss))
L = len(Newposs)
L2 = len(NewAA_poss)
print("nucleotide position mutated = " + str(L), "\n" + "amina acid position mutated = " + str(L2))

df1 = pd.DataFrame(columns=['wt NT', 'pos', 'A', 'T', 'G', 'C'], index=range(0, L))
df1['pos'] = Newposs
df2 = pd.DataFrame(columns=['wt AA', "AApos", 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X'], index=range(0, L2))
df2['AApos'] = NewAA_poss
infile2 = open('/Users/tmp_flu3MC/HCMut','r')
header = infile2.readline()
for line in infile2.readlines():
    posNT = line.rstrip().rsplit('\t')[0]
    wtNT = line.rstrip().rsplit('\t')[1][0]
    MTNT = line.rstrip().rsplit('\t')[1][-1]
    fit = line.rstrip().rsplit('\t')[15] ##modify the number

    posAA = line.rstrip().rsplit('\t')[2][1:-1]
    wtAA = line.rstrip().rsplit('\t')[2][0]
    MutAA = line.rstrip().rsplit('\t')[2][-1]

    po = df1.index[df1['pos'] == posNT][0]    #get row number by value in column
    df1.at[po, MTNT.upper()] = fit
    df1.at[po, "wt NT"] = wtNT

    po2 = df2.index[df2['AApos'] == posAA][0]
    df2.at[po2, 'wt AA'] = wtAA
    # if wtAA == MutAA: continue
    # if wtAA != MutAA and MutAA != '-':
    df2.at[po2, MutAA.upper()] = fit


infile2.close()
df1.to_excel('/Users/tmp_flu3MC/NTarray.xlsx',sheet_name='Sheet1')
df2.to_excel('/Users/tmp_flu3MC/AAarray.xlsx',sheet_name='Sheet1')

# if you want to list missense muts, you can add line 45 and 46 to filter out missense muts, and line47 should be indented.
#check list
#line 6, 27, 51, 52, input, output file location
#line33, modify the fit column number






