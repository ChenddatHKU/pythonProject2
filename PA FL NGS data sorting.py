import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile

## add phenotype info by merge AllMut file with flu3info file.
# wb1 = pd.read_excel(r"/Users/chendd/Desktop/PA FL serial passage analysis.xlsx", sheet_name='Sheet1')
# wb2 = pd.read_excel(r"/Users/chendd/Desktop/Book1.xlsx", sheet_name="Sheet2")
# print(wb1)
# print(wb2)
# outmerge = pd.merge(wb1, wb2, on="genotype", how="outer")
# print(type(outmerge))
# outmerge.to_excel(r"/Users/chendd/Desktop/Book2.xlsx", sheet_name="Sheet1 ")

fileID = []
genotype = []
df1 = pd.read_excel(r"/Users/chendd/Desktop/Book2.xlsx", sheet_name="Sheet1 ")
for ind in df1.index:
    if df1['depth'][ind] >= 1000:
        fileID.append(df1['ID'][ind])
        genotype.append(df1['genotype'][ind])
fileID= list(set(fileID))
genotypes = list(set(genotype))
L= len(genotypes)

#create output dataframe: df2
df2 = pd.DataFrame(columns=fileID, index=range(0,L))
df2['Nan'] = genotypes
#read in datafarame 1 data row by row
for ind in df1.index:
    genotype = df1['genotype'][ind]
    fre = df1['fre'][ind]
    fileid = df1['ID'][ind]
    pos = df1['pos'][ind]

    phenotype = df1['phenotype'][ind]
    position = df2.index[df2['Nan'] == genotype]
    df2.at[position, fileid] = fre
    df2.at[position, 'phenotype'] = phenotype
    df2.at[position, 'pos'] =pos
df2.to_excel(r"/Users/chendd/Desktop/Book3.xlsx", sheet_name="Sheet1 ")