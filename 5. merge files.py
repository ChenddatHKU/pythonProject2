import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile

wb1 = pd.read_excel(r"data loaction.xlsx", sheet_name='Sheet1')
wb2 = pd.read_excel(r"", sheet_name="Sheet2")
print(wb1)
print(wb2)
outmerge = pd.merge(wb1, wb2, on="merged name", how="outer")
print(outmerge)
outmerge.to_excel(r"data location", sheet_name="string ")

