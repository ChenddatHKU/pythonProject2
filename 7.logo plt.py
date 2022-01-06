# useful imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
# %matplotlib inline
import logomaker as lm
df1 = pd.read_excel(r"/Users/chendd/Desktop/AAarray.xlsx", sheet_name='Sheet2', index_col=0)
print(type((df1['R'][8])))
df2 = df1.replace(np.nan, 0.00001)    #replace NaN with 0
df3 = df2.loc[0:100]    # select specific rows in the dataframe, do not select all the data, otherwise the fig is too small
print(df3)
df4 = df2.loc[101:200]
logo1 = lm.Logo(df3,color_scheme='weblogo_protein')
logo2 = lm.Logo(df4,color_scheme='weblogo_protein')

print(logo1, logo2)

# add formatting statements
logo1.ax.set_ylabel('RF Score')
logo1.ax.set_xlabel('Site')

logo2.ax.set_ylabel('RF Score')
logo2.ax.set_xlabel('Site')

plt.show()