# Import package
import pandas as pd

# import data
fcst_wa = pd.read_csv("tmax_05.csv", index_col = 0)

# print table
print(fcst_wa.loc[["Tahoua"],["tmxdy3"]])
print(fcst_wa)
