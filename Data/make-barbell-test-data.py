import pandas as pd
import numpy as np
import os

script_dir = os.path.dirname(__file__) 
example_data_file = os.path.join(script_dir, "example_data.csv")

# make mock data for two sample comparison
df1 = pd.read_csv(example_data_file)
df2 = df1.copy()
df2['Entities pValue'] = df2['Entities pValue']  * (2 * np.random.rand(len(df2['Entities pValue'])))
df2['Entities FDR'] = df2['Entities FDR']  * (2 * np.random.rand(len(df2['Entities FDR'])))

# Subset modified sample data so only some overlapping pathways exists
df2 = df2[100:800]
df2.to_csv(os.path.join(script_dir, "mod_example_data.csv"), index=False)
