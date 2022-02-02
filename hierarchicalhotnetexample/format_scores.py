import pandas as pd
import numpy as np 
import sys

input_file, output_file = sys.argv[1:]

df = pd.read_csv(input_file,sep="\t")

df = df.drop_duplicates("gene")
df['score'] = df['pval'].apply(lambda x: -np.log10(x))
df.loc[df[df['pval'] > 0.1].index,"score"] = 0

df[['gene','score']].to_csv(output_file,index=False,header=None,sep="\t")
