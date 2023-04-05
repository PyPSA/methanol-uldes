

import pandas as pd

df = pd.DataFrame()

def add(df,s,name):
    df = df.reindex(df.index.union(s.index))
    df[name] = s
    return df

for fn in snakemake.input:
    print(fn)
    s = pd.read_csv(fn,
                    index_col=0,
                    header=None).squeeze("columns")

    name = fn[fn.rfind("/")+1:-4]
    df = add(df,
             s,
             name)

df.to_csv(snakemake.output[0])
