

import pandas as pd

df = pd.DataFrame()

def add(df,s,name):
    df = df.reindex(df.index.union(s.index))
    df[name] = s
    return df


if "balances" in snakemake.output[0]:
    index_col = (0,1)
else:
    index_col = 0

for fn in snakemake.input:
    print(fn)
    s = pd.read_csv(fn,
                    index_col=index_col,
                    header=None).squeeze("columns")

    fn = fn[fn.rfind("/")+1:-4]
    name = fn[fn.find("-")+1:]
    df = add(df,
             s,
             name)

df.to_csv(snakemake.output[0])
