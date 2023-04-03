

import pandas as pd, pypsa

n = pypsa.Network(snakemake.input[0])

s = pd.Series()

s.loc["total_load"] = n.loads_t.p.multiply(n.snapshot_weightings["generators"],axis=0)["load"].sum()
s.loc["objective"] = n.objective

s.loc["lcoe"] = s.loc["objective"]/s.loc["total_load"]


for name,attr in [("generators","p"),("links","p"),("stores","e")]:
    s = pd.concat((s,
                   getattr(n,name)[attr + "_nom_opt"]))

s = pd.concat((s,
               n.statistics().loc["Generator"]["Curtailment"].rename(lambda x: f"curtailment {x}")))

print(s)

s.to_csv(snakemake.output[0],
         header=False)
