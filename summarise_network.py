

import pandas as pd, pypsa, numpy as np

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

s = pd.concat((s,
               n.buses_t.marginal_price.mean().rename(lambda x: f"mean price {x}")))

s = pd.concat((s,
               n.generators_t.p.mean().rename(lambda x: f"mean generation {x}")))

s = pd.concat((s,
               n.generators_t.p_max_pu.multiply(n.generators.p_nom_opt).mean().rename(lambda x: f"mean available generation {x}")))

#market values
bus_map = (n.buses.carrier == "electricity")
bus_map.at[""] = False
for c in n.iterate_components(n.one_port_components):
    items = c.df.index[c.df.bus.map(bus_map).fillna(False)]
    if len(items) == 0:
        continue
    mv = (c.pnl.p[items].multiply(n.buses_t.marginal_price["electricity"], axis=0).sum()/c.pnl.p[items].sum()).groupby(c.df.loc[items,'carrier']).mean()
    s = pd.concat((s,
                   mv.rename(lambda x: x+ " mv").replace([np.inf, -np.inf], np.nan).dropna()))

for c in n.iterate_components(n.branch_components):
    for end in [col[3:] for col in c.df.columns if col[:3] == "bus"]:
        items = c.df.index[c.df["bus" + str(end)].map(bus_map,na_action=None)]
        if len(items) == 0:
            continue
        mv = (c.pnl["p"+end][items].multiply(n.buses_t.marginal_price["electricity"], axis=0).sum()/c.pnl["p"+end][items].sum()).groupby(c.df.loc[items,'carrier']).mean()
        s = pd.concat((s,
                       mv.rename(lambda x: x+ " mv").replace([np.inf, -np.inf], np.nan).dropna()))


s.loc["status"] = n.status

print(s)

s.to_csv(snakemake.output[0],
         header=False)
