

import pandas as pd, pypsa, numpy as np, sys

n = pypsa.Network(snakemake.input[0])

s = pd.Series()


if "load" in n.loads.index:
    s.loc["total_load"] = n.loads_t.p.multiply(n.snapshot_weightings["generators"],axis=0)["load"].sum()
elif "load" in n.generators.index:
    s.loc["total_load"] = -n.generators_t.p.multiply(n.snapshot_weightings["generators"],axis=0)["load"].sum()
else:
    print("cannot find load")
    sys.exit()
s.loc["objective"] = n.objective

s.loc["lcoe"] = s.loc["objective"]/s.loc["total_load"]


stats = n.statistics(aggregate_time="sum").groupby(level=1).sum()

stats["Total Expenditure"] = stats[["Capital Expenditure","Operational Expenditure"]].sum(axis=1)

for name,full_name in [("capex","Capital Expenditure"),("opex","Operational Expenditure"),("totex","Total Expenditure"),("capacity","Optimal Capacity")]:
    s = pd.concat((s,
                   stats[full_name].rename(lambda x: f"{x} {name}")))

s = pd.concat((s,
               n.statistics().loc["Generator"]["Curtailment"].rename(lambda x: f"{x} curtailment")))

s = pd.concat((s,
               n.statistics().groupby(level=1).sum()["Capacity Factor"].rename(lambda x: f"{x} capacity factor")))

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


#shares of elec
ts = pd.concat([(-1)**i*n.links_t["p"+str(i)][n.links.index[n.links["bus"+str(i)].map(n.buses.carrier) == "electricity"]] for i in range(3)]
               + [n.generators_t.p[n.generators.index[n.generators.bus.map(n.buses.carrier) == "electricity"]]],
               axis=1)
ts["load"] = n.loads_t.p["load"]

s = pd.concat((s,
               (ts.sum()/ts["load"].sum()).rename(lambda x: x+ " share")))


if "air separation unit" in n.links.index:
    s.loc["air separation unit oxygen capacity"] = n.links.at["air separation unit","p_nom_opt"]*n.links.at["air separation unit","efficiency"]

if "dac" in n.links.index:
    s.loc["dac co2 capacity"] = n.links.at["dac","p_nom_opt"]*n.links.at["dac","efficiency"]

if "Allam" in n.links.index:
    s.loc["Allam cycle electricity capacity"] = n.links.at["Allam","p_nom_opt"]*n.links.at["Allam","efficiency"]
    s.loc["methanol demand weighted price"] =  (n.buses_t.marginal_price["methanol"]*n.links_t.p0["Allam"]).sum()/n.links_t.p0["Allam"].sum()
    s.loc["oxygen demand weighted price"] =  (n.buses_t.marginal_price["oxygen"]*n.links_t.p0["Allam"]).sum()/n.links_t.p0["Allam"].sum()

if "CCGT" in n.links.index:
    s.loc["CCGT electricity capacity"] = n.links.at["CCGT","p_nom_opt"]*n.links.at["CCGT","efficiency"]
    s.loc["methanol demand weighted price"] =  (n.buses_t.marginal_price["methanol"]*n.links_t.p0["CCGT"]).sum()/n.links_t.p0["CCGT"].sum()

if "methanol synthesis" in n.links.index:
    s.loc["co2 demand weighted price"] =  (n.buses_t.marginal_price["co2"]*n.links_t.p0["methanol synthesis"]).sum()/n.links_t.p0["methanol synthesis"].sum()
    s.loc["hydrogen demand weighted price"] =  (n.buses_t.marginal_price["hydrogen"]*n.links_t.p0["methanol synthesis"]).sum()/n.links_t.p0["methanol synthesis"].sum()

if "hydrogen_turbine" in n.links.index:
    s.loc["hydrogen demand weighted price"] =  (n.buses_t.marginal_price["hydrogen"]*n.links_t.p0["hydrogen_turbine"]).sum()/n.links_t.p0["hydrogen_turbine"].sum()


s.loc["status"] = n.status

print(s)

s.to_csv(snakemake.output[0],
         header=False)
