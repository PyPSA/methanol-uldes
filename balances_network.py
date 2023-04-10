
import pypsa, pandas as pd

def calculate_balances(n):

    balances = pd.Series()

    bus_carriers = n.buses.carrier.unique()

    for i in bus_carriers:
        bus_map = (n.buses.carrier == i)
        bus_map.at[""] = False

        for c in n.iterate_components(n.one_port_components):

            items = c.df.index[c.df.bus.map(bus_map).fillna(False)]
            if len(items) == 0:
                continue

            s = c.pnl.p[items].multiply(n.snapshot_weightings.generators,axis=0).sum().multiply(c.df.loc[items, 'sign']).groupby(c.df.loc[items, 'carrier']).sum()
            s = pd.concat([s], keys=[i])
            balances = balances.reindex(s.index.union(balances.index))
            balances.loc[s.index] = s

        for c in n.iterate_components(n.branch_components):

            for end in [col[3:] for col in c.df.columns if col[:3] == "bus"]:

                items = c.df.index[c.df["bus" + str(end)].map(bus_map, na_action=False)]

                if len(items) == 0:
                    continue

                s = (-1)*c.pnl["p"+end][items].multiply(n.snapshot_weightings.generators,axis=0).sum().groupby(c.df.loc[items, 'carrier']).sum()
                s.index = s.index + end
                s = pd.concat([s], keys=[i])
                balances = balances.reindex(s.index.union(balances.index))
                balances.loc[s.index] = s
    return balances


if __name__ == "__main__":

    n = pypsa.Network(snakemake.input[0])
    balances = calculate_balances(n)
    balances.to_csv(snakemake.output[0],
                    header=False)
