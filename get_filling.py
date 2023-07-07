
import pypsa, pandas as pd

n = pypsa.Network(snakemake.input[0])
n.stores_t.e.to_csv(snakemake.output[0])
