#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pypsa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import seaborn as sns

sns.set_style("ticks")


# In[ ]:


if "snakemake" not in globals():
    # For runs outside snakemake, simple mock_snakemake
    from types import SimpleNamespace
    folder = "summaries/230509-71a-fixasu/"

    member = {
        "input": {
            "methanol": "networks/230509-71a-fixasu/DE-71a-1H-H2s-wm-nH2t-mflex0-ramp10.nc",
            "hydrogen": "networks/230509-71a-fixasu/DE-71a-1H-H2u.nc",
            },
        "output": {"figure": "filling_level.pdf"}
    }
    snakemake = SimpleNamespace(**member)


# In[ ]:


n = {}
n["m"] = pypsa.Network(snakemake.input["methanol"])
n["h"] = pypsa.Network(snakemake.input["hydrogen"])


# In[ ]:


fig, ax = plt.subplots()
fig.set_size_inches((6,4))


(n["h"].stores_t.e["hydrogen_energy"]/100/24).plot(ax=ax,label="underground hydrogen")
(n["m"].stores_t.e["methanol"]/100/24).plot(ax=ax,label="methanol")

ax.set_ylabel("storage filling level [days of demand]")
ax.set_ylim([0,np.ceil(ax.get_ylim()[1]/20)*20])
ax.set_xlabel("")
ax.legend(loc="upper right")


fig.savefig(snakemake.output["figure"],
            transparent=True,
            bbox_inches='tight')


# In[ ]:


n["m"].stores_t.e["oxygen storage"].plot()


# In[ ]:


n["m"].stores_t.e["co2"].plot()


# In[ ]:


s = n["m"].stores_t.e["methanol"]

fft = abs(pd.Series(np.fft.rfft(s - s.mean()),   # subtract mean!
                    index=np.fft.rfftfreq(len(s), d=1/len(s)))**2)


# In[ ]:


fft.plot()
#plt.xlim(0, 768)
plt.xlabel('1 / 71a')
plt.ylim(0,1e20)
plt.xlim(0,100)

