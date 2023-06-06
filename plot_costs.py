#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import matplotlib.pyplot as plt
import os
import yaml


# In[ ]:


# Plotting order (from bottom to top)
preferred_order = pd.Index([
    "solar",
    "wind",
    "battery",
    "hydrogen electrolyser",
    "hydrogen storage",
    "hydrogen turbine",
    "methanol storage",
    "methanol synthesis",
    "air separation unit",
    "oxygen storage",
    #"liquid oxygen storage",
    "Allam cycle",
    "heat pump",
    "direct air capture",
    "co2 storage",
    "CCGT",
    "methanol source",
    ])

# with open("../../efuels-server/config.yaml", "r") as f:
#     config = yaml.safe_load(f)
# colors = config["colors"]
# colors["OCGT"] = "orange"
# colors["CCGT"] = "brown"
# colors["CCGT+CC"] = "#5C4033"
# colors["biogenic co2"] = "green"
# colors["Allam cycle"] = "green"
# colors["heat pump"] = "b"
# colors["oxygen storage"] = "k"
# colors["air separation unit"] = "r"
# colors["dac"] = "#40e0d0"
# colors["direct air capture"] = colors["dac"]
# colors["hydrogen electrolyser"] = "#ffcdff"
# colors["hydrogen turbine"] = "#90EE90"

# manually setup colors
# colors = {
#     "CCGT": "#000000",
#     "Allam cycle": "#000000",
#     "hydrogen turbine": "#000000",
#     "oxygen storage": "#b4baea",
#     "air separation unit": "#586da6",
#     "heat pump": "#b4baea",
#     "co2 storage": "#586da6",
#     "direct air capture": "#8591c8",
#     "methanol storage": "#cfb93d",
#     "methanol synthesis": "#948b63",
#     "hydrogen electrolyser": "#586da6",
#     "hydrogen storage": "#8591c8",
#     "battery": "#cfb93d",
#     "solar": "#fae159",
#     "wind": "#084a96",
# }

# Automatically setup colors using seaborn colorblind-friendly colormap
import seaborn as sns
sns.set_style("ticks")
from itertools import cycle

colormap = sns.color_palette("colorblind", as_cmap=True)
colormap = colormap[-2:] + colormap[:-2] # start with yellow for solar

# map colors to labels
colors = {}
for label, color in zip(preferred_order, cycle(colormap)):
    colors[label] = color


# In[ ]:


if "snakemake" not in globals():
    # For runs outside snakemake, simple mock_snakemake
    from types import SimpleNamespace
    folder = "summaries/230601-71a-liquidco2allamo2fix/"

    member = {
        "input": {"statistics": folder+"statistics.csv"},
        "output": {"costs": folder+"costs.pdf"}
    }
    snakemake = SimpleNamespace(**member)


# In[ ]:


df = pd.read_csv(snakemake.input["statistics"], index_col=0)


# In[ ]:


df.loc["status"]


# In[ ]:


df = df.drop("status").astype(float)


# In[ ]:


df.loc[df.index[df.index.str.contains("totex")]].sum()/df.loc["total_load"] - df.loc["mean price electricity"]


# In[ ]:


fig, ax = plt.subplots()
fig.set_size_inches((14,4))


def rename(name):
    if "battery" in name:
        return "battery"
    elif "hydrogen stor" in name:
        return "hydrogen storage"
    elif name == "dac":
        return "direct air capture"
    elif name in ["oxygen liquefaction", "oxygen evaporation", "oxygen storage"]:
        return "oxygen storage"
    elif name in ["co2 liquefaction", "co2 evaporation", "co2 storage"]:
        return "co2 storage"
    else:
        return name

def rename_col(name):
    name = name.replace("-3a","").replace("-71a","").replace("-10a","").replace("H2s-wm-nH2t","MeOH").replace("-1H","").replace("mflex50-ramp5","lowflex").replace("mflex0-ramp10","highflex")
    name = name.replace("highflex-ccgt","CCGT")
    name = name.replace("-","\n")
    return name

costs = df.loc[df.index[df.index.str.contains("totex")]].multiply(1/df.loc["total_load"],axis=1)



costs.rename(lambda x: x[:-6],
             inplace=True)

costs.rename(rename_col,
             axis=1,
             inplace=True)

costs = costs.drop(costs.index.intersection(["co2 vent", "load", "oxygen vent", "oxygen storage standing losses"]))

rename_s = pd.Series(index=costs.index,
                     data=[rename(i) for i in costs.index])


costs = costs.groupby(rename_s).sum()


new_index = preferred_order.intersection(costs.index).append(costs.index.difference(preferred_order))
costs = costs.loc[new_index]

#print(costs)

costs.T.plot(kind="bar",stacked=True,color=[colors[i] for i in costs.index],
             linewidth=0,
             ax=ax,
             rot=0)

ax.set_ylim([0,190])

handles,labels = ax.get_legend_handles_labels()

handles.reverse()
labels.reverse()


ax.set_ylabel("average system electricity cost [€/MWh]")

#https://stackoverflow.com/questions/4700614/how-to-put-the-legend-outside-the-plot

# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
legend = ax.legend(handles, labels, loc="center left", bbox_to_anchor=(1, 0.5))


# annotate the top containers with the cumulative sum
#ax.bar_label(ax.containers[-1],padding=3,
#             fmt="%.1f")

y_offset = 4
for i, total in enumerate(costs.sum()):
    ax.text(i, total + y_offset, round(total), ha='center')
          #weight='bold')

#fig.tight_layout()

fig.savefig(snakemake.output["costs"],
            transparent=True,
            bbox_extra_artists=(legend,),
            bbox_inches='tight')


# %%
