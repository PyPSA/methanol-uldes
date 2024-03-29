## Copyright 2018-2020 Tom Brown

## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU Affero General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU Affero General Public License for more details.

## License and more information at:
## https://github.com/PyPSA/whobs-server



import pypsa

import sys

import pandas as pd
from pyomo.environ import Constraint

import json, os, hashlib, yaml

import xarray as xr

import scipy as sp

import numpy as np

import logging
logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

defaults = pd.read_csv("defaults.csv",index_col=[0,1],na_filter=False)

for (n,t) in [("f",float),("i",int)]:
    defaults.loc[defaults["type"] == n, "value"] = defaults.loc[defaults["type"] == n,"value"].astype(t)
#work around fact bool("False") returns True
defaults.loc[defaults.type == "b","value"] = (defaults.loc[defaults.type == "b","value"] == "True")

defaults_t = {str(year): defaults.swaplevel().loc[str(year)] for year in config["tech_years"]}
defaults_nt = defaults.swaplevel().loc[""]

default_assumptions = pd.concat((defaults_nt,defaults_t[str(config["tech_years_default"])])).sort_index()

if "snakemake" in globals():
    datasets = {
        "onwind0": snakemake.input["onwind0"],
        "onwind1": snakemake.input["onwind1"],
        "solar":  snakemake.input["solar"],
    }
else:
    # same files, but snakemake is not loaded
    # from https://doi.org/10.17864/1947.000321
    datasets = {
        "onwind0": "data/NUTS_0_wp_ons_sim_0_historical_loc_weighted.nc",
        "onwind1": "data/NUTS_0_wp_ons_sim_1_historical_loc_weighted.nc",
        "solar": "data/NUTS_0_sp_historical.nc",
    }


df = {}

snapshots = pd.date_range("1950-01-01","2020-12-31 23:00",
                          freq="H")

for key,value in datasets.items():
    ds = xr.open_dataset(value)
    df[key] = pd.DataFrame(data=ds["timeseries_data"].T,
                           index=snapshots,
                           columns=ds["NUTS_keys"].values)

ds = xr.open_dataset(snakemake.input.temperature)
df["temperature"] = pd.DataFrame(data=ds["detrended_data"][:,4,:].T,
                                 index=snapshots,
                                 columns=ds["NUTS_keys"].values)


current_version = config["current_version"]

octant_folder = config["octant_folder"]

#based on mean deviation against renewables.ninja capacity factors for European countries for 2011-2013
solar_correction_factor = 0.926328


override_component_attrs = pypsa.descriptors.Dict(
        {k: v.copy() for k, v in pypsa.components.component_attrs.items()}
    )
override_component_attrs["Link"].loc["bus2"] = [
        "string",
        np.nan,
        np.nan,
        "2nd bus",
        "Input (optional)",
    ]
override_component_attrs["Link"].loc["bus3"] = [
        "string",
        np.nan,
        np.nan,
        "3rd bus",
        "Input (optional)",
    ]
override_component_attrs["Link"].loc["efficiency2"] = [
        "static or series",
        "per unit",
        1.0,
        "2nd bus efficiency",
        "Input (optional)",
    ]
override_component_attrs["Link"].loc["efficiency3"] = [
        "static or series",
        "per unit",
        1.0,
        "3rd bus efficiency",
        "Input (optional)",
    ]
override_component_attrs["Link"].loc["p2"] = [
        "series",
        "MW",
        0.0,
        "2nd bus output",
        "Output",
    ]
override_component_attrs["Link"].loc["p3"] = [
        "series",
        "MW",
        0.0,
        "3rd bus output",
        "Output",
    ]



def annuity(lifetime,rate):
    if rate == 0.:
        return 1/lifetime
    else:
        return rate/(1. - 1. / (1. + rate)**lifetime)


assumptions_df = pd.DataFrame(columns=["FOM","fixed","discount rate","lifetime","investment"],
                              dtype=float)

threshold = 0.1




def run_optimisation(assumptions, pu, scenario_opts):
    """Needs cleaned-up assumptions and pu.
    return results_overview, results_series, error_msg"""


    year_start = assumptions['year_start']
    year_end = assumptions['year_end']

    Nyears = year_end - year_start + 1

    print(Nyears,"years considered")

    techs = [tech[:-5] for tech in assumptions if tech[-5:] == "_cost" and tech[-14:] != "_marginal_cost" and tech != "co2_cost"]

    print("calculating costs for",techs)


    for item in techs:
        assumptions_df.at[item,"discount rate"] = assumptions[item + "_discount"]/100.
        assumptions_df.at[item,"investment"] = assumptions[item + "_cost"]*1e3 if "EUR/kW" in defaults.loc[item + "_cost"]["unit"][0] else assumptions[item + "_cost"]
        assumptions_df.at[item,"FOM"] = assumptions[item + "_fom"]
        assumptions_df.at[item,"lifetime"] = assumptions[item + "_lifetime"]

    assumptions_df["fixed"] = [(annuity(v["lifetime"],v["discount rate"])+v["FOM"]/100.)*v["investment"]*Nyears for i,v in assumptions_df.iterrows()]

    print('Starting task for {} with assumptions {}'.format(assumptions["location"],assumptions_df))

    network = pypsa.Network(override_component_attrs=override_component_attrs)

    snapshots = pd.date_range("{}-01-01".format(year_start),"{}-12-31 23:00".format(year_end),
                              freq=str(assumptions["frequency"])+"H")

    network.set_snapshots(snapshots)

    network.snapshot_weightings = pd.Series(float(assumptions["frequency"]),index=network.snapshots)

    network.add("Bus","electricity",
                carrier="electricity")

    load = assumptions["load"]
    if assumptions["temperature_demand"]:
        tdep_demand = config["degree_cutoff"] - df["temperature"][country]
        tdep_demand[tdep_demand < 0.] = 0.
        load += (assumptions["temperature_demand_mean"]/tdep_demand.mean())*tdep_demand


    if assumptions["voll"]:
        network.add("Generator","load",
                    bus="electricity",
                    carrier="load",
                    marginal_cost=assumptions["voll_price"],
                    p_max_pu=0,
                    p_min_pu=-1,
                    p_nom=assumptions["load"])
    elif assumptions["elastic"]:
        #create inverse demand curve where elastic_intercept is price p where demand d
        #vanishes and load is demand d for zero p
        #inverse demand curve: p(d) = intercept - intercept/load*d
        #utility: U(d) = intercept*d - intercept/(2*load)*d^2
        #since demand is negative generator, take care with signs!
        network.add("Generator","load",
                    bus="electricity",
                    carrier="load",
                    marginal_cost=assumptions["elastic_intercept"],
                    marginal_cost_quadratic=assumptions["elastic_intercept"]/(2*assumptions["load"]),
                    p_max_pu=0,
                    p_min_pu=-1,
                    p_nom=assumptions["load"])
    else:
        network.add("Load","load",
                    bus="electricity",
                    carrier="load",
                    p_set=load)

    if assumptions["solar"]:
        network.add("Generator","solar",
                    bus="electricity",
                    carrier="solar",
                    p_max_pu = pu["solar"],
                    p_nom_extendable = True,
                    p_nom_min = assumptions["solar_min"],
                    p_nom_max = assumptions["solar_max"],
                    marginal_cost = 0.1, #Small cost to prefer curtailment to destroying energy in storage, solar curtails before wind
                    capital_cost = assumptions_df.at['solar','fixed'])

    if assumptions["wind"]:
        network.add("Generator","wind",
                    bus="electricity",
                    carrier="wind",
                    p_max_pu = pu["onwind"],
                    p_nom_extendable = True,
                    p_nom_min = assumptions["wind_min"],
                    p_nom_max = assumptions["wind_max"],
                    marginal_cost = 0.2, #Small cost to prefer curtailment to destroying energy in storage, solar curtails before wind
                    capital_cost = assumptions_df.at['wind','fixed'])



    for i in range(1,3):
        name = "dispatchable" + str(i)
        if assumptions[name]:
            network.add("Carrier",name,
                        co2_emissions=assumptions[name+"_emissions"])
            network.add("Generator",name,
                        bus="electricity",
                        carrier=name,
                        p_nom_extendable=True,
                        marginal_cost=assumptions[name+"_marginal_cost"],
                        capital_cost=assumptions_df.at[name,'fixed'])

    if assumptions["battery"]:

        network.add("Bus","battery",
                    carrier="battery")

        network.add("Store","battery_energy",
                    bus = "battery",
                    carrier="battery storage",
                    e_nom_extendable = True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at['battery_energy','fixed'])

        network.add("Link","battery_power",
                    bus0 = "electricity",
                    bus1 = "battery",
                    carrier="battery inverter",
                    efficiency = assumptions["battery_power_efficiency_charging"]/100.,
                    p_nom_extendable = True,
                    capital_cost=assumptions_df.at['battery_power','fixed'])

        network.add("Link","battery_discharge",
                    bus0 = "battery",
                    bus1 = "electricity",
                    carrier="battery discharger",
                    p_nom_extendable = True,
                    efficiency = assumptions["battery_power_efficiency_discharging"]/100.)


    network.add("Bus",
                "hydrogen",
                carrier="hydrogen")

    network.add("Bus",
                "oxygen",
                carrier="oxygen")

    network.add("Generator",
                "oxygen vent",
                bus="oxygen",
                p_max_pu=0,
                p_min_pu=-1,
                p_nom=1e6,
                carrier="oxygen vent")

    if assumptions["hydrogen_load"] != 0:
        network.add("Load","hydrogen_load",
                    bus="hydrogen",
                    carrier="hydrogen",
                    p_set=assumptions["hydrogen_load"])

    network.add("Link",
                "hydrogen_electrolyser",
                bus0="electricity",
                bus1="hydrogen",
                bus2="oxygen",
                carrier="hydrogen electrolyser",
                p_nom_extendable=True,
                efficiency=assumptions["hydrogen_electrolyser_efficiency"]/100.,
                efficiency2=8*assumptions["hydrogen_electrolyser_efficiency"]/100./33,  # divide by 33 to get tH2, multiply by 8 to get tO2
                capital_cost=assumptions_df.at["hydrogen_electrolyser","fixed"])

    network.add("Bus",
                "compressed hydrogen",
                carrier="compressed hydrogen")

    network.add("Link",
                "hydrogen_decompressor",
                carrier="hydrogen storing decompressor",
                bus0="compressed hydrogen",
                bus1="hydrogen",
                p_nom=1e6)

    # TODO e.g. "location" -> use coco based on 2letter country code
    # Depending on storage technology for H2 use different assumptions for compression and storage
    if "H2s" in scenario_opts:
        # H2 storage steel tanks (type 1)
        network.add("Link",
                "hydrogen_compressor",
                carrier="hydrogen storing compressor",
                bus0="hydrogen",
                bus1="compressed hydrogen",
                bus2="electricity",
                p_nom_extendable=True,
                efficiency=1,
                efficiency2=-assumptions["hydrogen_storage_tank_compressor_electricity"],
                capital_cost=assumptions_df.at["hydrogen_storage_tank_compressor","fixed"])

        network.add("Store",
                    "hydrogen_energy",
                    bus="compressed hydrogen",
                    carrier="hydrogen storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["hydrogen_storage_tank","fixed"],
                    )

    elif "H2u" in scenario_opts:
        # H2 storage caverns underground
        network.add("Link",
                "hydrogen_compressor",
                carrier="hydrogen storing compressor",
                bus0="hydrogen",
                bus1="compressed hydrogen",
                bus2="electricity",
                p_nom_extendable=True,
                efficiency=1,
                efficiency2=-assumptions["hydrogen_storage_cavern_compressor_electricity"],
                capital_cost=assumptions_df.at["hydrogen_storage_cavern_compressor","fixed"])

        network.add("Store",
                    "hydrogen_energy",
                    bus="compressed hydrogen",
                    carrier="hydrogen storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["hydrogen_storage_cavern","fixed"],
                    )

    elif "H2l" in scenario_opts:
        # H2 storage steel tanks (type 1)
        network.add("Link",
                "hydrogen_compressor",
                carrier="hydrogen storing compressor",
                bus0="hydrogen",
                bus1="compressed hydrogen",
                bus2="electricity",
                p_nom_extendable=True,
                efficiency=1,
                efficiency2=-assumptions["hydrogen_storage_tank_compressor_electricity"],
                capital_cost=assumptions_df.at["hydrogen_storage_tank_compressor","fixed"])

        network.add("Store",
                    "hydrogen_energy",
                    bus="compressed hydrogen",
                    carrier="hydrogen storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["hydrogen_storage_tank","fixed"],
                    )

        network.add("Bus",
                    "liquid hydrogen",
                    carrier="liquid hydrogen")

        # LH2 tanks
        network.add("Link",
                    "hydrogen_liquefaction",
                    carrier="hydrogen liquefaction",
                    bus0="compressed hydrogen",
                    bus1="liquid hydrogen",
                    bus2="electricity",
                    p_nom_extendable=True,
                    efficiency=assumptions["hydrogen_liquefaction_efficiency"],
                    efficiency2=-assumptions["hydrogen_liquefaction_electricity"]*assumptions["hydrogen_liquefaction_efficiency"],
                    capital_cost=assumptions_df.at["hydrogen_liquefaction","fixed"])

        network.add("Store",
                    "liquid hydrogen",
                    bus="liquid hydrogen",
                    carrier="hydrogen storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["liquid_hydrogen_storage","fixed"],
                    )

        network.add("Link",
                    "hydrogen_deliquefier",
                    carrier="hydrogen storing deliquefier",
                    bus0="liquid hydrogen",
                    bus1="compressed hydrogen",
                    p_nom=1e6)

        # Later implemented via custom constraint
        network.add("Link",
                    "liquid hydrogen storage standing losses",
                    bus0="liquid hydrogen",
                    bus1="compressed hydrogen",
                    carrier="liquid hydrogen storage standing losses",
                    efficiency=1,
                    p_nom=1e6 #dummy value instead of capacity expansion
        )

    if assumptions["hydrogen"]:

        network.add("Link",
                     "hydrogen_turbine",
                     bus0="hydrogen",
                     bus1="electricity",
                     carrier="hydrogen turbine",
                     p_nom_extendable=True,
                     efficiency=assumptions["hydrogen_turbine_efficiency"]/100.,
                     capital_cost=assumptions_df.at["hydrogen_turbine","fixed"]*assumptions["hydrogen_turbine_efficiency"]/100.)  #NB: fixed cost is per MWel

    if assumptions["methanol"] or assumptions["methanol_load"] != 0:

        network.add("Bus",
                    "co2",
                    carrier="co2")

        network.add("Bus",
                    "heat",
                    carrier="heat")

        network.add("Link",
                    "heat pump",
                    bus0="electricity",
                    bus1="heat",
                    carrier="heat pump",
                    p_nom_extendable=True,
                    capital_cost=assumptions_df.at["heat_pump","fixed"]*assumptions["heat_pump_efficiency"]/100.,
                    efficiency=assumptions["heat_pump_efficiency"]/100.)

        network.add("Link",
                    "dac",
                    bus0="electricity",
                    bus1="co2",
                    bus2="heat",
                    carrier="dac",
                    p_nom_extendable=True,
                    capital_cost=assumptions["dac_factor"]*assumptions_df.at["dac","fixed"]/assumptions["dac_electricity"],
                    efficiency=1/assumptions["dac_electricity"],
                    efficiency2=-assumptions["dac_heat"]/assumptions["dac_electricity"])

        network.add("Generator",
                    "co2 vent",
                    bus="co2",
                    p_max_pu=0,
                    p_min_pu=-1,
                    p_nom=1e6,
                    carrier="co2 vent")

        if assumptions["biogenic_co2"]:
            network.add("Generator",
                        "biogenic co2",
                        bus="co2",
                        p_nom=1e6,
                        marginal_cost=assumptions["biogenic_co2_price"],
                        carrier="biogenic co2")

        # Intermediary bus for pressure conversion before CO2 liquefaction
        network.add("Bus",
                    "co2 high pressure",
                    carrier="co2 high pressure")

        # Intermediary link; assume difference in liquefaction between low and high pressure for efficiency and cost
        # Reason: DAC CO2 output is assumed to be at low pressure and requires pressure increase
        # Allam cycle CO2 output is assumed to be at high pressure and does not require pressure increase
        network.add("Link",
                    "co2 compression",
                    carrier="co2 compression",
                    bus0="co2",
                    bus1="co2 high pressure",
                    bus2="electricity",
                    efficiency=1,
                    efficiency2=-assumptions["co2_liquefaction_efficiency"]+assumptions["co2_liquefaction_high_pressure_efficiency"],
                    p_nom_extendable=True,
                    capital_cost=assumptions_df.at["co2_liquefaction","fixed"] - assumptions_df.at["co2_liquefaction_high_pressure","fixed"],
                    )

        network.add("Bus",
                    "co2 liquid storage",
                    carrier="co2 liquid storage",
                    )

        network.add("Link",
                    "co2 liquefaction",
                    bus0="co2 high pressure",
                    bus1="co2 liquid storage",
                    bus2="electricity",
                    carrier="co2 liquefaction",
                    p_nom_extendable=True,
                    efficiency=1,
                    efficiency2=-assumptions["co2_liquefaction_high_pressure_efficiency"],
                    capital_cost=assumptions_df.at["co2_liquefaction_high_pressure","fixed"],
                    )

        network.add("Store",
                    "co2 liquid storage",
                    bus="co2 liquid storage",
                    carrier="co2 liquid storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["co2_storage","fixed"],
        )

        network.add("Link",
                    "co2 evaporation",
                    bus0="co2 liquid storage",
                    bus1="co2",
                    carrier="co2 evaporation",
                    p_nom_extendable=False,
                    efficiency=1,
                    p_nom=1e6 #dummy value instead of capacity expansion
                    )

        network.add("Bus",
                     "methanol",
                     carrier="methanol")

        network.add("Store",
                    "methanol",
                    bus="methanol",
                    carrier="methanol storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["liquid_carbonaceous_storage","fixed"]/config["mwh_per_m3"]["methanol"])

        network.add("Link",
                    "methanol synthesis",
                    bus0="hydrogen",
                    bus1="methanol",
                    bus2="electricity",
                    bus3="co2",
                    carrier="methanol synthesis",
                    p_nom_extendable=True,
                    p_min_pu=assumptions["methanolisation_min_part_load"]/100,
                    ramp_limit_up=assumptions["ramp"],
                    ramp_limit_down=assumptions["ramp"],
                    efficiency=assumptions["methanolisation_efficiency"],
                    efficiency2=-assumptions["methanolisation_electricity"]*assumptions["methanolisation_efficiency"],
                    efficiency3=-assumptions["methanolisation_co2"]*assumptions["methanolisation_efficiency"],
                    capital_cost=assumptions_df.at["methanolisation","fixed"]*assumptions["methanolisation_efficiency"]) #NB: cost is EUR/kW_MeOH

        if assumptions["meohsource"]:
            network.add("Generator",
                        "methanol source",
                        bus="methanol",
                        p_nom=1e6,
                        marginal_cost=assumptions["meohsource_marginal_cost"],#60 mimics 50 EUR/MWh LNG + 50 EUR/tCO2 CCS OR v. cheap clean MeOH
                        carrier="methanol source")

    if assumptions["methanol"] and not assumptions["ccgt"]:
        network.add("Link",
                    "Allam",
                    bus0="methanol",
                    bus1="electricity",
                    bus2="co2 high pressure",
                    bus3="oxygen",
                    carrier="Allam cycle",
                    p_nom_extendable=True,
                    efficiency=assumptions["allam_cycle_efficiency"]/100.,
                    efficiency2=(assumptions["allam_cycle_co2_capture_efficiency"]/100.)*config["co2_per_mwh"]["methanol"],
                    efficiency3=(-1)*assumptions["allam_cycle_o2"],
                    capital_cost=assumptions["allam_factor"]*assumptions_df.at["allam_cycle","fixed"]*(assumptions["allam_cycle_efficiency"]/100.))

        #ASU only necessary with Allam to top-up oxygen from electrolyser
        network.add("Link",
                    "air separation unit",
                    bus0="electricity",
                    bus1="oxygen",
                    carrier="air separation unit",
                    p_nom_extendable=True,
                    capital_cost=assumptions_df.at["air_separation_unit","fixed"]*assumptions["air_separation_unit_efficiency"],
                    efficiency=assumptions["air_separation_unit_efficiency"])


        # Add oxygen storage in case of Allam cycle
        network.add("Bus",
                "liquid oxygen",
                carrier="oxygen storage")

        network.add("Link",
                    "oxygen liquefaction",
                    bus0="oxygen",
                    bus1="liquid oxygen",
                    bus2="electricity",
                    carrier="oxygen liquefaction",
                    efficiency=1, # Perfect liquefaction of O2
                    efficiency2=-assumptions["oxygen_storage_liquefaction_efficiency"],
                    p_nom_extendable=True,
                    capital_cost=1, #Prevent storage cycling
                    marginal_cost=0.1, #Prevent storage cycling
        )

        network.add("Link",
                    "oxygen evaporation",
                    bus0="liquid oxygen",
                    bus1="oxygen",
                    carrier="oxygen evaporation",
                    efficiency=1, # Perfect evaporation; don't assume additional energy for compression for Allam cycle feed
                    p_nom=1e6 #dummy value instead of capacity expansion
        )

        # Later implemented via custom constraint
        network.add("Link",
                    "oxygen storage standing losses",
                    bus0="liquid oxygen",
                    bus1="oxygen",
                    carrier="oxygen storage standing losses",
                    efficiency=1,
                    p_nom=1e6 #dummy value instead of capacity expansion
        )

        network.add("Store",
                    "oxygen storage",
                    bus="liquid oxygen",
                    carrier="oxygen storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["oxygen_storage","fixed"],
                    )

    if assumptions["ccgt"]:
        network.add("Link",
                    "CCGT",
                    bus0="methanol",
                    bus1="electricity",
                    carrier="CCGT",
                    p_nom_extendable=True,
                    efficiency=0.6,
                    capital_cost=assumptions_df.at["hydrogen_turbine","fixed"]*0.6)

    if assumptions["ocgt"]:
        network.add("Link",
                    "OCGT",
                    bus0="methanol",
                    bus1="electricity",
                    carrier="OCGT",
                    p_nom_extendable=True,
                    efficiency=0.35,
                    capital_cost=0.5*assumptions_df.at["hydrogen_turbine","fixed"]*0.35)

    if assumptions["methanol_load"] != 0:
        network.add("Load","methanol_load",
                    bus="methanol",
                    carrier="methanol",
                    p_set=assumptions["methanol_load"])

    if assumptions["reformer"]:
        network.add("Link",
                    "reformer",
                    bus0="methanol",
                    bus1="hydrogen",
                    bus2="co2",
                    carrier="reformer",
                    p_nom_extendable=True,
                    efficiency=assumptions["reformer_efficiency"]/100.,
                    efficiency2=assumptions["reformer_capture_rate"],
                    capital_cost=assumptions_df.at["reformer","fixed"])


    if assumptions["co2_limit"]:
        network.add("GlobalConstraint","co2_limit",
                    sense="<=",
                    constant=assumptions["co2_emissions"]*assumptions["load"]*network.snapshot_weightings.objective.sum())

    network.consistency_check()

    solver_name = config["solver"]["name"]
    solver_options = config["solver_options"][config["solver"]["options"]]
    solver_logfile = snakemake.log["solver"]

    network.optimize.create_model()

    if assumptions["battery"]:
        network.model.add_constraints(network.model["Link-p_nom"].loc["battery_power"]
                                      -network.links.loc["battery_discharge", "efficiency"]*
                                      network.model["Link-p_nom"].loc["battery_discharge"] == 0,
                                      name='charger_ratio')

    if "oxygen storage standing losses" in network.links.index:
        # standing losses in %/day, convert to p.u./hour and round to 6 decimals because higher accuracy irrelevant
        hourly_standing_losses=np.round(1-np.power(1-assumptions["oxygen_storage_standing_loss"]/100, 1/24), 6)

        # Standing losses equivalent to state of charge of previous timestep
        network.model.add_constraints(
            hourly_standing_losses
            * network.model["Store-e"].loc[np.concatenate((network.snapshots[-1:], network.snapshots[:-1])), "oxygen storage"]
            - network.model["Link-p"].loc[:, "oxygen storage standing losses"]
            == 0,
            name="liquid-oxygen_standing-losses",
        )

    if "liquid hydrogen storage standing losses" in network.links.index:
        # standing losses in %/day, convert to p.u./hour and round to 6 decimals because higher accuracy irrelevant
        hourly_standing_losses=np.round(1-np.power(1-assumptions["liquid_hydrogen_storage_standing_loss"]/100, 1/24), 6)

        # Standing losses equivalent to state of charge of previous timestep
        network.model.add_constraints(
            hourly_standing_losses
            * network.model["Store-e"].loc[np.concatenate((network.snapshots[-1:], network.snapshots[:-1])), "liquid hydrogen"]
            - network.model["Link-p"].loc[:, "liquid hydrogen storage standing losses"]
            == 0,
            name="liquid-hydrogen_standing-losses",
        )

    if "H2u" in scenario_opts:
        # For hydrogen caverns compression and evaporation, i.e. injection and withdrawal, are constrained by ratio to each other
        # important as long as evaporation is free of cost
        logger.info("Adding hydrogen cavern withdrawal/injection ratio constraint")
        network.model.add_constraints(
            network.model["Link-p_nom"].loc["hydrogen_compressor"]
            - network.model["Link-p_nom"].loc["hydrogen_decompressor"] / assumptions["hydrogen_storage_cavern_withdrawal_injection_ratio"]
            == 0,
            name="hydrogen_cavern_withdrawal_injection_ratio",
        )

    if "clip_p_max_pu" in config["solver_options"]:
        network.generators_t.p_max_pu.where(lambda x: x > config["solver_options"]["clip_p_max_pu"], other=0.0, inplace=True)

    status, termination_condition = network.optimize.solve_model(solver_name=solver_name,
                                                                 solver_options=solver_options,
                                                                 log_fn=solver_logfile)

    print(status,termination_condition)

    if termination_condition in ["infeasible","infeasible or unbounded"]:
        return None, "Problem was infeasible"
    elif termination_condition in ["numeric"]:
        return None, "Numerical trouble encountered, problem could be infeasible"
    elif status == "ok" and termination_condition == "optimal":
        return network, "OK"
    elif status == "ok" and termination_condition == "suboptimal":
        return network, "suboptimal"
    elif status == "warning" and termination_condition == "suboptimal":
        return network, "suboptimal"
    else:
        return None, "Job failed to optimise correctly"

if __name__ == "__main__":


    # Detect running outside of snakemake and mock up snakemake for testing
    if 'snakemake' not in globals():
        from pypsa.descriptors import Dict
        import yaml
        from types import SimpleNamespace

        snakemake = SimpleNamespace()

        with open('config.yaml') as f:
            snakemake.config = yaml.safe_load(f)

        snakemake.wildcards = Dict({"country" : "DE",
                                    # "scenario" : "71a-1H-H2s",
                                    # "scenario" : "71a-1H-H2s-wm-nH2t-mflex0-ramp10",
                                    "scenario" : "1a-1H-H2s-wm-nH2t-mflex0-ramp10",
                                    })

        snakemake.output = ["results/{}-{}.nc".format(snakemake.wildcards.country,
                                                         snakemake.wildcards.scenario)]
        snakemake.log = {
            "python": f"logs/{snakemake.wildcards['country']}-{snakemake.wildcards['scenario']}-python.log",
            "solver": f"logs/{snakemake.wildcards['country']}-{snakemake.wildcards['scenario']}-solver.log",
            }

    country = snakemake.wildcards.country
    scenario = snakemake.wildcards.scenario

    logging.basicConfig(filename=snakemake.log["python"],
                        level=snakemake.config['logging_level'])


    logger.info(f"computing country {country} and scenario {scenario}")

    pu = pd.DataFrame()
    pu["onwind"] = df["onwind0"][country]
    pu["solar"] = df["solar"][country]


    assumptions = default_assumptions["value"].to_dict()
    assumptions["ramp"] = np.nan
    assumptions["meohsource"] = False
    assumptions["temperature_demand"] = False
    assumptions["dac_factor"] = 1.
    assumptions["allam_factor"] = 1.

    opts = scenario.split("-")
    if "wm" in opts:
        assumptions["methanol"] = True

        if not assumptions["ccgt"]:
            # Default is methanol + Allam cycle
            assumptions["air_separation_unit"] = True
            assumptions["oxygen_storage"] = True
            assumptions["allam_cycle_turbine"] = True
    if "wref" in opts:
        assumptions["reformer"] = True
    if "nowind" in opts:
        assumptions["wind"] = False
    if "nH2t" in opts:
        assumptions["hydrogen"] = False
    if "bioco2" in opts:
        assumptions["biogenic_co2"] = True
    if "ccgt" in opts:
        assumptions["ccgt"] = True
    if "ocgt" in opts:
        assumptions["ocgt"] = True
    if "pessimisticH2l" in opts:
        assumptions["hydrogen_liquefaction_electricity"] = 0.24
        assumptions["liquid_hydrogen_storage_standing_loss"] = 3

    for opt in opts:
        if opt[-1:] == "H":
            assumptions["frequency"] = int(opt[:-1])
        if opt[:2] == "ed":
            assumptions["load"] = float(opt[2:])
        if opt[:2] == "hd":
            assumptions["hydrogen_load"] = float(opt[2:])
        if opt[:2] == "md":
            assumptions["methanol_load"] = float(opt[2:])
        if opt[:5] == "mflex":
            assumptions["methanolisation_min_part_load"] = float(opt[5:])
            print("Methanol min part load set to",assumptions["methanolisation_min_part_load"])
        if opt[:4] == "voll":
            assumptions["voll"] = True
            assumptions["voll_price"] = float(opt[4:])
        if opt[:7] == "elastic":
            assumptions["elastic"] = True
            assumptions["elastic_intercept"] = float(opt[7:])
        if opt[:4] == "ramp":
            assumptions["ramp"] = float(opt[4:])/100
        if opt[:10] == "meohsource":
            assumptions["meohsource"] = True
            assumptions["meohsource_marginal_cost"] = float(opt[10:])
        if opt[:7] == "tdemand":
            assumptions["temperature_demand"] = True
            assumptions["temperature_demand_mean"] = float(opt[7:])
        if opt[:3] == "dac":
            assumptions["dac_factor"] = float(opt[3:])
        if opt[:5] == "allam":
            assumptions["allam_factor"] = float(opt[5:].replace("p","."))

    years = int(opts[0][:-1])
    print(years,"years to optimise")
    assumptions["year_start"] = 2020 - years + 1
    assumptions["year_end"] = 2020

    print("optimising from",assumptions["year_start"],"to",assumptions["year_end"])

    n, message = run_optimisation(assumptions,pu,opts)
    n.status = message

    n.export_to_netcdf(snakemake.output[0],
                       # compression of network
                       float32=True, compression={'zlib': True, "complevel":9, "least_significant_digit":5}
                       )
