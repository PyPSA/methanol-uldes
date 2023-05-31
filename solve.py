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




def run_optimisation(assumptions, pu):
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

    if assumptions["voll"]:
        network.add("Generator","load",
                    bus="electricity",
                    carrier="load",
                    marginal_cost=assumptions["voll_price"],
                    p_max_pu=0,
                    p_min_pu=-1,
                    p_nom=assumptions["load"])
    elif assumptions["elastic"]:
        network.add("Generator","load",
                    bus="electricity",
                    carrier="load",
                    marginal_cost=assumptions["elastic_intercept"],
                    p_max_pu=0,
                    p_min_pu=-1,
                    p_nom=assumptions["load"])
    else:
        network.add("Load","load",
                    bus="electricity",
                    carrier="load",
                    p_set=assumptions["load"])

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

    network.add("Bus",
                "liquid oxygen",
                carrier="liquid oxygen")

    network.add("Generator",
                "oxygen vent",
                bus="oxygen",
                p_max_pu=0,
                p_min_pu=-1,
                p_nom=1e6,
                carrier="oxygen vent")

    network.add("Link",
                "air separation unit",
                bus0="electricity",
                bus1="oxygen",
                carrier="air separation unit",
                p_nom_extendable=True,
                capital_cost=assumptions_df.at["air_separation_unit","fixed"]*assumptions["air_separation_unit_efficiency"],
                efficiency=assumptions["air_separation_unit_efficiency"])

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
                p_nom_extendable=True,
                capital_cost=1, #Prevent storage cycling
                marginal_cost=0.1, #Prevent storage cycling
    )

    network.add("Store",
                "oxygen storage",
                bus="liquid oxygen",
                carrier="liquid oxygen storage",
                e_nom_extendable=True,
                e_cyclic=True,
                capital_cost=assumptions_df.at["oxygen_storage","fixed"],
                #standing_loss=np.power(assumptions["oxygen_storage_standing_loss"]/100, 1/24), # standing losses in %/day, convert to p.u./hour
                )


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
                "hydrogen_compressor",
                carrier="hydrogen storing compressor",
                bus0="hydrogen",
                bus1="compressed hydrogen",
                bus2="electricity",
                p_nom_extendable=True,
                efficiency=1,
                efficiency2=-assumptions["hydrogen_compressor_electricity"],
                capital_cost=assumptions_df.at["hydrogen_compressor","fixed"])

    network.add("Link",
                "hydrogen_decompressor",
                carrier="hydrogen storing decompressor",
                bus0="compressed hydrogen",
                bus1="hydrogen",
                p_nom_extendable=True)

    network.add("Store",
                "hydrogen_energy",
                bus="compressed hydrogen",
                carrier="hydrogen storage",
                e_nom_extendable=True,
                e_cyclic=True,
                capital_cost=assumptions_df.at["hydrogen_energy","fixed"])


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
                    capital_cost=assumptions_df.at["dac","fixed"]/assumptions["dac_electricity"],
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


        network.add("Store",
                    "co2",
                    bus="co2",
                    carrier="co2 storage",
                    e_nom_extendable=True,
                    e_cyclic=True,
                    capital_cost=assumptions_df.at["co2_storage","fixed"])

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

    if assumptions["methanol"] and not assumptions["ccgt"]:
        network.add("Link",
                    "Allam",
                    bus0="methanol",
                    bus1="electricity",
                    bus2="co2",
                    bus3="oxygen",
                    carrier="Allam cycle",
                    p_nom_extendable=True,
                    efficiency=assumptions["allam_cycle_efficiency"]/100.,
                    efficiency2=(assumptions["allam_cycle_co2_capture_efficiency"]/100.)*assumptions["methanolisation_co2"],
                    efficiency3=(-1)*assumptions["allam_cycle_o2"],
                    capital_cost=assumptions_df.at["allam_cycle","fixed"]*(assumptions["allam_cycle_efficiency"]/100.))

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

    if assumptions["elastic"]:
        network.generators["quadratic_coefficient"] = 0.
        #create inverse demand curve where elastic_intercept is price p where demand d
        #vanishes and load is demand d for zero p
        #inverse demand curve: p(d) = intercept - intercept/load*d
        #utility: U(d) = intercept*d - intercept/(2*load)*d^2
        #since demand is negative generator, take care with signs!
        network.generators.at["load","quadratic_coefficient"] = assumptions["elastic_intercept"]/(2*assumptions["load"])
    
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
    if "nH2t" in opts:
        assumptions["hydrogen"] = False
    if "bioco2" in opts:
        assumptions["biogenic_co2"] = True
    if "ccgt" in opts:
        assumptions["ccgt"] = True
    if "ocgt" in opts:
        assumptions["ocgt"] = True
    if "H2s" in opts:
        assumptions["hydrogen_energy_cost"] = 13
        assumptions["hydrogen_energy_fom"] = 2
        assumptions["hydrogen_energy_lifetime"] = 20
    elif "H2u" in opts:
        pass
    else:
        print("no H2 storage defined")
        sys.exit()

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

    years = int(opts[0][:-1])
    print(years,"years to optimise")
    assumptions["year_start"] = 2020 - years + 1
    assumptions["year_end"] = 2020

    print("optimising from",assumptions["year_start"],"to",assumptions["year_end"])
    
    n, message = run_optimisation(assumptions,pu)
    n.status = message
    
    n.export_to_netcdf(snakemake.output[0])
