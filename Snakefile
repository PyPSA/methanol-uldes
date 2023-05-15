# Download files using Snakemake
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

# Move files around using Python
from shutil import move

# Extract zip files
import zipfile

configfile: "config.yaml"

wildcard_constraints:
    scenario="[a-zA-Z0-9\.\-\_]+",
    country="[a-zA-Z]+"

# Download timeseries Bloomfield and Brayshaw (2021) fomr https://doi.org/10.17864/1947.000321
rule download_timeseries:
    input:
        HTTP.remote(
            "researchdata.reading.ac.uk/321/4/ERA5_data_1950-2020.zip",
        ),
    output:
        solar="data/NUTS_0_sp_historical.nc",
        onwind0="data/NUTS_0_wp_ons_sim_0_historical_loc_weighted.nc",
        onwind1="data/NUTS_0_wp_ons_sim_1_historical_loc_weighted.nc",
    run:
        # Files to extract
        fps = [
            "ERA5_data_1950-2020/solar_power_capacity_factor/NUTS_0_sp_historical.nc",
            "ERA5_data_1950-2020/wp_onshore/NUTS_0_wp_ons_sim_0_historical_loc_weighted.nc",
            "ERA5_data_1950-2020/wp_onshore/NUTS_0_wp_ons_sim_1_historical_loc_weighted.nc",
        ]
        # extract and move files to this dir
        output_dir = "data/"
        with zipfile.ZipFile(input[0], "r") as zf:
            for fp in fps:
                zf.extract(fp, output_dir)
                move(output_dir + fp, output_dir + str(Path(fp).name))

rule solve_all:
    input:
        statistics="summaries/" + config['run'] + "/statistics.csv",
        balances="summaries/" + config['run'] + "/balances.csv",
        config="summaries/" + config['run'] + "/config.yaml"

rule copy_config:
    output: "summaries/" + config['run'] + "/config.yaml"
    threads: 1
    resources: mem_mb=1000
    script: "copy_config.py"

rule summarise_networks:
    input:
        expand("summaries/" + config['run'] + "/statistics-{country}-{scenario}.csv",
            **config['run_settings'])
    output: "summaries/" + config['run'] + "/statistics.csv"
    threads: 1
    resources:
        mem_mb=2000
    script: "summarise_networks.py"


rule summarise_network:
    input: "networks/" + config['run'] + "/{country}-{scenario}.nc"
    output: "summaries/" + config['run'] + "/statistics-{country}-{scenario}.csv"
    threads: 1
    resources:
        mem_mb=2000
    script: "summarise_network.py"

rule balances_networks:
    input:
        expand("summaries/" + config['run'] + "/balances-{country}-{scenario}.csv",
            **config['run_settings'])
    output: "summaries/" + config['run'] + "/balances.csv"
    threads: 1
    resources:
        mem_mb=2000
    script: "summarise_networks.py"

rule balances_network:
    input: "networks/" + config['run'] + "/{country}-{scenario}.nc"
    output: "summaries/" + config['run'] + "/balances-{country}-{scenario}.csv"
    threads: 1
    resources:
        mem_mb=2000
    script: "balances_network.py"


rule solve:
    input:
        solar=rules.download_timeseries.output.solar,
        onwind0=rules.download_timeseries.output.onwind0,
        onwind1=rules.download_timeseries.output.onwind1,
    output:
        "networks/" + config['run'] + "/{country}-{scenario}.nc"
    threads: 4
    resources:
        mem_mb=50000,
        walltime="08:00:00"
    log:
        solver="logs/" + config['run'] + "/{country}-{scenario}-solver.log",
        python="logs/" + config['run'] + "/{country}-{scenario}-python.log",
        memory="logs/" + config['run'] + "/{country}-{scenario}-memory.log",
    script: "solve.py"
