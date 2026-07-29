
# Filesystem operations
import shutil

# Extract zip files
import zipfile

configfile: "config.yaml"

wildcard_constraints:
    scenario="[a-zA-Z0-9\.\-\_]+",
    country="[a-zA-Z]+"

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
    group: "results"
    threads: 1
    resources:
        mem_mb=2000
    script: "summarise_networks.py"

rule summarise_network:
    input: "networks/" + config['run'] + "/{country}-{scenario}.nc"
    output: "summaries/" + config['run'] + "/statistics-{country}-{scenario}.csv"
    group: "scenario"
    resources:
        mem_mb=50000
    script: "summarise_network.py"

rule balances_networks:
    input:
        expand("summaries/" + config['run'] + "/balances-{country}-{scenario}.csv",
            **config['run_settings'])
    output: "summaries/" + config['run'] + "/balances.csv"
    group: "results"
    threads: 1
    resources:
        mem_mb=2000
    script: "summarise_networks.py"

rule balances_network:
    input: "networks/" + config['run'] + "/{country}-{scenario}.nc"
    output: "summaries/" + config['run'] + "/balances-{country}-{scenario}.csv"
    group: "scenario"
    resources:
        mem_mb=50000
    script: "balances_network.py"

rule get_fillings:
    input:
        expand("filling/" + config['run'] + "/filling-{country}-{scenario}.csv",
            **config['run_settings'])

rule get_filling:
    input: "networks/" + config['run'] + "/{country}-{scenario}.nc"
    output: "filling/" + config['run'] + "/filling-{country}-{scenario}.csv"
    group: "scenario"
    resources:
        mem_mb=50000
    script: "get_filling.py"

rule solve:
    input:
        solar="data/NUTS_0_sp_historical.nc",
        onwind0="data/NUTS_0_wp_ons_sim_0_historical_loc_weighted.nc",
        onwind1="data/NUTS_0_wp_ons_sim_1_historical_loc_weighted.nc",
        temperature="data/NUTS_0_t2m_detrended_timeseries_historical_pop_weighted.nc"
    output:
        "networks/" + config['run'] + "/{country}-{scenario}.nc"
    group: "scenario"
    threads: 4
    resources:
        mem_mb=50000
    benchmark:
        "benchmarks/" + config['run'] + "/{country}-{scenario}.tsv"
    log:
        solver="logs/" + config['run'] + "/{country}-{scenario}-solver.log",
        python="logs/" + config['run'] + "/{country}-{scenario}-python.log",
        memory="logs/" + config['run'] + "/{country}-{scenario}-memory.log",
    script: "solve.py"

rule plot_cost:
    input:
        statistics="summaries/" + config['run'] + "/statistics.csv",
    output:
        costs="summaries/" + config['run'] + "/costs.pdf",
    group: "results"
    script:
        "plot_costs.py"

rule plot_cascade:
    input:
        statistics="summaries/" + config['run'] + "/statistics.csv",
    output:
        scenarios=[
        'summaries/'+config['run']+'/cascade-DE-71a-1H-H2u.pdf',
        'summaries/'+config['run']+'/cascade-DE-71a-1H-H2s.pdf',
        'summaries/'+config['run']+'/cascade-DE-71a-1H-H2s-wm-nH2t-mflex0-ramp10.pdf',
        'summaries/'+config['run']+'/cascade-DE-71a-1H-H2s-wm-nH2t-mflex50-ramp5.pdf',
        'summaries/'+config['run']+'/cascade-DE-71a-1H-H2s-wm-nH2t-mflex0-ramp10-ccgt.pdf']
    group: "results"
    script:
        "plot_cascade.py"

rule plot_storage_filling_level:
    input:
        methanol="networks/" + config['run'] + "/DE-71a-1H-H2s-wm-nH2t-mflex0-ramp10.nc",
        hydrogen="networks/" + config['run'] + "/DE-71a-1H-H2u.nc",
    output:
        figure="summaries/" + config['run'] + "/storage_filling_level.pdf"
    group: "results"
    script:
        "plot_storage_filling_level.py"

rule plot_all:
    input:
        rules.plot_cost.output,
        rules.plot_cascade.output,
        rules.plot_storage_filling_level.output,
    group: "results"

rule all:
    input:
        rules.solve_all.input,
        rules.plot_all.input,
#    default_target: True
#    localrule: True
