configfile: "config.yaml"

wildcard_constraints:
    scenario="[a-zA-Z0-9\.\-\_]+",
    country="[a-zA-Z]+"

rule solve_all:
    input:
        expand("summaries/" + config['run'] + "/{country}-{scenario}.csv",
	        **config['run_settings'])
    output: "summaries/" + config['run'] + "/all.csv"
    threads: 2
    resources:
        mem_mb=2000
    script: "summarise_networks.py"


rule summarise_network:
    input: "networks/" + config['run'] + "/{country}-{scenario}.nc"
    output: "summaries/" + config['run'] + "/{country}-{scenario}.csv"
    threads: 4
    resources:
        mem_mb=2000
    script: "summarise_network.py"

rule solve:
    output: "networks/" + config['run'] + "/{country}-{scenario}.nc"
    threads: 4
    resources:
        mem_mb=20000
    script: "solve.py"
