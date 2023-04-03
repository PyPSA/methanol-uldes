configfile: "config.yaml"

wildcard_constraints:
    scenario="[a-zA-Z0-9\.\-\_]+",
    country="[a-zA-Z]+"

rule solve_all:
    input:
        expand("networks/{country}-{scenario}.nc",
	        **config['run_settings'])

rule solve:
    output: "networks/{country}-{scenario}.nc"
    threads: 4
    resources:
        mem_mb=20000
    script: "solve.py"
