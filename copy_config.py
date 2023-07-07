
from shutil import copy

files = [
    "config.yaml",
    "Snakefile",
    "solve.py",
    "defaults.csv",
]

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake('copy_config')

    for f in files:
        copy(f,snakemake.output[0].replace("config.yaml",f))
