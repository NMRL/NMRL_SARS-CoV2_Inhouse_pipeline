import sys, os
from pathlib import Path
sys.path.insert(0, f'{os.path.dirname(str(Path(__file__).absolute()))}/subscripts/')
from modules import run_all, run_assembly, run_downstream
from utilities import install_snakemake, parse_arguments, read_json_dict

"""
This is a wrapper script of COVIPIPE pipeline.
Date: 2022-06-13
Version: 0.0
"""

if __name__ == "__main__":
    args = parse_arguments(read_json_dict('./config_files/json/argument_data.json'))
    num_jobs = args.num_jobs
    if args.install_snakemake:
        install_snakemake()
    if args.mode == "all":
        run_all(args, num_jobs)
    elif args.mode == 'assembly':
        run_assembly(args,num_jobs)
    elif args.mode == "downstream":
        run_downstream(args, num_jobs)