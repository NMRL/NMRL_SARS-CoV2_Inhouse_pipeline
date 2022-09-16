from subscripts.covipipe_classes import run_all, run_assembly, run_downstream
from subscripts.covipipe_utilities import covipipe_housekeeper as hk 

"""
This is a wrapper script of COVIPIPE pipeline.
Date: 2022-06-13
Version: 0.0
"""

if __name__ == "__main__":
    args = hk.parse_arguments(hk.read_json_dict('./config_files/json/argument_data.json'))
    num_jobs = args.num_jobs
    if args.install_snakemake:
        hk.install_snakemake()
    if args.mode == "all":
        run_all(args, num_jobs)
    elif args.mode == 'assembly':
        run_assembly(args,num_jobs)
    elif args.mode == "downstream":
        run_downstream(args)