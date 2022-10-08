import os, subprocess, re, time, sys, argparse
from collections import deque
from pathlib import Path

full_path = Path(os.path.dirname(os.path.realpath(__file__))).parents[0] #fastq processing folder
parser = argparse.ArgumentParser(description='A script to estimate biological composition of sequencing library from obtained fastq file.') #argparser object to provide command-line functionality
req_arg_grp = parser.add_argument_group('required arguments') #add optional argument group
req_arg_grp.add_argument('-f', '--fastq', metavar='\b', help = 'Path to the folder containing fastq files to process (should be relative to ~/ and found in it). All files from folder and subfolders will be included in the analysis.', required=True) #required fastq file argument
parser.add_argument('-n', '--num_files', metavar='\b', help = 'Number of samples to be processed in-parallel on different nodes of the cluster (8 by default)', default=8, required=False)
parser.add_argument('-o', '--output_dir', metavar='\b', help = 'Folder where the fastq_screen output should be copied after processing is finished.', required=False, default=None)

if len(sys.argv)==1: #if no command-line arguments provided - display help and stop script excecution
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args() #args list from command-line input
home_path = os.path.expanduser('~/')
fastq_path = args.fastq #path to fastq files from command-line argument
report_time = time.strftime("%m_%d_%Y")
output_path = args.output_dir

os.chdir(home_path)
if os.path.isdir(f"{home_path}/{fastq_path}"):
    print("Fastq path valid")
    os.system(f"touch {home_path}/fastq_screen_lock.file")
else:
    sys.exit("Fastq path invalid. \n Please provide path to fastq files under user's home directory.")

file_list = []
for (root,dirs,files) in os.walk(fastq_path, topdown=True): # iterate over all subdir under fastq_path and add path to fastq files to list
    for name in files: 
        if 'fastq.gz' in name: file_list.append(os.path.join(root, name)) 

if len(file_list) == 0: sys.exit('No fastq files were found in the specified folder and subfolders.')
    

os.chdir(f'{full_path}/fastq_processing/subprocess')
if not os.path.isdir(f'{home_path}/fastq_screen_output/'): os.mkdir(f'{home_path}/fastq_screen_output/')
    

def submit_job(file_path,que):
    print("\n", file_path,
    f'{len(que)} samples left to process from total of {len(file_list)}', 
    f'qsub -F {file_path} {home_path} run_fastq_screen.sh', sep='\n')
    subprocess.check_call(['qsub', '-F', f'{home_path}/{file_path} {home_path}', 'run_fastq_screen.sh'])

sample_limit = int(args.num_files) #number of samples to be processed at the same time on different nodes
sleep_time = 180 #Sleep time in seconds to wait before adding new samples to que
sample_que = deque(file_list) #Create a que of samples
processing_over=False
while len(sample_que) > 0: #Until there are samples to be processed
    sample_path = sample_que.pop() #Pop right-most sample id from que (que is shortened by on element
    while len([re.search(".fastq.gz$",file).group(0) for file in os.listdir(f'{home_path}/fastq_screen_output/') if re.search(".fastq.gz$",file)]) >= sample_limit:
        print(f'Job submission halted for {sleep_time} seconds - {sample_limit} samples being processed at the moment.')
        time.sleep(sleep_time) #Wait for {sleep_time seconds}
    os.system(f'touch {home_path}/fastq_screen_output/$(basename {sample_path})')
    submit_job(sample_path,sample_que)

#if the flow of excecution is here then all samples in the que are submitted for processing
while len([report for report in os.listdir(f'{home_path}/fastq_screen_output/') if '_screen.html' in report]) < len(file_list): #waiting for the processing of final submitted samples to complete
    print(f'Waiting for the processing to finish - next check in {sleep_time} seconds')
    time.sleep(sleep_time)

os.system(f"rm -f {full_path}/fastq_processing/subprocess/read_library*")
if output_path:
    os.system(f'mv -n {home_path}/fastq_screen_output/* {output_path}/')

os.system(f'rm {home_path}/fastq_screen_lock.file')