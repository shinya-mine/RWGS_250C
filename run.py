#!/usr/bin/env python
import numpy as np
from datetime import datetime as dt
import sys, os
sys.dont_write_bytecode = True
import subprocess
import conditions, gen_data
import warnings
warnings.filterwarnings('ignore')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-s', '--job_size', type=int, default=50)
args = vars(parser.parse_args())

job_size = args['job_size']

desc_path = 'data/Descriptors.xlsx'
if os.path.exists(desc_path):
    print(f'Confirmed the existence of {desc_path}.')
else:
    subprocess.call("mv data/Descriptors_* data/Descriptors.xlsx", shell=True)
    print(f'{desc_path} has been created.')

condition = conditions.calc_condition()
computer = condition['computer']
core, hyper_threading = condition['core'], condition['hyper_threading']
Reaction = condition['Reaction']
Search_method = condition['Search_method']
data_cols = conditions.data_columns(condition)

if Search_method == 'all':
    job_title = 'A'
elif Search_method == 'all':
    job_title = 'L'

if hyper_threading == True:
    threads = 2 * core
else:
    threads = 1 * core

add_model = condition['add_model']

num_elem_wt = gen_data.elem_wt(condition)
print('num_elem_wt =', num_elem_wt)

def divisors_list(num):
    divisors = []
    for i in range(1, num+1):
        if num % i == 0:
            divisors.append(i)
    return divisors

def num_jobs(job_num_list):
    calc_cost = []
    for i in range(len(job_num_list)):
        calc_cost.append(num_elem_wt/job_num_list[i])
        if calc_cost[i] > job_size:
            num_of_jobs = job_num_list[i]
    return num_of_jobs

job_num_list = divisors_list(num_elem_wt)
if len(job_num_list) == 1:
    num_of_jobs = job_num_list[0]
else:
    num_of_jobs = int(num_jobs(job_num_list))
print('num_of_jobs =', num_of_jobs)
print('job_size =', int(num_elem_wt/num_of_jobs))

subprocess.call("mkdir job log out err", shell=True)

if Reaction == 'rwgs_250':
    Reaction = 'RWGS'
else:
    pass

if computer == 'A':
    for i, j in enumerate(zip(
        np.arange(0, num_elem_wt+(num_elem_wt/num_of_jobs), num_elem_wt/num_of_jobs, dtype=int),
        np.arange(num_elem_wt/num_of_jobs, num_elem_wt+(num_elem_wt/num_of_jobs), num_elem_wt/num_of_jobs, dtype=int)), 1):
        with open(f'run{computer}_{i}.sh', 'w') as f:
            print(
f'#!/bin/bash\n\
#============ PBS Options ============\n\
#QSUB -q gr10414a\n\
#QSUB -ug gr10414\n\
#QSUB -W 336:00\n\
#QSUB -A p=1:t={core}:c={core}:m=92160M\n\
#QSUB -N {Reaction}_{job_title}_p{add_model}_{i}\n\
#QSUB -r n\n\
#QSUB -o out/out_{i}\n\
#QSUB -e err/err_{i}\n\
#============ Shell Script ============\n\
ulimit -s unlimited\n\
cd $PBS_O_WORKDIR\n\
\n\
aprun -n $QSUB_PROCS -d $QSUB_THREADS -N $QSUB_PPN python {Search_method}_search.py --from {j[0]} --to {j[1]} --workers {core} --data_dump data_dump --cat_dump cat_dump --csv_name cand_{i} > log/log_cand_{i}',
            file=f
            )
        subprocess.call(f"qsub run{computer}_{i}.sh", shell=True)

elif computer == 'B':
    for i, j in enumerate(zip(
        np.arange(0, num_elem_wt+(num_elem_wt/num_of_jobs), num_elem_wt/num_of_jobs, dtype=int),
        np.arange(num_elem_wt/num_of_jobs, num_elem_wt+(num_elem_wt/num_of_jobs), num_elem_wt/num_of_jobs, dtype=int)), 1):
        with open(f'run{computer}_{i}.sh', 'w') as f:
            print(
f'#!/bin/bash\n\
#============ PBS Options ============\n\
#QSUB -q gr10414b\n\
#QSUB -ug gr10414\n\
#QSUB -W 336:00\n\
#QSUB -A p=1:t={threads}:c={core}:m=122880M\n\
#QSUB -N {Reaction}_{job_title}_p{add_model}_{i}\n\
#QSUB -r n\n\
#QSUB -o out/out_{i}\n\
#QSUB -e err/err_{i}\n\
#============ Shell Script ============\n\
#ulimit -s unlimited\n\
cd $PBS_O_WORKDIR\n\
\n\
mpiexec.hydra python {Search_method}_search.py --from {j[0]} --to {j[1]} --workers {core} --data_dump data_dump --cat_dump cat_dump --csv_name cand_{i} > log/log_cand_{i}',
            file=f
            )
        subprocess.call(f"qsub run{computer}_{i}.sh", shell=True)

elif computer == 'I':
    for i, j in enumerate(zip(
        np.arange(0, num_elem_wt+(num_elem_wt/num_of_jobs), num_elem_wt/num_of_jobs, dtype=int),
        np.arange(num_elem_wt/num_of_jobs, num_elem_wt+(num_elem_wt/num_of_jobs), num_elem_wt/num_of_jobs, dtype=int)), 1):
        with open(f'run{computer}_{i}.sh', 'w') as f:
            print(
f'#!/bin/bash\n\
#============ pjsub Options ===================\n\
#PJM -L "rscunit=ito-a"\n\
#PJM -L "rscgrp=ito-a-oc180150"\n\
#PJM -L "vnode=1"\n\
#PJM -L "vnode-core={core}"\n\
#PJM -L "elapse=168:00:00"\n\
#PJM -N {Reaction}_{i}_p{add_model}\n\
#PJM --no-stging\n\
#PJM --spath out/stat_{i}\n\
#PJM -o out/out_{i}\n\
#PJM -e err/err_{i}\n\
#PJM -S\n\
#PJM --restart\n\
#============ Shell Script ===================\n\
newgrp  oc180150\n\
LANG=C\n\
\n\
cd $PJM_O_WORKDIR\n\
export I_MPI_HYDRA_BOOTSTRAP_EXEC=pjrsh\n\
export I_MPI_HYDRA_HOST_FILE=$PJM_O_NODEINF\n\
export I_MPI_DEVICE=rdma\n\
export I_MPI_PERHOST={core}\n\
\n\
python {Search_method}_search.py --from {j[0]} --to {j[1]} --workers {core} --data_dump data_dump --cat_dump cat_dump --csv_name cand_{i} > log/log_cand_{i}',
            file=f
            )
        subprocess.call(f"pjsub run{computer}_{i}.sh", shell=True)

elif computer == 'local':
    for i, j in enumerate(zip(
        np.arange(0, num_elem_wt+(num_elem_wt/num_of_jobs), num_elem_wt/num_of_jobs, dtype=int),
        np.arange(num_elem_wt/num_of_jobs, num_elem_wt+(num_elem_wt/num_of_jobs), num_elem_wt/num_of_jobs, dtype=int)), 1):
        subprocess.call(f"python {Search_method}_search.py --from {j[0]} --to {j[1]} --workers {core} --data_dump data_dump --cat_dump cat_dump --csv_name cand_{i} > log/log_cand_{i}", shell=True)

subprocess.call(f"mv run{computer}_* job/", shell=True)