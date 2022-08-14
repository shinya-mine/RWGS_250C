#!/usr/bin/env python
import numpy as np
import subprocess

import sys, os
sys.dont_write_bytecode = True
import conditions

import random
random.seed(1107)
np.random.seed(1107)

import warnings
warnings.filterwarnings('ignore')

PATH = 'extraction'
os.makedirs(PATH, exist_ok = True)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-ex', '--extract', type=str, default='H')
parser.add_argument('-wt', '--extract_wt', type=float, default=-1)
parser.add_argument('-r', '--extract_range', type=str, default='equal')

args = vars(parser.parse_args())

extract = args['extract']
extract_wt = args['extract_wt']
extract_range = args['extract_range']

### Basic Settings ###
condition = conditions.calc_condition()
computer = condition['computer']
core = condition['core']

threads = 1 * core
"""
hyper_threading = condition['hyper_threading']
if hyper_threading == True:
    threads = 2 * core
else:
    threads = 1 * core
"""

if computer == 'A':
    with open(f'run{computer}_ex.sh', 'w') as f:
        print(
f'#!/bin/bash\n\
#============ PBS Options ============\n\
#QSUB -q gr10414a\n\
#QSUB -ug gr10414\n\
#QSUB -W 336:00\n\
#QSUB -A p=1:t={threads}:c={core}:m=92160M\n\
#QSUB -N ex_{extract}\n\
#QSUB -r n\n\
#QSUB -o out/out_ex_{extract}\n\
#QSUB -e err/err_ex_{extract}\n\
#============ Shell Script ============\n\
ulimit -s unlimited\n\
cd $PBS_O_WORKDIR\n\
\n\
aprun -n $QSUB_PROCS -d $QSUB_THREADS -N $QSUB_PPN python extract.py -ex {extract} -wt {extract_wt} -r {extract_range} > log/log_ex_{extract}',
            file=f
            )
    subprocess.call(f"qsub run{computer}_ex.sh", shell=True)

elif computer == 'B':
    with open(f'run{computer}_ex.sh', 'w') as f:
        print(
f'#!/bin/bash\n\
#============ PBS Options ============\n\
#QSUB -q gr10414b\n\
#QSUB -ug gr10414\n\
#QSUB -W 336:00\n\
#QSUB -A p=1:t={threads}:c={core}:m=122880M\n\
#QSUB -N ex_{extract}\n\
#QSUB -r n\n\
#QSUB -o out/out_ex_{extract}\n\
#QSUB -e err/err_ex_{extract}\n\
#============ Shell Script ============\n\
#ulimit -s unlimited\n\
cd $PBS_O_WORKDIR\n\
\n\
mpiexec.hydra python extract.py -ex {extract} -wt {extract_wt} -r {extract_range} > log/log_ex_{extract}',
            file=f
            )
    subprocess.call(f"qsub run{computer}_ex.sh", shell=True)

elif computer == 'I':
    with open(f'run{computer}_ex.sh', 'w') as f:
        print(
f'#!/bin/bash\n\
#============ pjsub Options ===================\n\
#PJM -L "rscunit=ito-a"\n\
#PJM -L "rscgrp=ito-a-oc180150"\n\
#PJM -L "vnode=1"\n\
#PJM -L "vnode-core={core}"\n\
#PJM -L "elapse=168:00:00"\n\
#PJM -N ex_{extract}\n\
#PJM --no-stging\n\
#PJM --spath out/stat_ex_{extract}\n\
#PJM -o out/out_ex_{extract}\n\
#PJM -e err/err_ex_{extract}\n\
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
python extract.py -ex {extract} -wt {extract_wt} -r {extract_range} > log/log_ex_{extract}',
            file=f
            )
    subprocess.call(f"pjsub run{computer}_ex.sh", shell=True)
subprocess.call(f"mv run{computer}_ex.sh job/", shell=True)