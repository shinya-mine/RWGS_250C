#!/usr/bin/env python
import numpy as np
import sys, os
sys.dont_write_bytecode = True
import conditions
import subprocess
import warnings
warnings.filterwarnings('ignore')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-id', '--job_id', type=int, default=1)
parser.add_argument('-n', '--num_of_jobs', type=int, default=100)
args = vars(parser.parse_args())
Start_job_id = args['job_id']
Num_of_qdel = args['num_of_jobs'] * 1.5

condition = conditions.calc_condition()
date, Reaction, computer = condition['date'], condition['Reaction'], condition['computer']

if computer == 'A' or computer == 'B' or computer == 'Shimizu':
	for job_id in np.arange(Start_job_id, Start_job_id+Num_of_qdel+1, 1, dtype=int):
		subprocess.call(f"qdel {job_id}", shell=True)
elif computer == 'I':
	for job_id in np.arange(Start_job_id, Start_job_id+Num_of_qdel+1, 1, dtype=int):
		subprocess.call(f"pjdel {job_id}", shell=True)

subprocess.call("rm -rf dump out err log job cand Figures", shell=True)
subprocess.call("rm -rf dump out err log job cand Figures", shell=True)
