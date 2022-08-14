#!/usr/bin/env python
import numpy as np
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.model_selection import ShuffleSplit

import os, sys
sys.dont_write_bytecode = True

import pickle
import subprocess

import conditions, gen_data, analysis

import random
random.seed(1107)
np.random.seed(1107)

import warnings
warnings.filterwarnings('ignore')

desc_path = 'data/Descriptors.xlsx'
if os.path.exists(desc_path):
    print(f'Confirmed the existence of {desc_path}.')
else:
    subprocess.call("mv data/Descriptors_* data/Descriptors.xlsx", shell=True)
    print(f'{desc_path} has been created.')

condition = conditions.calc_condition()
Reaction = condition['Reaction']
print(Reaction)

PATH1 = 'dump'
os.makedirs(PATH1, exist_ok = True)
data_pkl = f'dump/data_dump.pkl'
pkl_out = f'dump/catal_dump.pkl'

if os.path.exists(data_pkl) and os.path.exists(pkl_out):
    print(f'gen_data recovered from {data_pkl}')
    with open(data_pkl, 'rb') as f1:
        converted = pickle.load(f1)
    print(f'cand_data recovered from {pkl_out}')
    with open(pkl_out, 'rb') as f2:
        catal_cand = pickle.load(f2)
else:
    converted = gen_data.data_convert(condition)
    catal_cand = gen_data.make_catal_cand(condition, converted)
    with open(data_pkl, 'wb') as f1:
        pickle.dump(converted, f1)
        print(f'data saved as {data_pkl}')
    with open(pkl_out, 'wb') as f2:
        pickle.dump(catal_cand, f2)
        print(f'data saved as {pkl_out}')

add_model = condition['add_model']
cand_add_num = condition['cand_add_num']
target_name = condition['target_name']

save_depth = condition['save_depth']
catal_comb, elem_wt = catal_cand['catal_comb'], catal_cand['elem_wt']
add_wt = elem_wt[:, 0: cand_add_num]

add_com = catal_cand['add_com']
add_vect = catal_cand['add_vect']

feat, target = converted['feat'], converted['target']
cur_max = target.max()

print('-' * 50)
model  = ExtraTreesRegressor(n_estimators = 100, random_state = 1107, n_jobs = -1)
cvf = ShuffleSplit(n_splits = 100, random_state = 1107, test_size = 0.2)
error = analysis.crossvalid(feat, target, model, cvf)
print(error)
print('-' * 50)
subprocess.call(f"rm -rf dump", shell=True)