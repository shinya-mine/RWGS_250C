#!/usr/bin/env python
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.model_selection import ShuffleSplit
from numba import jit
from humanfriendly import format_size
import scipy.sparse as sparse
from scipy.stats import norm
import os, sys
sys.dont_write_bytecode = True

import pickle
import subprocess
from contextlib import contextmanager
import time, datetime

import conditions, gen_data, analysis

import random
random.seed(1107)
np.random.seed(1107)

import warnings
warnings.filterwarnings('ignore')

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--workers', type=int, default=1)
parser.add_argument('--from', type=int, default=0)
parser.add_argument('--to', type=int, default=1000)
parser.add_argument('--pkl_name', type=str, default='catal_dump')
parser.add_argument('--csv_name', type=str, default="cand_0")

args = vars(parser.parse_args())

num_workers = args['workers']
idx_from = args['from']
idx_to = args['to']
csv_name = args['csv_name']
pkl_name = args['pkl_name']

@contextmanager
def timer(name, alpha=None):
    t0 = time.time()
    yield
    dt = time.time() - t0
    print(f'[{name}] done in {dt:.2f} s')
    if alpha:
        delta = datetime.timedelta(seconds=dt * alpha)
        print(f'(estimate {delta} s for processing all)')

def get_dir_size(path='.'):
    total = 0
    with os.scandir(path) as it:
        for entry in it:
            if entry.is_file():
                total += entry.stat().st_size
            elif entry.is_dir():
                total += get_dir_size(entry.path)
    return total

def get_size(path='.'):
    if os.path.isfile(path):
        return os.path.getsize(path)
    elif os.path.isdir(path):
        return get_dir_size(path)

condition = conditions.calc_condition()
Reaction = condition['Reaction']
print(Reaction)

def EI(mu, sigma, cur_min, cur_max):
    if Reaction == 'N2O' or Reaction == 'H2SCR' or Reaction == 'NH3SCR' or Reaction == 'CH4' or Reaction == 'EtOH_CO2' or Reaction == 'EtOH_CO':
        Z = (cur_min - mu) / sigma
        ei = (cur_min - mu) * norm.cdf(Z) + sigma*norm.pdf(Z)
    elif Reaction == 'rwgs_250' or Reaction == 'rwgs_250_1wt' or Reaction == 'rwgs_300' or Reaction == 'CH3OH':
        Z = (mu - cur_max) / sigma
        ei = (mu - cur_max) * norm.cdf(Z) + sigma*norm.pdf(Z)
    return ei

def encode(comps, ele):
    if CalT_num == 0:
        if cand_pgm_num == 0 and cand_add_num == 2 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]}"
        elif cand_pgm_num == 0 and cand_add_num == 2 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}"
        elif cand_pgm_num == 0 and cand_add_num == 3 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]}"
        elif cand_pgm_num == 0 and cand_add_num == 3 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}"
        elif cand_pgm_num == 0 and cand_add_num == 4 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]}"
        elif cand_pgm_num == 0 and cand_add_num == 4 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}"
        elif cand_pgm_num == 0 and cand_add_num == 5 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]}"
        elif cand_pgm_num == 0 and cand_add_num == 5 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}"
        elif cand_pgm_num == 0 and cand_add_num == 6 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]}"
        elif cand_pgm_num == 0 and cand_add_num == 6 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}"
        
        elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]}"
        elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}"
        elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]}"
        elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}"
        elif cand_pgm_num == 1 and cand_add_num == 4 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]}"
        elif cand_pgm_num == 1 and cand_add_num == 4 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}"
        elif cand_pgm_num == 1 and cand_add_num == 5 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]}"
        elif cand_pgm_num == 1 and cand_add_num == 5 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}"
        
        elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]}"
        elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}"
        elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]}"
        elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}"
        elif cand_pgm_num == 2 and cand_add_num == 4 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]}"
        elif cand_pgm_num == 2 and cand_add_num == 4 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}"
        elif cand_pgm_num == 2 and cand_add_num == 5 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]} {comps[6]}, {ele[7]}"
        elif cand_pgm_num == 2 and cand_add_num == 5 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]} {comps[6]}"
        
        elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]}"
        elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}"
        elif cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]}"
        elif cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}"
        elif cand_pgm_num == 3 and cand_add_num == 4 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]} {comps[6]}, {ele[7]}"
        elif cand_pgm_num == 3 and cand_add_num == 4 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]} {comps[6]}"
        elif cand_pgm_num == 3 and cand_add_num == 5 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]} {comps[6]}, {ele[7]} {comps[7]}, {ele[8]}"
        elif cand_pgm_num == 3 and cand_add_num == 5 and cand_supp_num == 0:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]} {comps[6]}, {ele[7]} {comps[7]}"
        
    elif CalT_num == 1:
        if cand_pgm_num == 0 and cand_add_num == 2 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]}, {ele[3]}"
        elif cand_pgm_num == 0 and cand_add_num == 3 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]}, {ele[4]}"
        elif cand_pgm_num == 0 and cand_add_num == 4 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]}, {ele[5]}"
        elif cand_pgm_num == 0 and cand_add_num == 5 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]}, {ele[6]}"
        elif cand_pgm_num == 0 and cand_add_num == 6 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]}, {ele[7]}"

        elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]}, {ele[4]}"
        elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]}, {ele[5]}"
        elif cand_pgm_num == 1 and cand_add_num == 4 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]}, {ele[6]}"
        elif cand_pgm_num == 1 and cand_add_num == 5 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]}, {ele[7]}"
        
        elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]}, {ele[5]}"
        elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]}, {ele[6]}"
        elif cand_pgm_num == 2 and cand_add_num == 4 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]}, {ele[7]}"
        elif cand_pgm_num == 2 and cand_add_num == 5 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]} {comps[6]}, {ele[7]}, {ele[8]}"
        
        elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]}, {ele[6]}"
        elif cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]}, {ele[7]}"
        elif cand_pgm_num == 3 and cand_add_num == 4 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]} {comps[6]}, {ele[7]}, {ele[8]}"
        elif cand_pgm_num == 3 and cand_add_num == 5 and cand_supp_num == 1:
            return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}, {ele[6]} {comps[6]}, {ele[7]} {comps[7]}, {ele[8]}, {ele[9]}"

def get_encode(idx, comps):
    return [encode(comps, ele) for ele in catal_comb[idx]]

def pred_std_csr(X, trees, predictions, min_variance=0.0):
    std = np.zeros(X.shape[0])
    for tree in trees:
        var_tree = tree.tree_.impurity[tree.apply(X)]
        var_tree[var_tree < min_variance] = min_variance
        mean_tree = tree.predict(X)
        std += var_tree + mean_tree ** 2
    std /= len(trees)
    std -= predictions ** 2.0
    std[std < 0.0] = 0.0
    std = std ** 0.5
    return std

def calc_ei(elem_wt, pgm_wt, add_wt):
    if pgm_model == 0 and cand_pgm_num != 0:
        pgm_pattern = sum([d * w for d, w in zip(pgm_vect, pgm_wt)])
    elif pgm_model == 1 and cand_pgm_num != 0:
        pgm_com_pattern = sum([d * w for d, w in zip(pgm_com, pgm_wt)])
        pgm_vect_pattern = sum([d * w for d, w in zip(pgm_vect, pgm_wt)])
        pgm_pattern = sparse.hstack([pgm_com_pattern, pgm_vect_pattern])
        del pgm_vect_pattern
    elif pgm_model == 2 and cand_pgm_num != 0:
        pgm_pattern = sum([d * w for d, w in zip(pgm_com, pgm_wt)])
    else:
        print('PGM passed...')

    if add_model == 0 and cand_add_num != 0:
        add_pattern = sum([d * w for d, w in zip(add_vect, add_wt)])
    elif add_model == 1 and cand_add_num != 0:
        add_com_pattern = sum([d * w for d, w in zip(add_com, add_wt)])
        add_vect_pattern = sum([d * w for d, w in zip(add_vect, add_wt)])
        add_pattern = sparse.hstack([add_com_pattern, add_vect_pattern])
        del add_vect_pattern
    elif add_model == 2 and cand_add_num != 0:
        add_pattern = sum([d * w for d, w in zip(add_com, add_wt)])
    else:
        print('Additive passed...')
    
    if supp_model == 0 and cand_supp_num != 0:
        supp_pattern = sum([d * w for d, w in zip(supp_vect, supp_wt)])
    elif supp_model == 1 and cand_supp_num != 0:
        supp_com_pattern = sum([d * w for d, w in zip(supp_com, supp_wt)])
        supp_vect_pattern = sum([d * w for d, w in zip(supp_vect, supp_wt)])
        supp_pattern = sparse.hstack([supp_com_pattern, supp_vect_pattern])
        del supp_vect_pattern
    elif supp_model == 2 and cand_supp_num != 0:
        supp_pattern = sum([d * w for d, w in zip(supp_com, supp_wt)])
    else:
        print('Support passed...')

    if cand_pgm_num != 0 and cand_add_num != 0 and cand_supp_num != 0 and CalT_num == 1:
        elem_pattern = sparse.hstack([pgm_pattern, add_pattern, supp_pattern, CalT_pattern])
    elif cand_pgm_num == 0 and cand_add_num != 0 and cand_supp_num != 0 and CalT_num == 1:
        elem_pattern = sparse.hstack([add_pattern, supp_pattern, CalT_pattern])
    elif cand_pgm_num == 0 and cand_add_num != 0 and cand_supp_num == 0 and CalT_num == 0:
        elem_pattern = add_pattern
    else:
        elem_pattern = sparse.hstack([pgm_pattern, add_pattern, supp_pattern])

    if cand_pgm_num != 0 and cand_add_num != 0 and cand_supp_num != 0:
        print(f'elem_pattern_P{pgm_model}A{add_model}S{supp_model}:', elem_pattern.shape, format_size(elem_pattern.data.nbytes))
    elif cand_pgm_num == 0 and cand_add_num != 0 and cand_supp_num != 0:
        print(f'elem_pattern_A{add_model}S{supp_model}:', elem_pattern.shape, format_size(elem_pattern.data.nbytes))
    elif cand_pgm_num != 0 and cand_add_num == 0 and cand_supp_num != 0:
        print(f'elem_pattern_P{pgm_model}S{supp_model}:', elem_pattern.shape, format_size(elem_pattern.data.nbytes))
    elif cand_pgm_num == 0 and cand_add_num != 0 and cand_supp_num == 0:
        print(f'elem_pattern_A{add_model}:', elem_pattern.shape, format_size(elem_pattern.data.nbytes))
        
    mu = model.predict(elem_pattern)
    sigma = pred_std_csr(elem_pattern, model.estimators_, mu)
    idx_s = sigma == 0
    sigma[idx_s] = 1e-5
    scores = EI(mu, sigma, cur_min, cur_max)
    if save_depth == -1:
        idx = scores.argsort()[::-1]
    else:
        idx = scores.argsort()[::-1][:save_depth]
    return pd.DataFrame({'ei': scores[idx]}, index=get_encode(idx, elem_wt))

print('# cpu cores', os.cpu_count())
print('num_workers', num_workers)

PATH1 = 'dump'
os.makedirs(PATH1, exist_ok = True)
data_pkl = f'dump/data_dump.pkl'
pkl_out = f'dump/{pkl_name}.pkl'

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


data_size = get_size('dump')*1e-6
print(int(data_size),'MB')
computer, core = condition['computer'], condition['core']

if computer == 'A' and data_size < 92160/core:
    print(int(data_size),'MB < ', int(92160/core), 'MB')
    print('The data size of this job is small enough for the RAM in the computer.')
    print('Continue calculation...')
elif computer == 'A' and data_size >= 92160/core:
    print(int(data_size),'MB > ', int(92160/core), 'MB')
    print('!!!!!!!!!! >>>>>>>>>> ERROR <<<<<<<<<< !!!!!!!!!!')
    print('The data size of this job is too large for the RAM installed in the computer.')
    print('Edit the condition.py file and change the calculation conditions.')
    print('Shut down the computation...')
    sys.exit()
elif computer == 'B' and data_size < 122880/core:
    print(int(data_size),'MB < ', int(122880/core), 'MB')
    print('The data size of this job is small enough for the RAM in the computer.')
    print('Continue calculation...')
elif computer == 'B' and data_size >= 122880/core:
    print(int(data_size),'MB > ', int(122880/core), 'MB')
    print('!!!!!!!!!! >>>>>>>>>> ERROR <<<<<<<<<< !!!!!!!!!!')
    print('The data size of this job is too large for the RAM installed in the computer.')
    print('Edit the condition.py file and change the calculation conditions.')
    print('Shut down the computation...')
    sys.exit()
elif computer == 'I' and data_size < 192000/core:
    print(int(data_size),'MB < ', int(192000/core), 'MB')
    print('The data size of this job is small enough for the RAM in the computer.')
    print('Continue calculation...')
elif computer == 'I' and data_size >= 192000/core:
    print(int(data_size),'MB > ', int(192000/core), 'MB')
    print('!!!!!!!!!! >>>>>>>>>> ERROR <<<<<<<<<< !!!!!!!!!!')
    print('The data size of this job is too large for the RAM installed in the computer.')
    print('Edit the condition.py file and change the calculation conditions.')
    print('Shut down the computation...')
    sys.exit()
else:
    print('!!!!!!!!!! >>>>>>>>>> ERROR <<<<<<<<<< !!!!!!!!!!')
    sys.exit()

pgm_model, add_model, supp_model = condition['pgm_model'], condition['add_model'], condition['supp_model']
cand_pgm_num, cand_add_num, cand_supp_num, CalT_num  = condition['cand_pgm_num'], condition['cand_add_num'], condition['cand_supp_num'], condition['CalT_num']

target_temp = condition['target_temp']
if Reaction == 'rwgs_250' or Reaction == 'rwgs_300' or Reaction == 'CH3OH':
    target_name = condition['target_name']
else:
    pass

save_depth = condition['save_depth']
catal_comb, elem_wt = catal_cand['catal_comb'], catal_cand['elem_wt']
if CalT_num != 0:
    CalT_pattern = catal_cand['CalT_com'][0]
else:
    pass

pgm_wt = elem_wt[:, 0: cand_pgm_num]
add_wt = elem_wt[:, cand_pgm_num: cand_pgm_num+cand_add_num]
if cand_supp_num == 1:
    supp_wt = elem_wt[:, -cand_supp_num]
else:
    pass

if cand_pgm_num == 0 and cand_supp_num != 0:
    add_com, supp_com =  catal_cand['add_com'], catal_cand['supp_com']
    add_vect, supp_vect =  catal_cand['add_vect'], catal_cand['supp_vect']

elif cand_add_num == 0 and cand_supp_num != 0:
    pgm_com,  supp_com = catal_cand['pgm_com'], catal_cand['supp_com']
    pgm_vect, supp_vect = catal_cand['pgm_vect'], catal_cand['supp_vect']

elif  cand_pgm_num == 0 and cand_supp_num == 0:
    add_com, add_vect = catal_cand['add_com'], catal_cand['add_vect']

elif  cand_add_num == 0 and cand_supp_num == 0:
    pgm_com, pgm_vect = catal_cand['pgm_com'], catal_cand['pgm_vect']

else:
    pgm_com, add_com, supp_com = catal_cand['pgm_com'], catal_cand['add_com'], catal_cand['supp_com']
    pgm_vect, add_vect, supp_vect = catal_cand['pgm_vect'], catal_cand['add_vect'], catal_cand['supp_vect']

feat, target = converted['feat'], converted['target']
cur_max, cur_min = target.max(), target.min()

if idx_from == 0:
    print('-' * 50)
    model  = ExtraTreesRegressor(n_estimators = 100, random_state = 1107, n_jobs = -1)
    cvf = ShuffleSplit(n_splits = 100, random_state = 1107, test_size = 0.2)
    analysis.crossvalid(feat, target, model, cvf)
    print('-' * 50)
else:
    pass

total_num = len(elem_wt)
elem_wt = elem_wt[idx_from:idx_to]
pgm_wt = pgm_wt[idx_from:idx_to]
add_wt = add_wt[idx_from:idx_to]

alpha = len(elem_wt)/total_num
print('alpha = ', alpha)
print(f'processing elem_wt {idx_from}-{idx_to} ({len(elem_wt)}/{total_num}, {alpha:3.1%})\n')

with timer(f'fitting surrogate : n_jobs=1'):
    model = ExtraTreesRegressor(n_estimators=100, random_state=1107, n_jobs=1)
    model.fit(feat, target)

with timer('prosess_pool', 1.0/alpha):
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        res = executor.map(calc_ei, elem_wt, pgm_wt, add_wt)
        print('EI value calculations were complete !')

PATH2 = 'cand'
os.makedirs(PATH2, exist_ok = True)
cand = pd.concat(res)
print(cand.shape)
cand.to_csv(f'cand/{csv_name}.csv')

if idx_to == total_num:
    if Reaction == 'rwgs_250_1wt':
        test_data_sheet = condition['data_sheet_name']
        train_data_sheet = condition['train_data_sheet_name']

        converted_test = analysis.analysis_data_convert(condition, test_data_sheet, use_models=[pgm_model, add_model, supp_model], idx=None)
        converted_train = analysis.analysis_data_convert(condition, train_data_sheet, use_models=[pgm_model, add_model, supp_model], idx=None)
        feat_test, target_test = converted_test['feat'], converted_test['target']
        feat_train, target_train = converted_train['feat'], converted_train['target']
        feat_all = pd.concat([feat_train, feat_test], axis=0).reset_index(drop=True)
        target_all = pd.concat([target_train, target_test], axis=0).reset_index(drop=True)
        feat_cols = feat_test.columns
        model  = ExtraTreesRegressor(n_estimators = 100, random_state = 1107, n_jobs = -1)
        cvf = ShuffleSplit(n_splits = 100, random_state = 1107, test_size = 0.2)
        #cvf = KFold(n_splits=5, shuffle=True)

        print('1. CV: test_data (New dataset)')
        analysis.crossvalid(feat_test, target_test, model, cvf)
        print('2. CV: training_data (396 data points)')
        analysis.crossvalid(feat_train, target_train, model, cvf)
        print('3. CV: ALL_data')
        analysis.crossvalid(feat_all, target_all, model, cvf)
        print('4. CV: training_data + test_data (Only a part of test_data was used for verification)')
        analysis.crossvalid_concat(feat_train, feat_test, target_train, target_test, model, cvf)
    else:
        pass
    K_cluster = condition['K_cluster']
    subprocess.call(f"python sum_k_clus.py -k {K_cluster}", shell=True)
    print('All calculations and K-cluster calculations have been successfully completed !!')
    
    extractions = condition['extractions']
    if len(extractions) != 0:
        for i in range(len(extractions)):
            subprocess.call(f"python extract.py -ex {extractions[i]}", shell=True)
            print(f'K-clustering was attempted by extracting candidate catalysts that contained {extractions[i]}.')
    else:
        print('Extractions were NOT performed...')
    
    print('ML analysis start...')
    PATH = 'Figures'
    os.makedirs(PATH, exist_ok = True)

    if Reaction == 'rwgs_250_1wt':
        converted_test = analysis.analysis_data_convert(condition, condition['data_sheet_name'], use_models=[pgm_model, add_model, supp_model], idx=None)
        converted_train = analysis.analysis_data_convert(condition, condition['train_data_sheet_name'], use_models=[pgm_model, add_model, supp_model], idx=None)
        feat_test, target_test = converted_test['feat'], converted_test['target']
        feat_train, target_train = converted_train['feat'], converted_train['target']
        feat = pd.concat([feat_train, feat_test], axis=0).reset_index(drop=True)
        target = pd.concat([target_train, target_test], axis=0).reset_index(drop=True)
    else:
        converted = analysis.analysis_data_convert(condition, condition['data_sheet_name'], use_models=[pgm_model, add_model, supp_model], idx=None)
        feat, target = converted['feat'], converted['target']

    feat_cols = feat.columns
    model = ExtraTreesRegressor(n_estimators = 100, random_state = 1107, n_jobs = -1)
    analysis.shap_summary_plot(condition, model, feat, target, save=True)
    analysis.one_shot_plot(condition, feat, target, model, Reaction, test_size=0.1, random_state=1107, save=True)
    analysis.plot_importance(condition, model, labels=feat_cols, topk=20, save=True)
    analysis.loo_plot(condition, feat, target, feat_cols, model, save=True)
    print('All ML analysis have been successfully completed !!')
else:
    print(f'The calculation {csv_name} has been successfully completed !!')
