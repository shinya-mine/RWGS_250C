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
import os.path
sys.dont_write_bytecode = True
from pathlib import Path

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
parser.add_argument('--data_dump', type=str, default='data_dump')
parser.add_argument('--cat_dump', type=str, default='cat_dump')
parser.add_argument('--csv_name', type=str, default="cand_0")

args = vars(parser.parse_args())

num_workers = args['workers']
idx_from = args['from']
idx_to = args['to']
csv_name = args['csv_name']
data_dump_name = args['data_dump']
cat_dump_name = args['cat_dump']

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

def EI(mu, sigma, cur_max):
    Z = (mu - cur_max) / sigma
    ei = (mu - cur_max) * norm.cdf(Z) + sigma*norm.pdf(Z) # pdf == Probability Density Function
    return ei

def encode(comps, ele):
    if cand_add_num == 2:
        return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}"
    elif cand_add_num == 3:
        return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}"
    elif cand_add_num == 4:
        return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}"
    elif cand_add_num == 5:
        return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}"
    elif cand_add_num == 6:
        return f"{ele[0]} {comps[0]}, {ele[1]} {comps[1]}, {ele[2]} {comps[2]}, {ele[3]} {comps[3]}, {ele[4]} {comps[4]}, {ele[5]} {comps[5]}"

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

def calc_ei(elem_wt):
    if add_model == 0 and cand_add_num != 0:
        elem_pattern = sum([d * w for d, w in zip(add_vect, elem_wt)])
    elif add_model == 1 and cand_add_num != 0:
        elem_com_pattern = sum([d * w for d, w in zip(add_com, elem_wt)])
        elem_vect_pattern = sum([d * w for d, w in zip(add_vect, elem_wt)])
        elem_pattern = sparse.hstack([elem_com_pattern, elem_vect_pattern])
        del elem_vect_pattern
    elif add_model == 2 and cand_add_num != 0:
        elem_pattern = sum([d * w for d, w in zip(add_com, elem_wt)])
    else:
        print('Additive passed...')
    print(f'elem_pattern_A{add_model}:', elem_pattern.shape, format_size(elem_pattern.data.nbytes))
    mu = model.predict(elem_pattern)
    sigma = pred_std_csr(elem_pattern, model.estimators_, mu)
    idx_s = sigma == 0
    sigma[idx_s] = 1e-5
    scores = EI(mu, sigma, cur_max)
    if save_depth == -1:
        idx = scores.argsort()[::-1]
    else:
        idx = scores.argsort()[::-1][:save_depth]
    return pd.DataFrame({'ei': scores[idx]}, index=get_encode(idx, elem_wt))

print('# cpu cores', os.cpu_count())
print('num_workers', num_workers)

PATH1 = 'dump'
os.makedirs(PATH1, exist_ok = True)
data_dump_out = f'dump/{data_dump_name}.pkl'
cat_dump_out = f'dump/{cat_dump_name}.pkl'

if idx_from == 0:
    converted = gen_data.data_convert(condition)
    catal_cand = gen_data.make_catal_cand(condition, converted)
    with open(data_dump_out, 'wb') as f1:
        pickle.dump(converted, f1)
        print(f'data saved as {data_dump_out}')
    with open(cat_dump_out, 'wb') as f2:
        pickle.dump(catal_cand, f2)
        print(f'data saved as {cat_dump_out}')
    
    print('Creation of the dump file is complete !!!')
    myfile = Path(f'{PATH1}/dump_finish.txt')
    myfile.touch(exist_ok=True)
else:
    print('Now Waiting ...')
    dump_fin_path = f'{PATH1}/dump_finish.txt'
    while( not( os.path.isfile(dump_fin_path))):
        time.sleep(1)
    print('Checked output of dump files !!!')
    #os.remove(dump_fin_path)
    print(f'gen_data recovered from {data_dump_out}')
    with open(data_dump_out, 'rb') as f1:
        converted = pickle.load(f1)
    print(f'cand_data recovered from {cat_dump_out}')
    with open(cat_dump_out, 'rb') as f2:
        catal_cand = pickle.load(f2)

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

add_model = condition['add_model']
cand_add_num = condition['cand_add_num']
target_name = condition['target_name']

save_depth = condition['save_depth']
catal_comb, elem_wt = catal_cand['catal_comb'], catal_cand['elem_wt']

add_com = catal_cand['add_com']
add_vect = catal_cand['add_vect']

feat, target = converted['feat'], converted['target']
cur_max = target.max()

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

alpha = len(elem_wt)/total_num
print('alpha = ', alpha)
print(f'processing elem_wt {idx_from}-{idx_to} ({len(elem_wt)}/{total_num}, {alpha:3.1%})\n')

with timer(f'fitting surrogate : n_jobs=1'):
    model = ExtraTreesRegressor(n_estimators=100, random_state=1107, n_jobs=1)
    model.fit(feat, target)

with timer('prosess_pool', 1.0/alpha):
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        res = executor.map(calc_ei, elem_wt)
        print('EI value calculations were complete !')

PATH2 = 'cand'
os.makedirs(PATH2, exist_ok = True)
cand = pd.concat(res)
print(cand.shape)
cand.to_csv(f'cand/{csv_name}.csv')

if idx_to == total_num:
    K_cluster = condition['K_cluster']
    subprocess.call(f"python sum_k_clus.py -k {K_cluster}", shell=True)
    print('All calculations and K-cluster calculations have been successfully completed !!')
    
    print('ML analysis start...')
    PATH = 'Figures'
    os.makedirs(PATH, exist_ok = True)
    converted = analysis.analysis_data_convert(condition, condition['data_sheet_name'], use_models=[add_model], idx=None)
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
