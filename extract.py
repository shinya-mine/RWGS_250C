#!/usr/bin/env python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from skopt.learning import ExtraTreesRegressor as opt_ETR
from skopt.learning import RandomForestRegressor as opt_RFR
from sklearn.cluster import KMeans

import sys, os
sys.dont_write_bytecode = True
import conditions
import gen_data

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
date = condition['date']
Reaction = condition['Reaction']
ML_model = condition['ML_model']
Search_method = condition['Search_method']

catal_number, K_cluster = condition['catal_number'], condition['K_cluster']

data_file_name = f'data/{date}_{Reaction}_dataset.xlsx'
data_sheet_name = condition['data_sheet_name']
skip_rows = condition['skip_rows']
desc_file_name = f'data/Descriptors.xlsx'
elem_desc_sheet = 'Descriptors_elem'
supp_desc_sheet = 'Descriptors_supp'

data_cols = conditions.data_columns(Reaction, condition)
elem_cols, wt_cols = data_cols['elem'], data_cols['wt']
data_labels = data_cols['labels']

if Reaction == 'rwgs_250' or Reaction == 'rwgs_300' or Reaction == 'CH3OH':
    target_name = condition['target_name']
elif Reaction == 'N2O' or Reaction == 'H2SCR' or Reaction == 'NH3SCR' or Reaction == 'CH4':
    target_temp = condition['target_temp']

pgm_model, add_model, supp_model = condition['pgm_model'], condition['add_model'], condition['supp_model']
cand_pgm_num, cand_add_num, cand_supp_num, CalT_num  = condition['cand_pgm_num'], condition['cand_add_num'], condition['supp_num'], condition['CalT_num']
fix_pgm_num, fix_add_num, fix_supp_num = condition['fix_pgm_num'], condition['fix_add_num'], condition['fix_supp_num']

if CalT_num == 1:
    CalT_list = condition['CalT_list']

### Define feat and target ###
converted = gen_data.data_convert(condition)
feat, target = converted['feat'], converted['target']
pgm_desc, add_desc, supp_desc = converted['pgm_desc'], converted['add_desc'], converted['supp_desc']
pgm_list, add_list, supp_list = list(pgm_desc.index), list(add_desc.index), list(supp_desc.index)
pgm_list.remove('H')
add_list.remove('H')

### USE Descriptors ###
desc = conditions.use_desc()
pgm_use_desc, add_use_desc, supp_use_desc = desc['pgm'], desc['add'], desc['supp']

### INPUT FILES (About Cand.) ###
cand_data = f'results/{date}_{Reaction}_cand_sum_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}.csv'
cand_elem_cols, cand_wt_cols, cand_labels = data_cols['cand_elem'], data_cols['cand_wt'], data_cols['cand_labels']

### Top catalysts OUTPUT Files ###
out_csv_file_top = f'{PATH}/{date}_{Reaction}_cand_sum_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}_{extract}.csv'
fig_title_top = f'{date}_{Reaction}_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}_top{catal_number}_{extract}'
out_fig_name_top = f'{PATH}/{date}_{Reaction}_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}_top{catal_number}_{extract}.png'

### K-mean catalysts OUTPUT Files ###
out_csv_file_K = f"{PATH}/{date}_{Reaction}_K{K_cluster}_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}_{extract}.csv"
fig_title_K = f'{date}_{Reaction}_K{K_cluster}_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}_{extract}'
out_fig_name_K = f"{PATH}/{date}_{Reaction}_K{K_cluster}_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}_{extract}.png"

cand = gen_data.cand_separation(cand_data)

if extract in pgm_list:
    print(f'Extract {extract}-containing catalyst candidates in PGMs...')
    if cand_pgm_num == 0:
        print('ERROR: There are no PGMs in the candidate catalysts!')
        
    elif cand_pgm_num == 1:
        cand_pgm1 = cand[cand['PGM1'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] != 0.0]
        cand = pd.DataFrame()
        cand = cand_pgm1
        
    elif cand_pgm_num == 2:
        cand_pgm1 = cand[cand['PGM1'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] != 0.0]
            
        cand_pgm2 = cand[cand['PGM2'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_pgm2 = cand_pgm2[cand_pgm2['PGM2_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_pgm2 = cand_pgm2[cand_pgm2['PGM2_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_pgm2 = cand_pgm2[cand_pgm2['PGM2_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_pgm2 = cand_pgm2[cand_pgm2['PGM2_wt%'] != 0.0]
        cand = pd.DataFrame()
        cand = pd.concat([cand_pgm1, cand_pgm2])
        
    elif cand_pgm_num == 3:
        cand_pgm1 = cand[cand['PGM1'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_pgm1 = cand_pgm1[cand_pgm1['PGM1_wt%'] != 0.0]
            
        cand_pgm2 = cand[cand['PGM2'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_pgm2 = cand_pgm2[cand_pgm2['PGM2_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_pgm2 = cand_pgm2[cand_pgm2['PGM2_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_pgm2 = cand_pgm2[cand_pgm2['PGM2_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_pgm2 = cand_pgm2[cand_pgm2['PGM2_wt%'] != 0]
            
        cand_pgm3 = cand[cand['PGM3'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_pgm3 = cand_pgm3[cand_pgm3['PGM3_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_pgm3 = cand_pgm3[cand_pgm3['PGM3_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_pgm3 = cand_pgm3[cand_pgm3['PGM3_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_pgm3 = cand_pgm3[cand_pgm3['PGM3_wt%'] != 0.0]
        cand = pd.DataFrame()
        cand = pd.concat([cand_pgm1, cand_pgm2, cand_pgm3])
    
elif extract in add_list:
    print(f'Extract {extract}-containing catalyst candidates in Additivess...')
    if cand_add_num == 0:
        print('ERROR: There are no Additives in the candidate catalysts!')
    
    elif cand_add_num == 1:
        cand_Ad1 = cand[cand['Ad1'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] != 0.0]
        cand = pd.DataFrame()
        cand = cand_Ad1
    
    elif cand_add_num == 2:
        cand_Ad1 = cand[cand['Ad1'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] != 0.0]
            
        cand_Ad2 = cand[cand['Ad2'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] != 0]
        cand = pd.DataFrame()
        cand = pd.concat([cand_Ad1, cand_Ad2])
        
    elif cand_add_num == 3:
        cand_Ad1 = cand[cand['Ad1'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] != 0.0]
            
        cand_Ad2 = cand[cand['Ad2'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] != 0]
        
        cand_Ad3 = cand[cand['Ad3'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] != 0.0]
        cand = pd.DataFrame()
        cand = pd.concat([cand_Ad1, cand_Ad2, cand_Ad3])
    
    elif cand_add_num == 4:
        cand_Ad1 = cand[cand['Ad1'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] != 0.0]
            
        cand_Ad2 = cand[cand['Ad2'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] != 0]
        
        cand_Ad3 = cand[cand['Ad3'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] != 0.0]
        
        cand_Ad4 = cand[cand['Ad4'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad4 = cand_Ad4[cand_Ad4['Ad4_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad4 = cand_Ad4[cand_Ad4['Ad4_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad4 = cand_Ad4[cand_Ad4['Ad4_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad4 = cand_Ad4[cand_Ad4['Ad4_wt%'] != 0.0]
        cand = pd.DataFrame()
        cand = pd.concat([cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4])
    
    elif cand_add_num == 5:
        cand_Ad1 = cand[cand['Ad1'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad1 = cand_Ad1[cand_Ad1['Ad1_wt%'] != 0.0]
        
        cand_Ad2 = cand[cand['Ad2'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad2 = cand_Ad2[cand_Ad2['Ad2_wt%'] != 0]
        
        cand_Ad3 = cand[cand['Ad3'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad3 = cand_Ad3[cand_Ad3['Ad3_wt%'] != 0.0]
        
        cand_Ad4 = cand[cand['Ad4'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad4 = cand_Ad4[cand_Ad4['Ad4_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad4 = cand_Ad4[cand_Ad4['Ad4_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad4 = cand_Ad4[cand_Ad4['Ad4_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad4 = cand_Ad4[cand_Ad4['Ad4_wt%'] != 0.0]
        
        cand_Ad5 = cand[cand['Ad5'] == extract]
        if extract_wt > 0.0 and extract_range == 'small':
            cand_Ad5 = cand_Ad5[cand_Ad5['Ad5_wt%'] <= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'large':
            cand_Ad5 = cand_Ad5[cand_Ad5['Ad5_wt%'] >= extract_wt]
        elif extract_wt > 0.0 and extract_range == 'equal':
            cand_Ad5 = cand_Ad5[cand_Ad5['Ad5_wt%'] == extract_wt]
        elif extract_wt == -1:
            cand_Ad5 = cand_Ad5[cand_Ad5['Ad5_wt%'] != 0.0]
        cand = pd.DataFrame()
        cand = pd.concat([cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5])
        
elif extract in supp_list:
    if cand_supp_num == 0:
        print('The support is fixed...')
    elif cand_supp_num == 1:
        print(f'Extract the candidate catalysts containing {extract} as a support...')
        cand = cand[cand['Support'] == extract]

elif CalT_num == 1 and extract in CalT_list:
    print(f'Extraction of catalyst candidates to be calcined at {extract}°C...')
    cand = cand[cand['Cal. T.'] == extract]

else:
    print('!!!!!!!!!! >>>>>>>>>> ERROR <<<<<<<<<< !!!!!!!!!!')
    print(f'{extract} is not present in the catalyst candidates. Please check for your input errors.')
    print('extract.py is shut down...')
    sys.exit()

cand = gen_data.cand_str(cand, cand_pgm_num, cand_add_num, cand_supp_num, CalT_num, cand_wt_cols)

if K_cluster <= 50:
    plt.figure(facecolor='white', figsize = (6,12))
elif K_cluster <= 100 and K_cluster > 50:
    plt.figure(facecolor='white', figsize = (6,24))
elif K_cluster > 100:
    plt.figure(facecolor='white', figsize = (6,36))

top_catal = cand.sort_values("ei")
top_catal = top_catal.iloc[-catal_number:, :]
x = np.arange(len(top_catal))

plt.barh(x,top_catal.loc[:,"ei"], color = 'blue')
plt.yticks(x, top_catal.loc[:, 'Top catal.'], fontsize=8)
plt.title(fig_title_top)
plt.xlabel('Expected improvement (EI)')
plt.ylabel('Candidate catalyst composition (wt%)')
plt.savefig(out_fig_name_top, dpi=600, bbox_inches="tight")

if len(cand) > K_cluster:
    model = opt_ETR(n_estimators = 100, n_jobs = -1)
    model.fit(feat, target)
    k = KMeans(n_clusters=K_cluster, random_state = 1107)
    cluster = k.fit_predict(cand.loc[:, cand_wt_cols])
    cluster = pd.Series(cluster, index=cand.index, name='cluster')
    cand_clus = pd.concat([cand, cluster], axis=1)
    clus_high = cand_clus.sort_values(by=['cluster','ei']).drop_duplicates(subset=['cluster'],keep='last')
    clus_high = clus_high.sort_values(by='ei', ascending=False)
    clus_high.to_csv(out_csv_file_K)
    print('clus_high_shape:', clus_high.shape)
    clus_high = gen_data.clus_high_str(clus_high, cand_pgm_num, cand_add_num, cand_supp_num, CalT_num, cand_wt_cols)

    if K_cluster <= 50:
        plt.figure(facecolor='white', figsize = (6,12))
    elif K_cluster <= 100 and K_cluster > 50:
        plt.figure(facecolor='white', figsize = (6,24))
    elif K_cluster > 100:
        plt.figure(facecolor='white', figsize = (6,35))

    top_catal = clus_high.sort_values("ei")
    top_catal = clus_high.iloc[-K_cluster:, :]
    x = np.arange(len(top_catal))
    plt.yticks(x, top_catal.loc[:, 'Top catal.'], fontsize=8)
    plt.barh(x, top_catal.loc[:,'ei'], color = 'blue')
    plt.gca().invert_yaxis()
    plt.title(fig_title_K)
    plt.xlabel('Expected improvement (EI)')
    plt.ylabel('Candidate catalyst composition (wt%)')
    plt.savefig(out_fig_name_K, dpi=600, bbox_inches="tight")
elif len(cand) <= K_cluster:
    print(f'The number of extracted catalyst candidates is less than {K_cluster}, so clustering was not possible.')
