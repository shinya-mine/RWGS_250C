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
desc_file_name = f'data/Descriptors.xlsx'
elem_desc_sheet = 'Descriptors_elem'

data_cols = conditions.data_columns(condition)
elem_cols, wt_cols = data_cols['elem'], data_cols['wt']
data_labels = data_cols['labels']
target_name = condition['target_name']

add_model = condition['add_model']
cand_add_num = condition['cand_add_num']
fix_add_num = condition['fix_add_num']

### Define feat and target ###
converted = gen_data.data_convert(condition)
feat, target = converted['feat'], converted['target']
add_desc = converted['add_desc']
add_list = list(add_desc.index)
add_list.remove('H')

### USE Descriptors ###
desc = conditions.use_desc()
add_use_desc = desc['add']

### INPUT FILES (About Cand.) ###
cand_data = f'results/{date}_{Reaction}_cand_sum_{ML_model}_prop{add_model}_{Search_method}.csv'
cand_elem_cols, cand_wt_cols, cand_labels = data_cols['cand_elem'], data_cols['cand_wt'], data_cols['cand_labels']

### Top catalysts OUTPUT Files ###
out_csv_file_top = f'{PATH}/{date}_{Reaction}_cand_sum_{ML_model}_prop{add_model}_{Search_method}_{extract}.csv'
fig_title_top = f'{date}_{Reaction}_{ML_model}_prop{add_model}_{Search_method}_top{catal_number}_{extract}'
out_fig_name_top = f'{PATH}/{date}_{Reaction}_{ML_model}_prop{add_model}_{Search_method}_top{catal_number}_{extract}.png'

### K-mean catalysts OUTPUT Files ###
out_csv_file_K = f"{PATH}/{date}_{Reaction}_{ML_model}_prop{add_model}_{Search_method}_K{K_cluster}_{extract}.csv"
fig_title_K = f'{date}_{Reaction}_{ML_model}_prop{add_model}_{Search_method}_K{K_cluster}_{extract}'
out_fig_name_K = f"{PATH}/{date}_{Reaction}_{ML_model}_prop{add_model}_{Search_method}_K{K_cluster}_{extract}.png"

cand = gen_data.cand_separation(cand_data)

if extract in add_list:
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

else:
    print('!!!!!!!!!! >>>>>>>>>> ERROR <<<<<<<<<< !!!!!!!!!!')
    print(f'{extract} is not present in the catalyst candidates. Please check for your input errors.')
    print('extract.py is shut down...')
    sys.exit()

cand = gen_data.cand_str(cand, cand_add_num, cand_wt_cols)

top_catal = cand.sort_values("ei")
top_catal = top_catal.iloc[-catal_number:, :]
x_top = np.arange(len(top_catal))

if catal_number <= 50:
    plt.figure(facecolor='white', figsize = (6,12))
    plt.barh(x_top,top_catal.loc[:,"ei"], color='blue')
    plt.yticks(x_top, top_catal.loc[:, 'Top catal.'], fontsize=8)
    plt.title(fig_title_top)
    plt.xlabel('Expected improvement (EI)')
    plt.ylabel('Candidate catalyst composition (wt%)')
    plt.savefig(out_fig_name_top, dpi=600, bbox_inches="tight")

elif catal_number <= 100 and catal_number > 50:
    fig = plt.figure(figsize = (16,12))
    ax1 = fig.add_subplot(1,2,1)
    ax2 = fig.add_subplot(1,2,2)
    fig.subplots_adjust(bottom=0, left=0, top=1, right=1, wspace=1)
    
    top_catal1 = top_catal.iloc[0:int(1/2*catal_number),:]
    top_catal2 = top_catal.iloc[int(1/2*catal_number):catal_number,:]
    x_top1 = np.arange(int(1/2*len(top_catal)))

    fig.suptitle(fig_title_top, x=0.35, y=1.03, fontsize='x-large')
    ax1.barh(x_top1, top_catal1.loc[:,'ei'], color = 'blue')
    ax1.invert_yaxis()
    ax1.set_xlabel('Expected improvement (EI)')
    ax1.set_yticks(x_top1)
    ax1.set_ylabel('Candidate catalyst composition (wt%)')
    ax1.set_yticklabels(top_catal1.loc[:, 'Top catal.'])

    ax2.barh(x_top1, top_catal2.loc[:,'ei'], color = 'blue')
    ax2.invert_yaxis()
    ax2.set_xlabel('Expected improvement (EI)')
    ax2.set_yticks(x_top1)
    ax2.set_ylabel('Candidate catalyst composition (wt%)')
    ax2.set_yticklabels(top_catal2.loc[:, 'Top catal.'])
    plt.savefig(out_fig_name_top, dpi=600, bbox_inches="tight")

elif catal_number <= 150 and catal_number > 100:
    fig = plt.figure(figsize = (25,12))
    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)
    ax3 = fig.add_subplot(1,3,3)
    fig.subplots_adjust(bottom=0, left=0, top=1, right=1, wspace=1)
    
    top_catal1 = top_catal.iloc[0:int(1/3*catal_number),:]
    top_catal2 = top_catal.iloc[int(1/3*catal_number):int(2/3*catal_number),:]
    top_catal3 = top_catal.iloc[int(2/3*catal_number):catal_number,:]
    x_top1 = np.arange(int(1/3*len(top_catal)))

    fig.suptitle(fig_title_top, x=0.4, y=1.03, fontsize='x-large')
    ax1.barh(x_top1, top_catal1.loc[:,'ei'], color = 'blue')
    ax1.invert_yaxis()
    ax1.set_xlabel('Expected improvement (EI)')
    ax1.set_yticks(x_top1)
    ax1.set_ylabel('Candidate catalyst composition (wt%)')
    ax1.set_yticklabels(top_catal1.loc[:, 'Top catal.'])

    ax2.barh(x_top1, top_catal2.loc[:,'ei'], color = 'blue')
    ax2.invert_yaxis()
    ax2.set_xlabel('Expected improvement (EI)')
    ax2.set_yticks(x_top1)
    ax2.set_ylabel('Candidate catalyst composition (wt%)')
    ax2.set_yticklabels(top_catal2.loc[:, 'Top catal.'])

    ax3.barh(x_top1, top_catal3.loc[:,'ei'], color = 'blue')
    ax3.invert_yaxis()
    ax3.set_xlabel('Expected improvement (EI)')
    ax3.set_yticks(x_top1)
    ax3.set_ylabel('Candidate catalyst composition (wt%)')
    ax3.set_yticklabels(top_catal3.loc[:, 'Top catal.'])
    plt.savefig(out_fig_name_top, dpi=600, bbox_inches="tight")

elif catal_number > 150:
    fig = plt.figure(figsize = (25,12))
    ax1 = fig.add_subplot(1,4,1)
    ax2 = fig.add_subplot(1,4,2)
    ax3 = fig.add_subplot(1,4,3)
    ax4 = fig.add_subplot(1,4,4)
    fig.subplots_adjust(bottom=0, left=0, top=1, right=1, wspace=1)
    
    top_catal1 = top_catal.iloc[0:int(1/4*catal_number),:]
    top_catal2 = top_catal.iloc[int(1/4*catal_number):int(1/2*catal_number),:]
    top_catal3 = top_catal.iloc[int(1/2*catal_number):int(3/4*catal_number),:]
    top_catal4 = top_catal.iloc[int(3/4*catal_number):catal_number,:]
    x_top1 = np.arange(int(1/4*len(top_catal)))

    fig.suptitle(fig_title_top, x=0.43, y=1.03, fontsize='x-large')
    ax1.barh(x_top1, top_catal1.loc[:,'ei'], color = 'blue')
    ax1.invert_yaxis()
    ax1.set_xlabel('Expected improvement (EI)')
    ax1.set_yticks(x_top1)
    ax1.set_ylabel('Candidate catalyst composition (wt%)')
    ax1.set_yticklabels(top_catal1.loc[:, 'Top catal.'])

    ax2.barh(x_top1, top_catal2.loc[:,'ei'], color = 'blue')
    ax2.invert_yaxis()
    ax2.set_xlabel('Expected improvement (EI)')
    ax2.set_yticks(x_top1)
    ax2.set_ylabel('Candidate catalyst composition (wt%)')
    ax2.set_yticklabels(top_catal2.loc[:, 'Top catal.'])

    ax3.barh(x_top1, top_catal3.loc[:,'ei'], color = 'blue')
    ax3.invert_yaxis()
    ax3.set_xlabel('Expected improvement (EI)')
    ax3.set_yticks(x_top1)
    ax3.set_ylabel('Candidate catalyst composition (wt%)')
    ax3.set_yticklabels(top_catal3.loc[:, 'Top catal.'])

    ax4.barh(x_top1, top_catal3.loc[:,'ei'], color = 'blue')
    ax4.invert_yaxis()
    ax4.set_xlabel('Expected improvement (EI)')
    ax4.set_yticks(x_top1)
    ax4.set_ylabel('Candidate catalyst composition (wt%)')
    ax4.set_yticklabels(top_catal3.loc[:, 'Top catal.'])
    plt.savefig(out_fig_name_top, dpi=600, bbox_inches="tight")

if len(cand) > K_cluster:
    model = opt_ETR(n_estimators=100, random_state=1107, n_jobs=-1)
    model.fit(feat, target)
    k = KMeans(n_clusters=K_cluster, random_state=1107)
    cluster = k.fit_predict(cand.loc[:, cand_wt_cols])
    cluster = pd.Series(cluster, index=cand.index, name='cluster')

    cand_top = cand.copy()
    cand = pd.concat([cand, cluster], axis=1)
    cand = cand.sort_values(by=['cluster','ei']).drop_duplicates(subset=['cluster'],keep='last')
    cand = cand.sort_values(by='ei', ascending=False)
    cand.to_csv(out_csv_file_K)
    print(f'cand_K{K_cluster}_shape:', cand.shape)
    cand = gen_data.cand_str(cand, cand_add_num, cand_wt_cols)

    K_catal = cand.sort_values('ei')
    K_catal = cand.iloc[-K_cluster:, :]
    x_K = np.arange(len(K_catal))

    if K_cluster <= 50:
        plt.figure(facecolor='white', figsize = (6,12))
        plt.yticks(x_K, K_catal.loc[:, 'Top catal.'], fontsize=8)
        plt.barh(x_K, K_catal.loc[:,'ei'], color='blue')
        plt.gca().invert_yaxis()
        plt.title(fig_title_K)
        plt.xlabel('Expected improvement (EI)')
        plt.ylabel('Candidate catalyst composition (wt%)')
        plt.savefig(out_fig_name_K, dpi=600, bbox_inches="tight")

    elif K_cluster <= 100 and K_cluster > 50:
        fig = plt.figure(facecolor='white', figsize = (16,12))
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)
        fig.subplots_adjust(bottom=0, left=0, top=1, right=1, wspace=1)
        
        K_catal1 = K_catal.iloc[0:int(1/2*K_cluster),:]
        K_catal2 = K_catal.iloc[int(1/2*K_cluster):K_cluster,:]
        x_K1 = np.arange(int(1/2*len(K_catal)))

        fig.suptitle(fig_title_K, x=0.35, y=1.03, fontsize='x-large')
        ax1.barh(x_K1, K_catal1.loc[:,'ei'], color='blue')
        ax1.invert_yaxis()
        ax1.set_xlabel('Expected improvement (EI)')
        ax1.set_yticks(x_K1)
        ax1.set_ylabel('Candidate catalyst composition (wt%)')
        ax1.set_yticklabels(K_catal1.loc[:, 'Top catal.'])

        ax2.barh(x_K1, K_catal2.loc[:,'ei'], color='blue')
        ax2.invert_yaxis()
        ax2.set_xlabel('Expected improvement (EI)')
        ax2.set_yticks(x_K1)
        ax2.set_ylabel('Candidate catalyst composition (wt%)')
        ax2.set_yticklabels(K_catal2.loc[:, 'Top catal.'])
        plt.savefig(out_fig_name_K, dpi = 600, bbox_inches = "tight")
            
    elif K_cluster <= 150 and K_cluster > 100:
        fig = plt.figure(facecolor='white', figsize = (25,12))
        ax1 = fig.add_subplot(1,3,1)
        ax2 = fig.add_subplot(1,3,2)
        ax3 = fig.add_subplot(1,3,3)
        fig.subplots_adjust(bottom=0, left=0, top=1, right=1, wspace=1)
        
        K_catal1 = K_catal.iloc[0:int(1/3*K_cluster),:]
        K_catal2 = K_catal.iloc[int(1/3*K_cluster):int(2/3*K_cluster),:]
        K_catal3 = K_catal.iloc[int(2/3*K_cluster):K_cluster,:]
        x_K1 = np.arange(int(1/3*len(K_catal)))

        fig.suptitle(fig_title_K, x=0.4, y=1.03, fontsize='x-large')
        ax1.barh(x_K1, K_catal1.loc[:,'ei'], color='blue')
        ax1.invert_yaxis()
        ax1.set_xlabel('Expected improvement (EI)')
        ax1.set_yticks(x_K1)
        ax1.set_ylabel('Candidate catalyst composition (wt%)')
        ax1.set_yticklabels(K_catal1.loc[:, 'Top catal.'])

        ax2.barh(x_K1, K_catal2.loc[:,'ei'], color='blue')
        ax2.invert_yaxis()
        ax2.set_xlabel('Expected improvement (EI)')
        ax2.set_yticks(x_K1)
        ax2.set_ylabel('Candidate catalyst composition (wt%)')
        ax2.set_yticklabels(K_catal2.loc[:, 'Top catal.'])

        ax3.barh(x_K1, K_catal3.loc[:,'ei'], color='blue')
        ax3.invert_yaxis()
        ax3.set_xlabel('Expected improvement (EI)')
        ax3.set_yticks(x_K1)
        ax3.set_ylabel('Candidate catalyst composition (wt%)')
        ax3.set_yticklabels(K_catal3.loc[:, 'Top catal.'])
        plt.savefig(out_fig_name_K, dpi = 600, bbox_inches = "tight")

    elif K_cluster > 150:
        fig = plt.figure(facecolor='white', figsize = (32,12))
        ax1 = fig.add_subplot(1,4,1)
        ax2 = fig.add_subplot(1,4,2)
        ax3 = fig.add_subplot(1,4,3)
        ax4 = fig.add_subplot(1,4,4)
        fig.subplots_adjust(bottom=0, left=0, top=1, right=1, wspace=1)
        
        K_catal1 = K_catal.iloc[0:int(1/4*K_cluster),:]
        K_catal2 = K_catal.iloc[int(1/4*K_cluster):int(1/2*K_cluster),:]
        K_catal3 = K_catal.iloc[int(1/2*K_cluster):int(3/4*K_cluster),:]
        K_catal4 = K_catal.iloc[int(2/3*K_cluster):K_cluster,:]
        x_K1 = np.arange(int(1/4*len(K_catal)))

        fig.suptitle(fig_title_K, x=0.43, y=1.03, fontsize='x-large')
        ax1.barh(x_K1, K_catal1.loc[:,'ei'], color='blue')
        ax1.invert_yaxis()
        ax1.set_xlabel('Expected improvement (EI)')
        ax1.set_yticks(x_K1)
        ax1.set_ylabel('Candidate catalyst composition (wt%)')
        ax1.set_yticklabels(K_catal1.loc[:, 'Top catal.'])

        ax2.barh(x_K1, K_catal2.loc[:,'ei'], color='blue')
        ax2.invert_yaxis()
        ax2.set_xlabel('Expected improvement (EI)')
        ax2.set_yticks(x_K1)
        ax2.set_ylabel('Candidate catalyst composition (wt%)')
        ax2.set_yticklabels(K_catal2.loc[:, 'Top catal.'])

        ax3.barh(x_K1, K_catal3.loc[:,'ei'], color='blue')
        ax3.invert_yaxis()
        ax3.set_xlabel('Expected improvement (EI)')
        ax3.set_yticks(x_K1)
        ax3.set_ylabel('Candidate catalyst composition (wt%)')
        ax3.set_yticklabels(K_catal3.loc[:, 'Top catal.'])

        ax4.barh(x_K1, K_catal4.loc[:,'ei'], color='blue')
        ax4.invert_yaxis()
        ax4.set_xlabel('Expected improvement (EI)')
        ax4.set_yticks(x_K1)
        ax4.set_ylabel('Candidate catalyst composition (wt%)')
        ax4.set_yticklabels(K_catal4.loc[:, 'Top catal.'])
        plt.savefig(out_fig_name_K, dpi = 600, bbox_inches = "tight")

elif len(cand) <= K_cluster:
    print(f'The number of extracted catalyst candidates is less than {K_cluster}, so clustering was not possible.')
