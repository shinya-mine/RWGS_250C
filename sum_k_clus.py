#!/usr/bin/env python
import numpy as np
import pandas as pd
import glob
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

PATH = 'results'
os.makedirs(PATH, exist_ok = True)

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--catal_number', type=int, default=50)
parser.add_argument('-k', '--K_cluster', type=int, default=50)

args = vars(parser.parse_args())

catal_number = args['catal_number']
K_cluster = args['K_cluster']

### INPUT FILES ###
condition = conditions.calc_condition()
date = condition['date']
Reaction = condition['Reaction']
ML_model = condition['ML_model']
Search_method = condition['Search_method']

data_file_name = f'data/{date}_{Reaction}_dataset.xlsx'
data_sheet_name = condition['data_sheet_name']
skip_rows = condition['skip_rows']
desc_file_name = 'data/Descriptors.xlsx'
elem_desc_sheet = 'Descriptors'
supp_desc_sheet = 'Descriptors_supp'

data_cols = conditions.data_columns(Reaction, condition)
cand_elem_cols = data_cols['cand_elem']
cand_wt_cols = data_cols['cand_wt']
cand_labels = data_cols['cand_labels']

if Reaction == 'rwgs_250' or Reaction == 'rwgs_300' or Reaction == 'CH3OH':
        target_name = condition['target_name']
elif Reaction == 'N2O' or Reaction == 'H2SCR' or Reaction == 'NH3SCR' or Reaction == 'CH4':
        target_temp = condition['target_temp']

pgm_model, add_model, supp_model = condition['pgm_model'], condition['add_model'], condition['supp_model']
cand_pgm_num, cand_add_num, cand_supp_num, CalT_num  = condition['cand_pgm_num'], condition['cand_add_num'], condition['supp_num'], condition['CalT_num']
fix_pgm_num, fix_add_num, fix_supp_num = condition['fix_pgm_num'], condition['fix_add_num'], condition['fix_supp_num']
essential_pgms, essential_adds, essential_supp = condition['essential_pgms'], condition['essential_adds'], condition['essential_supp']

max_pgm_wt = condition['max_pgm_wt']
essential_pgm_wt = condition['essential_pgm_wt']
pgm_wt, add_wt = condition['pgm_wt'], condition['add_wt']

if CalT_num == 1:
        CalT_list = condition['CalT_list']
else:
        pass

desc_cols = conditions.desc_columns(Reaction)
pgm_plus_ReAu = desc_cols['pgm_plus_ReAu']
basic_desc_cols = desc_cols['basic_desc_columns']
noble_gas = desc_cols['noble_gas']
drop_elems = desc_cols['drop_elems']

desc = conditions.use_desc()
pgm_use_desc, add_use_desc, supp_use_desc = desc['pgm'], desc['add'], desc['supp']

### Top 40 catalysts OUTPUT Files ###
out_csv_file_top = f'results/{date}_{Reaction}_cand_sum_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}.csv'
fig_title_top = f'{date}_{Reaction}_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}_top{catal_number}'
out_fig_name_top = f'results/{date}_{Reaction}_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}_top{catal_number}.png'
cand_data = f'results/{date}_{Reaction}_cand_sum_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}.csv' # == results_csv_file_top

### K-mean 40 catalysts OUTPUT Files ###
out_csv_file_K = f"results/{date}_{Reaction}_K{K_cluster}_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}.csv"
fig_title_K = f'{date}_{Reaction}_K{K_cluster}_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}'
out_fig_name_K = f"results/{date}_{Reaction}_K{K_cluster}_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}.png"

csv_files = glob.glob('cand/*.csv')
data_list = []
for file in csv_files:
        print(file)
        data_list.append(pd.read_csv(file))

cand = pd.concat(data_list, axis=0, sort=True)
print("cand_shape")
print(cand.shape)
cand.to_csv(out_csv_file_top, index=False)

converted = gen_data.data_convert(condition)
feat, target = converted['feat'], converted['target']

cand = gen_data.cand_separation(cand_data)
print('cand_separation:', cand.shape)
cand = gen_data.cand_str(cand, cand_pgm_num, cand_add_num, cand_supp_num, CalT_num, cand_wt_cols)

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


### PLOT K-TOP CATALYST CANDIDATES ###
model = opt_ETR(n_estimators = 100, random_state = 1107, n_jobs = -1)
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

K_catal = clus_high.sort_values("ei")
K_catal = clus_high.iloc[-K_cluster:, :]
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
