#!/usr/bin/env python
import numpy as np
import pandas as pd
import itertools
import scipy.sparse as sparse
import os, sys
sys.dont_write_bytecode = True
import conditions

import warnings
warnings.filterwarnings('ignore')

def data_convert(condition):
    converted = {}
    condition = conditions.calc_condition()
    date = condition['date']
    data_sheet_name = condition['data_sheet_name']
    Reaction = condition['Reaction']

    data_file_name = f'data/{date}_{Reaction}_dataset.xlsx'
    desc_file_name = 'data/Descriptors.xlsx'
    elem_desc_sheet = 'Descriptors_elem'

    data_cols = conditions.data_columns(condition)
    elem_cols, wt_cols = data_cols['elem'], data_cols['wt']
    target_name = condition['target_name']

    add_model, add_num = condition['add_model'], condition['add_num']

    desc_cols = conditions.desc_columns()
    basic_desc_cols = desc_cols['basic_desc_columns']
    noble_gas = desc_cols['noble_gas']
    drop_elems = desc_cols['drop_elems']

    desc = conditions.use_desc()
    add_use_desc = desc['add']

	### Read dataset and descriptors ###
    data = pd.read_excel(data_file_name, sheet_name=data_sheet_name)
    data.loc[:, elem_cols] = data.loc[:, elem_cols].fillna('H')
    data.loc[:, wt_cols] = data.loc[:, wt_cols].fillna(0)
    
    desc = pd.read_excel(desc_file_name, sheet_name=elem_desc_sheet, index_col='Symbol')

    add_desc = desc.drop(noble_gas, axis=0)
    add_desc = add_desc.drop(drop_elems, axis=0)
    add_desc = add_desc[basic_desc_cols].fillna(add_desc.mean())
    add_desc = add_desc[add_use_desc]
    add_aw = add_desc.loc[:, 'AW']
    add_desc = add_desc.drop('AW', axis=1)

	### Define feat and target ###
    if add_model == 0 and add_num != 0:
        feat = pd.DataFrame(index=np.arange(len(data)), columns=add_desc.index)
        for i in range(len(data)):
            for j in range(add_num):
                feat.loc[i, data.loc[i, f'Ad{j+1}']] = data.loc[i, f'Ad{j+1}_wt%']
        feat = feat.fillna(0)
        feat = feat.drop('H', axis=1)

    elif add_model == 1 and add_num != 0:
        feat_sub = pd.DataFrame(index=np.arange(len(data)), columns=add_desc.index)
        for i in range(len(data)):
            for j in range(add_num):
                feat_sub.loc[i, data.loc[i, f'Ad{j+1}']] = data.loc[i, f'Ad{j+1}_wt%']
        feat_sub = feat_sub.fillna(0)
        feat_sub = feat_sub.drop('H', axis=1)

        feat = sum([
            np.multiply(np.array(add_desc.loc[data.loc[:, f"Ad{i+1}"]]),
            np.array(data[f"Ad{i+1}_wt%"]).reshape(-1, 1)) /
            np.array(add_aw.loc[data.loc[:, f"Ad{i+1}"]]).reshape(-1, 1) for i in range(add_num)
            ])
        feat = np.hstack((feat_sub, feat))

    elif add_model == 2 and add_num != 0:
        feat = sum([
            np.multiply(np.array(add_desc.loc[data.loc[:, f"Ad{i+1}"]]),
            np.array(data[f"Ad{i+1}_wt%"]).reshape(-1, 1)) /
            np.array(add_aw.loc[data.loc[:, f"Ad{i+1}"]]).reshape(-1, 1) for i in range(add_num)
            ])

    ### Define target ###
    target = data.loc[:, target_name]
    target = target.fillna(target.mean())

    converted = {
        'data': data,
        'desc': desc,
        'add_desc': add_desc,
        'add_aw': add_aw,
        'feat': feat,
        'target': target
    }
    return converted

def make_catal_cand(condition, converted):
    catal_cand = {}
    condition = conditions.calc_condition()
    add_desc = converted['add_desc']
    add_aw = converted['add_aw']

    add_model = condition['add_model']
    cand_add_num = condition['cand_add_num']
    fix_add_num = condition['fix_add_num']
    essential_adds = condition['essential_adds']
    add_wt = condition['add_wt']

    ### Create the elem and wt combination ###
    if cand_add_num != 0:
        elem_combination = np.array(
            list(itertools.combinations(add_desc.index[1:], cand_add_num))) #[1:]: drop　’H' in the 1st row.
        
    if fix_add_num != 0:
        print('Fix_add_num:', fix_add_num)
        add_fix = essential_adds
        residue_adds = list(set(add_desc.index) - set(essential_adds))
        residue_adds.remove('H')
        desc_add_fix = add_desc.loc[add_fix]
        desc_add_residue = add_desc.loc[residue_adds]
        elem_comb_fix = np.array(list(itertools.combinations(desc_add_fix.index, fix_add_num)), dtype='object')
        elem_comb_residue = np.array(list(itertools.combinations(desc_add_residue.index, cand_add_num-fix_add_num)), dtype='object')
        elem_comb = np.array(list(itertools.product(elem_comb_fix, elem_comb_residue)))
        
        elem_combination = []
        for i in range(len(elem_comb)):
            elems = np.hstack((elem_comb[i, 0], elem_comb[i, 1]))
            elem_combination.append(elems)
        elem_combination = np.array(elem_combination)

    elem_wt = list(itertools.product(add_wt, repeat=cand_add_num))

    if cand_add_num != 0 and add_model == 0:
        add_vect = [pd.DataFrame(columns=add_desc.index[1:], index=np.arange(len(elem_combination))).fillna(0) for i in range(cand_add_num)]
        for s in list(add_desc.index[1:]):
            for i in range(0, cand_add_num):
                idx = elem_combination[:, i] == s
                add_vect[i].loc[idx, s] = 1.0
        add_vect = [np.array(arr) for arr in add_vect]
        add_vect_csr = [sparse.csr_matrix(m) for m in add_vect]
        add_com_csr = add_vect_csr # Don't use

    elif cand_add_num != 0 and add_model == 1:
        add_vect = [pd.DataFrame(columns=add_desc.index[1:], index=np.arange(len(elem_combination))).fillna(0) for i in range(cand_add_num)]
        for s in list(add_desc.index[1:]):
            for i in range(0, cand_add_num):
                idx = elem_combination[:, i] == s
                add_vect[i].loc[idx, s] = 1.0
        add_vect = [np.array(arr) for arr in add_vect]
        add_vect_csr = [sparse.csr_matrix(m) for m in add_vect]

        add_com = [np.array(add_desc.loc[elem_combination[:, i]], dtype='float') /
            np.array(add_aw.loc[elem_combination[:, i]], dtype='float').reshape(-1, 1)
            for i in range(0, cand_add_num)]
        add_com_csr = [sparse.csr_matrix(m) for m in add_com]

    elif cand_add_num != 0 and add_model == 2:
        add_com = [np.array(add_desc.loc[elem_combination[:, i]], dtype='float') /
            np.array(add_aw.loc[elem_combination[:, i]], dtype='float').reshape(-1, 1)
            for i in range(0, cand_add_num)]
        add_com_csr = [sparse.csr_matrix(m) for m in add_com]
        add_vect_csr = add_com_csr # Don't use

    catal_cand = {
        'add_com': add_com_csr,
        'add_vect': add_vect_csr,
        'catal_comb': elem_combination,
        'elem_wt': elem_wt
        }

    return catal_cand

def elem_wt(condition):
    cand_add_num = condition['cand_add_num']
    add_wt = condition['add_wt']

    elem_wt = list(itertools.product(add_wt, repeat=cand_add_num))
    num_elem_wt = len(elem_wt)

    return num_elem_wt

def cand_separation(cand_data):
    condition = conditions.calc_condition()
    cand_add_num = condition['cand_add_num']
    data_cols = conditions.data_columns(condition)
    cand_elem_cols, cand_wt_cols, cand_labels = data_cols['cand_elem'], data_cols['cand_wt'], data_cols['cand_labels']

    data = pd.read_csv(cand_data, sep=',')
    comp = data.loc[:, 'Unnamed: 0']
    ei = data.loc[:, 'ei']
    
    sep_data = []
    for i in range(len(data)):
            string = comp[i].split()
            string.append(ei[i])
            sep_data.append(string)
    
    cand_sep_data = pd.DataFrame(sep_data)
    cand_sep_data.columns = cand_labels
    print('cand_sep_data:', cand_sep_data.shape)
    
    sep_data_wt = []
    if cand_add_num != 0:
        for i in range(cand_add_num):
            cand_sep_data_wt = cand_sep_data[f'Ad{i+1}_wt%'].str.strip(',')
            sep_data_wt.append(cand_sep_data_wt)
    
    cand_wt = pd.DataFrame(sep_data_wt).transpose()
    cand_wt.columns = cand_wt_cols
    
    cand_elem = cand_sep_data.loc[:, cand_elem_cols]
    cand_ei = cand_sep_data.loc[:, 'ei']
    cand = pd.concat([cand_elem, cand_wt, cand_ei], axis=1)
    cand = cand.reindex(columns=cand_labels)
    cand[cand_wt_cols] = cand[cand_wt_cols].astype(float)
    return cand

def cand_str(cand, cand_add_num, cand_wt_cols):
    if cand[cand_wt_cols[0]].dtypes == float:
        cand[cand_wt_cols] = cand[cand_wt_cols].astype(str)
    else:
        pass
    
    if cand_add_num == 2:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2],  sep=', ')

    elif cand_add_num == 3:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3],  sep=', ')

    elif cand_add_num == 4:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3, cand_Ad4],  sep=', ')

    elif cand_add_num == 5:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        cand_Ad5 = cand['Ad5'].str.cat([cand['Ad5_wt%']], sep=' ')
        cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5],  sep=', ')

    cand[cand_wt_cols] = cand[cand_wt_cols].astype(float)
    return cand
