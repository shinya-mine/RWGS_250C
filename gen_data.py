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
    skip_rows = condition['skip_rows']
    Reaction = condition['Reaction']

    data_file_name = f'data/{date}_{Reaction}_dataset.xlsx'
    desc_file_name = 'data/Descriptors.xlsx'
    elem_desc_sheet = 'Descriptors_elem'
    supp_desc_sheet = 'Descriptors_supp'

    data_cols = conditions.data_columns(Reaction, condition)
    elem_cols, wt_cols = data_cols['elem'], data_cols['wt']

    if Reaction == 'rwgs_250' or Reaction == 'rwgs_300' or Reaction == 'rwgs_250_1wt' or Reaction == 'CH3OH':
        target_name = condition['target_name']
    elif Reaction == 'N2O' or Reaction == 'H2SCR' or Reaction == 'NH3SCR' or Reaction == 'CH4':
        target_temp = condition['target_temp']

    pgm_model, add_model, supp_model = condition['pgm_model'], condition['add_model'], condition['supp_model']
    pgm_num, add_num, supp_num = condition['pgm_num'], condition['add_num'], condition['supp_num']

    desc_cols = conditions.desc_columns(Reaction)
    pgm_plus_ReAu = desc_cols['pgm_plus_ReAu']
    basic_desc_cols = desc_cols['basic_desc_columns']
    noble_gas = desc_cols['noble_gas']
    drop_elems = desc_cols['drop_elems']

    desc = conditions.use_desc()
    pgm_use_desc, add_use_desc, supp_use_desc = desc['pgm'], desc['add'], desc['supp']

	### Read dataset and descriptors ###
    print(f'pgm_model={pgm_model}', f'add_model={add_model}', f'supp_model={supp_model}',)
    if Reaction == 'rwgs_250_1wt':
        train_data_sheet_name = condition['train_data_sheet_name']
        test_data_sheet_name = condition['data_sheet_name']
        train_data = pd.read_excel(data_file_name, sheet_name=train_data_sheet_name, skiprows=skip_rows)
        test_data = pd.read_excel(data_file_name, sheet_name=test_data_sheet_name, skiprows=skip_rows)
        data = pd.concat([train_data, test_data], axis=0).reset_index(drop=True)
    else:
        data = pd.read_excel(data_file_name, sheet_name=data_sheet_name, skiprows=skip_rows)
        
    data.loc[:, elem_cols] = data.loc[:, elem_cols].fillna('H')
    data.loc[:, wt_cols] = data.loc[:, wt_cols].fillna(0)
    
    desc = pd.read_excel(desc_file_name, sheet_name=elem_desc_sheet, index_col='Symbol')
    if Reaction == 'N2O':
        pgm_desc = desc.loc[pgm_plus_ReAu].drop('Os', axis=0)
    else:
        pgm_desc = desc.loc[pgm_plus_ReAu].drop(['Re', 'Os'], axis=0)

    pgm_desc = pgm_desc[basic_desc_cols].fillna(pgm_desc.mean())
    pgm_desc = pgm_desc.rename(columns=lambda s: s+' (PGM)')
    pgm_desc = pgm_desc[pgm_use_desc]
    pgm_aw = pgm_desc.loc[:, 'AW (PGM)']
    pgm_desc = pgm_desc.drop('AW (PGM)', axis=1)

    add_desc = desc.drop(noble_gas, axis=0)
    add_desc = add_desc.drop(drop_elems, axis=0)
    add_desc = add_desc[basic_desc_cols].fillna(add_desc.mean())
    add_desc = add_desc.rename(columns=lambda s: s+' (Additive)')
    add_desc = add_desc[add_use_desc]
    add_aw = add_desc.loc[:, 'AW (Additive)']
    add_desc = add_desc.drop('AW (Additive)', axis=1)

    supp_desc = pd.read_excel(desc_file_name, sheet_name=supp_desc_sheet, index_col='Support_name')
    supp_desc = supp_desc[supp_desc[Reaction]==1]
    supp_desc = supp_desc[supp_use_desc].fillna(supp_desc.mean())
    supp_mw = supp_desc.loc[:, 'MW (Support)']
    supp_desc = supp_desc.drop('MW (Support)', axis=1)
    supp_desc = supp_desc.fillna(supp_desc.mean())
    print('pgm_desc:', len(pgm_desc.columns), pgm_desc.columns)
    print('add_desc:', len(add_desc.columns), add_desc.columns)
    print(f'pgm_num = {pgm_num}', f'add_num = {add_num}')

	### Define feat and target ###
    if pgm_model == 0 and pgm_num != 0:
        feat_pgm = pd.DataFrame(index=np.arange(len(data)), columns=pgm_desc.index)
        for i in range(len(data)):
            for j in range(pgm_num):
                feat_pgm.loc[i, data.loc[i, f'PGM{j+1}']] = data.loc[i, f'PGM{j+1}_wt%']
        feat_pgm = feat_pgm.fillna(0)
        feat_pgm = feat_pgm.drop('H', axis=1)
        print(f'pgm_model={pgm_model}', 'feat_pgm_shape:', feat_pgm.shape)

    elif pgm_model == 1 and pgm_num != 0:
        feat_sub_pgm = pd.DataFrame(index=np.arange(len(data)), columns=pgm_desc.index)
        for i in range(len(data)):
            for j in range(pgm_num):
                feat_sub_pgm.loc[i, data.loc[i, f'PGM{j+1}']] = data.loc[i, f'PGM{j+1}_wt%']
        feat_sub_pgm = feat_sub_pgm.fillna(0)
        feat_sub_pgm = feat_sub_pgm.drop('H', axis=1)

        feat_pgm = sum([
            np.multiply(np.array(pgm_desc.loc[data.loc[:, f"PGM{i+1}"]]),
            np.array(data[f"PGM{i+1}_wt%"]).reshape(-1, 1)) /
            np.array(pgm_aw.loc[data.loc[:, f"PGM{i+1}"]]).reshape(-1, 1) for i in range(pgm_num)
            ])
        feat_pgm = np.hstack((feat_sub_pgm, feat_pgm))
        print(f'pgm_model={pgm_model}', 'feat_shape_pgm:', feat_pgm.shape)

    elif pgm_model == 2 and pgm_num != 0:
        feat_pgm = sum([
            np.multiply(np.array(pgm_desc.loc[data.loc[:, f"PGM{i+1}"]]),
            np.array(data[f"PGM{i+1}_wt%"]).reshape(-1, 1)) /
            np.array(pgm_aw.loc[data.loc[:, f"PGM{i+1}"]]).reshape(-1, 1) for i in range(pgm_num)
            ])
        print(f'pgm_model={pgm_model}', 'feat_shape_pgm:', feat_pgm.shape)

    if add_model == 0 and add_num != 0:
        feat_add = pd.DataFrame(index=np.arange(len(data)), columns=add_desc.index)
        for i in range(len(data)):
            for j in range(add_num):
                feat_add.loc[i, data.loc[i, f'Ad{j+1}']] = data.loc[i, f'Ad{j+1}_wt%']
        feat_add = feat_add.fillna(0)
        feat_add = feat_add.drop('H', axis=1)
        print(f'add_model={add_model}', 'feat_add_shape:', feat_add.shape)

    elif add_model == 1 and add_num != 0:
        feat_sub_add = pd.DataFrame(index=np.arange(len(data)), columns=add_desc.index)
        for i in range(len(data)):
            for j in range(add_num):
                feat_sub_add.loc[i, data.loc[i, f'Ad{j+1}']] = data.loc[i, f'Ad{j+1}_wt%']
        feat_sub_add = feat_sub_add.fillna(0)
        feat_sub_add = feat_sub_add.drop('H', axis=1)

        feat_add = sum([
            np.multiply(np.array(add_desc.loc[data.loc[:, f"Ad{i+1}"]]),
            np.array(data[f"Ad{i+1}_wt%"]).reshape(-1, 1)) /
            np.array(add_aw.loc[data.loc[:, f"Ad{i+1}"]]).reshape(-1, 1) for i in range(add_num)
            ])
        feat_add = np.hstack((feat_sub_add, feat_add))
        print(f'add_model={add_model}', 'feat_shape_add:', feat_add.shape)

    elif add_model == 2 and add_num != 0:
        feat_add = sum([
            np.multiply(np.array(add_desc.loc[data.loc[:, f"Ad{i+1}"]]),
            np.array(data[f"Ad{i+1}_wt%"]).reshape(-1, 1)) /
            np.array(add_aw.loc[data.loc[:, f"Ad{i+1}"]]).reshape(-1, 1) for i in range(add_num)
            ])
        print(f'add_model={add_model}', 'feat_shape_add:', feat_add.shape)
    
    if pgm_num == 0:
        feat = feat_add
        print('feat_shape_Add:', feat.shape)
    elif add_num == 0:
        feat = feat_pgm
        print('feat_shape_PGM:', feat.shape)
    else:
        feat = np.hstack((feat_pgm, feat_add))
        print('feat_shape_PGM+Add:', feat.shape)
    
    supp_name_list = list(supp_desc.index)
    print('Kinds of support:', len(supp_name_list))
    print(supp_name_list)
    if len(supp_name_list) != 0:
        supp_data = supp_desc.loc[data.loc[:, 'Support_name']]
        print('Support_desc_columns:', len(supp_desc.columns), supp_desc.columns)
    else:
        pass

    if supp_num != 0:
        if supp_model == 0 and pgm_num != 0 and add_num != 0:
            feat_sub_pgm_conv = pd.DataFrame(index=np.arange(len(data)), columns=pgm_desc.index)
            for i in range(len(data)):
                for j in range(pgm_num):
                    feat_sub_pgm_conv.loc[i, data.loc[i, f'PGM{j+1}']] = data.loc[i, f'PGM{j+1}_wt%']
            feat_sub_pgm_conv = feat_sub_pgm_conv.fillna(0)
            feat_sub_pgm_conv = feat_sub_pgm_conv.drop('H', axis=1)

            feat_sub_add_conv = pd.DataFrame(index=np.arange(len(data)), columns=add_desc.index)
            for i in range(len(data)):
                for j in range(add_num):
                    feat_sub_add_conv.loc[i, data.loc[i, f'Ad{j+1}']] = data.loc[i, f'Ad{j+1}_wt%']
            feat_sub_add_conv = feat_sub_add_conv.fillna(0)
            feat_sub_add_conv = feat_sub_add_conv.drop('H', axis=1)
            feat_sub_conv = pd.concat([feat_sub_pgm_conv, feat_sub_add_conv], axis=1)

            feat_supp = pd.DataFrame(index=np.arange(len(data)), columns=supp_name_list)
            for i in range(len(data)):
                feat_supp.loc[i, data.loc[i, 'Support_name']] = 100 - feat_sub_conv.iloc[i].sum()
            feat_supp = feat_supp.fillna(0)
            print(f'supp_model={supp_model}', 'feat_shape_supp:', feat_supp.shape)
        
        elif supp_model == 0 and pgm_num == 0 and add_num != 0:
            feat_sub_conv = pd.DataFrame(index=np.arange(len(data)), columns=add_desc.index)
            for i in range(len(data)):
                for j in range(add_num):
                    feat_sub_conv.loc[i, data.loc[i, f'Ad{j+1}']] = data.loc[i, f'Ad{j+1}_wt%']
            feat_sub_conv = feat_sub_conv.fillna(0)
            feat_sub_conv = feat_sub_conv.drop('H', axis=1)
            
            feat_supp = pd.DataFrame(index=np.arange(len(data)), columns=supp_name_list)
            for i in range(len(data)):
                feat_supp.loc[i, data.loc[i, 'Support_name']] = 100 - feat_sub_conv.iloc[i].sum()
            feat_supp = feat_supp.fillna(0)
            print(f'supp_model={supp_model}', 'feat_shape_supp:', feat_supp.shape)

        elif supp_model == 1 and pgm_num != 0 and add_num != 0:
            feat_sub_pgm_prop1 = pd.DataFrame(index=np.arange(len(data)), columns=pgm_desc.index)
            for i in range(len(data)):
                for j in range(pgm_num):
                    feat_sub_pgm_prop1.loc[i, data.loc[i, f'PGM{j+1}']] = data.loc[i, f'PGM{j+1}_wt%']
            feat_sub_pgm_prop1 = feat_sub_pgm_prop1.fillna(0)
            feat_sub_pgm_prop1 = feat_sub_pgm_prop1.drop('H', axis=1)

            feat_sub_add_prop1 = pd.DataFrame(index=np.arange(len(data)), columns=add_desc.index)
            for i in range(len(data)):
                for j in range(add_num):
                    feat_sub_add_prop1.loc[i, data.loc[i, f'Ad{j+1}']] = data.loc[i, f'Ad{j+1}_wt%']
            feat_sub_add_prop1 = feat_sub_add_prop1.fillna(0)
            feat_sub_add_prop1 = feat_sub_add_prop1.drop('H', axis=1)

            feat_sub_prop1 = pd.concat([feat_sub_pgm_prop1, feat_sub_add_prop1], axis=1)

            feat_sub_supp = pd.DataFrame(index=np.arange(len(data)), columns=supp_name_list)
            for i in range(len(data)):
                feat_sub_supp.loc[i, data.loc[i, 'Support_name']] = 100 - feat_sub_prop1.iloc[i].sum()
            feat_sub_supp = feat_sub_supp.fillna(0)
            feat_supp = np.array(supp_data)
            feat_supp = np.hstack((feat_sub_supp, feat_supp))
            print(f'supp_model={supp_model}', 'feat_shape_supp:', feat_supp.shape)
        
        elif supp_model == 1 and pgm_num == 0 and add_num != 0:
            feat_sub_prop1 = pd.DataFrame(index=np.arange(len(data)), columns=add_desc.index)
            for i in range(len(data)):
                for j in range(add_num):
                    feat_sub_prop1.loc[i, data.loc[i, f'Ad{j+1}']] = data.loc[i, f'Ad{j+1}_wt%']
            feat_sub_prop1 = feat_sub_prop1.fillna(0)
            feat_sub_prop1 = feat_sub_prop1.drop('H', axis=1)

            feat_sub_supp = pd.DataFrame(index=np.arange(len(data)), columns=supp_name_list)
            for i in range(len(data)):
                feat_sub_supp.loc[i, data.loc[i, 'Support_name']] = 100 - feat_sub_prop1.iloc[i].sum()
            feat_sub_supp = feat_sub_supp.fillna(0)
            feat_supp = np.array(supp_data)
            feat_supp = np.hstack((feat_sub_supp, feat_supp))
            print(f'supp_model={supp_model}', 'feat_shape_supp:', feat_supp.shape)

        elif supp_model == 2:
            feat_supp = np.array(supp_data)
            print(f'supp_model={supp_model}', 'feat_shape_supp:', feat_supp.shape)

        feat = np.hstack((feat, feat_supp))
        print('feat_shape_PGM+Add+Supp:', feat.shape, '==> feat_shape')
    
    else:
        pass
    
    if Reaction == 'CH3OH' or Reaction == 'NH3SCR':
        feat_CalT = data.loc[:, 'Calc. temp. (℃)']
        feat_CalT = np.array(feat_CalT).reshape(-1,1)
        print('feat_shape_CalT:', feat_CalT.shape)

        feat = np.hstack((feat, feat_CalT))
        print('feat_shape_(PGM)+Add+Supp+CalT:', feat.shape, '===> feat_shape (True)')
    
    ### Define target ###
    if Reaction == 'rwgs_250' or Reaction == 'rwgs_250_1wt' or Reaction == 'rwgs_300':
        print(f'Reaction name is {Reaction}')
        target = data.loc[:, target_name]
        target = target.fillna(target.mean())
        print('target_shape:', target.shape)
    
    elif Reaction == 'CH3OH':
        print('target1: CH3OH rate (mmol g-1 h-1)')
        print('target2: CH4 rate (mmol g-1 h-1)')
        print('target3: CH3OH rate x selec')
        target_cols = data.loc[:, ["CH3OH rate (mmol g-1 h-1)", "CH4 rate (mmol g-1 h-1)", "CH3OH rate x selec"]]
        target_cols = target_cols.fillna(target_cols.mean())
        print('target_cols_shape:', target_cols.shape)
        
        target = target_cols.loc[:, target_name]
        print(f'Select {target_name} as Target')
    
    elif Reaction == 'N2O':
        print(f'target: T{target_temp} (℃)')
        target = data.loc[:, f'T{target_temp} (℃)']
        target = target.fillna(900) # sufficiently high temperature
        print('target_shape:', target.shape)
    
    elif Reaction == 'H2SCR':
        print('target:', f'N2 yield_T{target_temp}')
        target = data.loc[:, f'N2 yield_T{target_temp}']
        target = target.replace('-',250)
        target = target.astype(float)
        print('target_shape:', target.shape)
    
    elif Reaction == 'NH3SCR':
        print(f'target: T{target_temp}_N2 yield')
        target = data.loc[:, f'T{target_temp}_N2 yield']
        target = target.replace('-',400) # sufficiently high temperature
        target = target.astype(float)
        print('target_shape:', target.shape)
    
    elif Reaction == 'CH4':
        print('target:', f'T{target_temp} (℃)')
        target = data.loc[:, f'T{target_temp} (℃)']
        target = target.fillna(900) # sufficiently high temperature
        print('target_shape:', target.shape)
    
    elif Reaction == 'EtOH':
        print('Now preparation ...')
    
    converted = {
        'data': data,
        'desc': desc,
        'pgm_desc': pgm_desc,
        'add_desc': add_desc,
        'supp_desc': supp_desc,
        'pgm_aw': pgm_aw,
        'add_aw': add_aw,
        'supp_mw': supp_mw,
        'feat': feat,
        'target': target
    }
    return converted

def make_catal_cand(condition, converted):
    catal_cand = {}
    condition = conditions.calc_condition()
    Reaction = condition['Reaction']
    pgm_desc, add_desc, supp_desc = converted['pgm_desc'], converted['add_desc'], converted['supp_desc']
    pgm_aw, add_aw, supp_mw = converted['pgm_aw'], converted['add_aw'], converted['supp_mw']

    data_cols = conditions.data_columns(Reaction, condition)
    cand_elem_cols = data_cols['cand_elem']
    pgm_model, add_model, supp_model = condition['pgm_model'], condition['add_model'], condition['supp_model']
    cand_pgm_num, cand_add_num, cand_supp_num, CalT_num  = condition['cand_pgm_num'], condition['cand_add_num'], condition['supp_num'], condition['CalT_num']
    fix_pgm_num, fix_add_num, fix_supp_num = condition['fix_pgm_num'], condition['fix_add_num'], condition['fix_supp_num']
    essential_pgms, essential_adds, essential_supp = condition['essential_pgms'], condition['essential_adds'], condition['essential_supp']

    print('#' * 60)
    print('!!! 2022/03/10: modified...(by Mine): Ir droped from pgm_desc !!!')
    if pgm_model == 2:
        pgm_desc = pgm_desc.drop(['Ir'], axis=0)
    else:
        pass
    print('!!! When the Ir reagent arrives, this line needs to be removed completely !!!')
    print('#' * 60)

    # in the future ...
    specific_pgms = condition['specific_pgms']
    specific_pgm_wt = condition['specific_pgm_wt']
    specific_adds = condition['specific_adds']
    specific_add_wt = condition['specific_add_wt']

    max_pgm_wt = condition['max_pgm_wt']
    essential_pgm_wt = condition['essential_pgm_wt']
    pgm_wt, add_wt = condition['pgm_wt'], condition['add_wt']

    if CalT_num == 1:
        CalT_list = condition['CalT_list']

    ### Create the elem and wt combination ###
    if cand_pgm_num != 0:
        elem_combination_pgm = np.array(
            list(itertools.combinations(pgm_desc.index[1:], cand_pgm_num))) #[1:]: drop ’H' in the 1st row.
        print('elem_PGM_comb:', len(elem_combination_pgm))
        
    if cand_add_num != 0:
        elem_combination_add = np.array(
            list(itertools.combinations(add_desc.index[1:], cand_add_num))) #[1:]: drop　’H' in the 1st row.
        print('elem_Add_comb:', len(elem_combination_add))
    
    if fix_pgm_num != 0:
        print('Fix_pgm_num:', fix_pgm_num)
        pgm_fix = essential_pgms
        residue_pgms = list(set(pgm_desc.index) - set(essential_pgms))
        residue_pgms.remove('H')
        desc_pgm_fix = pgm_desc.loc[pgm_fix]
        desc_pgm_residue = pgm_desc.loc[residue_pgms]
        elem_comb_pgm_fix = np.array(list(itertools.combinations(desc_pgm_fix.index, fix_pgm_num)), dtype='object')
        print('elem_comb_pgm_fix_shape:', elem_comb_pgm_fix.shape)
        elem_comb_pgm_residue = np.array(list(itertools.combinations(desc_pgm_residue.index, cand_pgm_num-fix_pgm_num)), dtype='object')
        print('elem_comb_pgm_residue_shape:', elem_comb_pgm_residue.shape)
        elem_comb_pgm = np.array(list(itertools.product(elem_comb_pgm_fix, elem_comb_pgm_residue)))
        
        elem_combination_pgm = []
        for i in range(len(elem_comb_pgm)):
            elems = np.hstack((elem_comb_pgm[i, 0], elem_comb_pgm[i, 1]))
            elem_combination_pgm.append(elems)
        elem_combination_pgm = np.array(elem_combination_pgm)
        print('elem_comb_pgm_oroginal_shape:', elem_combination_pgm.shape)
        
    if fix_add_num != 0:
        print('Fix_add_num:', fix_add_num)
        add_fix = essential_adds
        residue_adds = list(set(add_desc.index) - set(essential_adds))
        residue_adds.remove('H')
        desc_add_fix = add_desc.loc[add_fix]
        desc_add_residue = add_desc.loc[residue_adds]
        elem_comb_add_fix = np.array(list(itertools.combinations(desc_add_fix.index, fix_add_num)), dtype='object')
        print('elem_comb_add_fix_shape:', elem_comb_add_fix.shape)
        elem_comb_add_residue = np.array(list(itertools.combinations(desc_add_residue.index, cand_add_num-fix_add_num)), dtype='object')
        print('elem_comb_add_residue_shape:', elem_comb_add_residue.shape)
        elem_comb_add = np.array(list(itertools.product(elem_comb_add_fix, elem_comb_add_residue)))
        
        elem_combination_add = []
        for i in range(len(elem_comb_add)):
            elems = np.hstack((elem_comb_add[i, 0], elem_comb_add[i, 1]))
            elem_combination_add.append(elems)
        elem_combination_add = np.array(elem_combination_add)
        print('elem_comb_add_shape:', elem_combination_add.shape)

    if cand_pgm_num != 0 and cand_add_num != 0:
        # Create a list of elements with one row and 'pgm_num+add_num' columns.
        elem_combination = np.array(list(itertools.product(elem_combination_pgm, elem_combination_add)))
        elem_combination = np.array(list(itertools.chain.from_iterable(list(itertools.chain.from_iterable(elem_combination)))), dtype='object').reshape(
            len(elem_combination_pgm)*len(elem_combination_add), -1)
        print('elem_combs_(PGMs+Adds):', len(elem_combination))
    elif cand_pgm_num == 0 and cand_add_num != 0:
        elem_combination = elem_combination_add
        print('elem_combs_(Adds):', len(elem_combination))
    elif cand_pgm_num != 0 and cand_add_num == 0:
        elem_combination = elem_combination_pgm
        print('elem_combs_(PGMs):', len(elem_combination))
    
    if cand_supp_num != 0:
        supp_name_list = list(supp_desc.index)
        if fix_supp_num != 0:
            supp_name_list = essential_supp
        print('essential_supp_list:', supp_name_list)

        elem_supp_combination = []
        for i in range(len(elem_combination)):
            for j in range(len(supp_name_list)):
                elem_supp = np.append(elem_combination[i], supp_name_list[j])
                elem_supp_combination.append(elem_supp)
        elem_supp_combination = np.array(elem_supp_combination)
        print('elem+supp_combs_(PGMs+Adds+Supp):', len(elem_supp_combination))
    else:
        elem_supp_combination = elem_combination
        print('elem_(+supp)_combs_(PGMs+Adds):', len(elem_supp_combination))

    if CalT_num == 1:
        Catal_combination = []
        for i in range(len(elem_supp_combination)):
            for j in range(len(CalT_list)):
                CalT = np.append(elem_supp_combination[i], CalT_list[j])
                Catal_combination.append(CalT)
        Catal_combination = np.array(Catal_combination)
        print('Total number of combinations of catalyst compositions:',len(Catal_combination))

    else:
        Catal_combination = elem_supp_combination
        print('Total number of combinations of catalyst compositions:',len(Catal_combination))

    print('number of pgm_wt:', len(pgm_wt))
    print('number of add_wt:', len(add_wt))
    print('search_pgm_num =', cand_pgm_num, 'search_add_num =', cand_add_num)

    ## itertools.product: Directly multiply the number of 'pgm_wt' datapoints by 'pgm_num' times. ##
    #1 Create wt% combinations of the six PGMs except essential_pgm, and convert them to numpy format.
    elem_wt_pgm = list(itertools.product(pgm_wt, repeat=cand_pgm_num-fix_pgm_num))
    elem_wt_pgm_np = np.array(elem_wt_pgm).reshape(len(elem_wt_pgm),-1)
    print('elem_wt_pgm_(without_essential_pgms)_shape:', len(elem_wt_pgm))
    #2 Create a combination of essential_pgm wt% and convert it to numpy format.
    if fix_pgm_num != 0:
        elem_wt_fix_pgm = list(itertools.product(essential_pgm_wt, repeat=fix_pgm_num))
        elem_wt_fix_pgm_np = np.array(elem_wt_fix_pgm).reshape(len(elem_wt_fix_pgm),-1)
        print('elem_wt_fix_pgm_shape:', len(elem_wt_fix_pgm))
    #3 Create Additive wt% combinations and convert to numpy format.
    elem_wt_add = list(itertools.product(add_wt, repeat=cand_add_num))
    elem_wt_add_np = np.array(elem_wt_add).reshape(len(elem_wt_add), -1)
    print('elem_wt_add_shape:', len(elem_wt_add))
    # Smooth and reshape elem_wt by acting on itertools.chain.from_iterable twice.
    if fix_pgm_num != 0:
        print('elem_wt_pgm * elem_wt_fix_pgm * elem_wt_add =', len(elem_wt_pgm)*len(elem_wt_fix_pgm)*len(elem_wt_add))
        elem_wt = np.array(list(itertools.product(elem_wt_fix_pgm_np, elem_wt_pgm_np, elem_wt_add_np)))
        elem_wt = np.array(list(itertools.chain.from_iterable(
            list(itertools.chain.from_iterable(elem_wt))))).reshape(
                len(elem_wt_pgm)*len(elem_wt_fix_pgm)*len(elem_wt_add), -1)
        print('elem_wt_PGM+Add_shape:', elem_wt.shape)
    else:
        print('elem_wt_pgm * elem_wt_add =', len(elem_wt_pgm)*len(elem_wt_add))
        elem_wt = np.array(list(itertools.product(elem_wt_pgm_np, elem_wt_add_np)))
        elem_wt = np.array(list(itertools.chain.from_iterable(
            list(itertools.chain.from_iterable(elem_wt))))).reshape(len(elem_wt_pgm)*len(elem_wt_add), -1)
        print('elem_wt_PGM+Add_shape:', elem_wt.shape)

    # (Essential_PGM_wt) + PGM_wt <= max_pgm_wt
    if 'Support' in cand_elem_cols:
        cand_elem_cols.remove('Support')
    else:
        pass
    elem_wt = pd.DataFrame(elem_wt, columns=cand_elem_cols)
    elem_wt_pgm_sum = elem_wt.loc[:, cand_elem_cols[:cand_pgm_num]].sum(axis=1)
    elem_wt_pgm_sum = pd.DataFrame(elem_wt_pgm_sum, columns=['PGM_sum'])
    elem_wt_sum = pd.concat([elem_wt, elem_wt_pgm_sum], axis=1)
    idx_PGM = elem_wt_sum.loc[:, 'PGM_sum'] <= max_pgm_wt
    elem_wt = elem_wt_sum[idx_PGM]
    elem_wt = elem_wt.drop('PGM_sum', axis=1)
    elem_wt = np.array(elem_wt)
    print(f'elem_wt_PGM+Add(PGM_wt <= {max_pgm_wt})_shape:', elem_wt.shape)
    
    if cand_supp_num != 0:
        elem_wt_total = []
        elem_wt_sum = np.sum(elem_wt, axis=1)
        elem_wt_sup = np.array(100-elem_wt_sum)
        print('elem_sup_shape:', elem_wt_sup.shape)
        for i in range(len(elem_wt_sum)):
            elem_wt_supp = np.append(elem_wt[i], 100-elem_wt_sum[i])
            elem_wt_total.append(elem_wt_supp)
        elem_wt_total = np.array(elem_wt_total)
        print('elem_wt_total_number (PGM+Add+Supp):', len(elem_wt_total))
        print('elem_wt_shape:', elem_wt_total.shape)
    else:
        elem_wt_total = elem_wt
        print('elem_wt_total_number (PGM+Add):', len(elem_wt_total))
        print('elem_wt_shape:', elem_wt_total.shape)
    
    if cand_pgm_num != 0 and pgm_model == 0:
        pgm_vect = [pd.DataFrame(columns=pgm_desc.index[1:], index=np.arange(len(Catal_combination))).fillna(0) for i in range(cand_pgm_num)]
        for s in list(pgm_desc.index[1:]):
            for i in range(cand_pgm_num):
                idx = Catal_combination[:, i] == s
                pgm_vect[i].loc[idx, s] = 1.0
        pgm_vect = [np.array(arr) for arr in pgm_vect]
        pgm_vect_csr = [sparse.csr_matrix(m) for m in pgm_vect]
        pgm_com_csr = pgm_vect_csr # Don't use
        print(f'pgm_model={pgm_model}', 'pgm_vect_shape:', pgm_vect[0].shape)

    elif cand_pgm_num != 0 and pgm_model == 1:
        pgm_vect = [pd.DataFrame(columns=pgm_desc.index[1:], index=np.arange(len(Catal_combination))).fillna(0) for i in range(cand_pgm_num)]
        for s in list(pgm_desc.index[1:]):
            for i in range(cand_pgm_num):
                idx = Catal_combination[:, i] == s
                pgm_vect[i].loc[idx, s] = 1.0
        pgm_vect = [np.array(arr) for arr in pgm_vect]
        pgm_vect_csr = [sparse.csr_matrix(m) for m in pgm_vect]
        print(f'pgm_model={pgm_model}', 'pgm_vect_shape:', pgm_vect[0].shape)

        pgm_com = [np.array(pgm_desc.loc[Catal_combination[:, i]], dtype='float') /
            np.array(pgm_aw.loc[Catal_combination[:, i]], dtype='float').reshape(-1, 1) for i in range(cand_pgm_num)]
        pgm_com_csr = [sparse.csr_matrix(m) for m in pgm_com]
        print(f'pgm_model={pgm_model}', 'pgm_com_shape:', pgm_com[0].shape)

    elif cand_pgm_num != 0 and pgm_model == 2:
        pgm_com = [np.array(pgm_desc.loc[Catal_combination[:, i]], dtype='float') /
            np.array(pgm_aw.loc[Catal_combination[:, i]], dtype='float').reshape(-1, 1) for i in range(cand_pgm_num)]
        pgm_com_csr = [sparse.csr_matrix(m) for m in pgm_com]
        pgm_vect_csr = pgm_com_csr # Don't use
        print(f'pgm_model={pgm_model}', 'pgm_com_shape:', pgm_com[0].shape)

    if cand_add_num != 0 and add_model == 0:
        add_vect = [pd.DataFrame(columns=add_desc.index[1:], index=np.arange(len(Catal_combination))).fillna(0) for i in range(cand_add_num)]
        for s in list(add_desc.index[1:]):
            for i in range(cand_pgm_num, cand_pgm_num+cand_add_num):
                idx = Catal_combination[:, i] == s
                add_vect[i-cand_pgm_num].loc[idx, s] = 1.0
        add_vect = [np.array(arr) for arr in add_vect]
        add_vect_csr = [sparse.csr_matrix(m) for m in add_vect]
        add_com_csr = add_vect_csr # Don't use
        print(f'add_model={add_model}', 'add_vect_shape:', add_vect[0].shape)

    elif cand_add_num != 0 and add_model == 1:
        add_vect = [pd.DataFrame(columns=add_desc.index[1:], index=np.arange(len(Catal_combination))).fillna(0) for i in range(cand_add_num)]
        for s in list(add_desc.index[1:]):
            for i in range(cand_pgm_num, cand_pgm_num+cand_add_num):
                idx = Catal_combination[:, i] == s
                add_vect[i-cand_pgm_num].loc[idx, s] = 1.0
        add_vect = [np.array(arr) for arr in add_vect]
        add_vect_csr = [sparse.csr_matrix(m) for m in add_vect]
        print(f'add_model={add_model}', 'add_vect_shape:', add_vect[0].shape)

        add_com = [np.array(add_desc.loc[Catal_combination[:, i]], dtype='float') /
            np.array(add_aw.loc[Catal_combination[:, i]], dtype='float').reshape(-1, 1)
            for i in range(cand_pgm_num, cand_pgm_num+cand_add_num)]
        add_com_csr = [sparse.csr_matrix(m) for m in add_com]
        print(f'add_model={add_model}', 'add_com_shape:', add_com[0].shape)

    elif cand_add_num != 0 and add_model == 2:
        add_com = [np.array(add_desc.loc[Catal_combination[:, i]], dtype='float') /
            np.array(add_aw.loc[Catal_combination[:, i]], dtype='float').reshape(-1, 1)
            for i in range(cand_pgm_num, cand_pgm_num+cand_add_num)]
        add_com_csr = [sparse.csr_matrix(m) for m in add_com]
        add_vect_csr = add_com_csr # Don't use
        print(f'add_model={add_model}', 'add_com_shape:', add_com[0].shape)

    if cand_supp_num != 0 and supp_model == 0:
        supp_vect = [pd.DataFrame(columns=supp_desc.index, index=np.arange(len(Catal_combination))).fillna(0) for i in range(cand_supp_num)]
        for s in list(supp_desc.index):
            for i in range(cand_pgm_num+cand_add_num, cand_pgm_num+cand_add_num+cand_supp_num):
                idx = Catal_combination[:, i] == s
                supp_vect[i-(cand_pgm_num+cand_add_num)].loc[idx, s] = 1.0
        supp_vect = [np.array(arr) for arr in supp_vect]
        supp_vect_csr = [sparse.csr_matrix(m) for m in supp_vect]
        supp_com_csr = supp_vect_csr # Don't use
        print(f'supp_model={supp_model}', 'supp_vect_shape:', supp_vect[0].shape)

    elif cand_supp_num != 0 and supp_model == 1:
        supp_vect = [pd.DataFrame(columns=supp_desc.index, index=np.arange(len(Catal_combination))).fillna(0) for i in range(cand_supp_num)]
        for s in list(supp_desc.index):
            for i in range(cand_pgm_num+cand_add_num, cand_pgm_num+cand_add_num+cand_supp_num):
                idx = Catal_combination[:, i] == s
                supp_vect[i-(cand_pgm_num+cand_add_num)].loc[idx, s] = 1.0
        supp_vect = [np.array(arr) for arr in supp_vect]
        supp_vect_csr = [sparse.csr_matrix(m) for m in supp_vect]
        print(f'supp_model={supp_model}', 'supp_vect_shape:', supp_vect[0].shape)
        if CalT_num == 0:
            supp_com = [np.array(supp_desc.loc[Catal_combination[:,-1]], dtype='float') /
                np.array(supp_mw.loc[Catal_combination[:,-1]], dtype='float').reshape(-1, 1)]
            supp_com_csr = [sparse.csr_matrix(m) for m in supp_com]
            print(f'supp_model={supp_model}', 'supp_com_shape:', supp_com[0].shape)
        elif CalT_num == 1:
            supp_com = [np.array(supp_desc.loc[Catal_combination[:,-2]], dtype='float') /
                np.array(supp_mw.loc[Catal_combination[:,-2]], dtype='float').reshape(-1, 1)]
            supp_com_csr = [sparse.csr_matrix(m) for m in supp_com]
            print(f'supp_model={supp_model}', 'supp_com_shape:', supp_com[0].shape)

    elif cand_supp_num != 0 and supp_model == 2:
        if CalT_num == 0:
            supp_com = [np.array(supp_desc.loc[Catal_combination[:,-1]], dtype='float') /
                np.array(supp_mw.loc[Catal_combination[:,-1]], dtype='float').reshape(-1, 1)]
            supp_com_csr = [sparse.csr_matrix(m) for m in supp_com]
            supp_vect_csr = supp_com_csr # Don't use
            print(f'supp_model={supp_model}', 'supp_com_shape:', supp_com[0].shape)
        elif CalT_num == 1:
            supp_com = [np.array(supp_desc.loc[Catal_combination[:,-2]], dtype='float') /
                np.array(supp_mw.loc[Catal_combination[:,-2]], dtype='float').reshape(-1, 1)]
            supp_com_csr = [sparse.csr_matrix(m) for m in supp_com]
            supp_vect_csr = supp_com_csr # Don't use
            print(f'supp_model={supp_model}', 'supp_com_shape:', supp_com[0].shape)

    if CalT_num == 1:
        CalT_com = [np.array(Catal_combination[:, -1], dtype='float').reshape(-1, 1)]
        CalT_com_csr = [sparse.csr_matrix(m) for m in CalT_com]
        print('CalT_com_shape:', CalT_com[0].shape)
    else:
        CalT_com = [] # Don't use
        CalT_com_csr = [] # Don't use
    
    if cand_pgm_num == 0 and cand_supp_num != 0:
        catal_cand = {
            'add_com': add_com_csr,
            'supp_com': supp_com_csr,
            'CalT_com': CalT_com_csr,
            'add_vect': add_vect_csr,
            'supp_vect': supp_vect_csr,
            'catal_comb': Catal_combination,
            'elem_wt': elem_wt_total
            }
        
    elif cand_add_num == 0 and cand_supp_num != 0:
        catal_cand = {
            'pgm_com': pgm_com_csr,
            'supp_com': supp_com_csr,
            'CalT_com': CalT_com_csr,
            'pgm_vect': pgm_vect_csr,
            'supp_vect': supp_vect_csr,
            'catal_comb': Catal_combination,
            'elem_wt': elem_wt_total
            }
        
    elif cand_pgm_num != 0 and cand_add_num != 0 and cand_supp_num == 0:
        catal_cand = {
            'pgm_com': pgm_com_csr,
            'add_com': add_com_csr,
            'CalT_com': CalT_com_csr,
            'pgm_vect': pgm_vect_csr,
            'add_vect': add_vect_csr,
            'catal_comb': Catal_combination,
            'elem_wt': elem_wt_total
            }
        
    elif cand_pgm_num == 0 and cand_supp_num == 0:
        catal_cand = {
            'add_com': add_com_csr,
            'CalT_com': CalT_com_csr,
            'add_vect': add_vect_csr,
            'catal_comb': Catal_combination,
            'elem_wt': elem_wt_total
            }
        
    elif cand_add_num == 0 and cand_supp_num == 0:
        catal_cand = {
            'pgm_com': pgm_com_csr,
            'CalT_com': CalT_com_csr,
            'pgm_vect': pgm_vect_csr,
            'catal_comb': Catal_combination,
            'elem_wt': elem_wt_total
            }
        
    else:
        catal_cand = {
            'pgm_com': pgm_com_csr,
            'add_com': add_com_csr,
            'supp_com': supp_com_csr,
            'CalT_com': CalT_com_csr,
            'pgm_vect': pgm_vect_csr,
            'add_vect': add_vect_csr,
            'supp_vect': supp_vect_csr,
            'catal_comb': Catal_combination,
            'elem_wt': elem_wt_total
            }

    return catal_cand

def elem_wt(condition, data_cols):
    cand_elem_cols = data_cols['cand_elem']
    fix_pgm_num = condition['fix_pgm_num']
    fix_add_num = condition['fix_add_num']
    cand_pgm_num, cand_add_num = condition['cand_pgm_num'], condition['cand_add_num']
    max_pgm_wt = condition['max_pgm_wt']
    essential_pgm_wt = condition['essential_pgm_wt']
    pgm_wt, add_wt = condition['pgm_wt'], condition['add_wt']
    # in the future ...
    specific_pgms = condition['specific_pgms']
    specific_pgm_wt = condition['specific_pgm_wt']
    specific_adds = condition['specific_adds']
    specific_add_wt = condition['specific_add_wt']

    ## itertools.product: Directly multiply the number of 'pgm_wt' datapoints by 'pgm_num' times. ##
    #1 Create wt% combinations of the six PGMs except essential_pgm, and convert them to numpy format.
    elem_wt_pgm = list(itertools.product(pgm_wt, repeat=cand_pgm_num-fix_pgm_num))
    elem_wt_pgm_np = np.array(elem_wt_pgm).reshape(len(elem_wt_pgm),-1)
    #2 Create a combination of essential_pgm wt% and convert it to numpy format.
    if fix_pgm_num != 0:
        elem_wt_fix_pgm = list(itertools.product(essential_pgm_wt, repeat=fix_pgm_num))
        elem_wt_fix_pgm_np = np.array(elem_wt_fix_pgm).reshape(len(elem_wt_fix_pgm),-1)
    #3 Create Additive wt% combinations and convert to numpy format.
    elem_wt_add = list(itertools.product(add_wt, repeat=cand_add_num))
    elem_wt_add_np = np.array(elem_wt_add).reshape(len(elem_wt_add), -1)
    # Smooth and reshape elem_wt by acting on itertools.chain.from_iterable twice.
    if fix_pgm_num != 0:
        elem_wt = np.array(list(itertools.product(elem_wt_fix_pgm_np, elem_wt_pgm_np, elem_wt_add_np)))
        elem_wt = np.array(list(itertools.chain.from_iterable(
            list(itertools.chain.from_iterable(elem_wt))))).reshape(
                len(elem_wt_pgm)*len(elem_wt_fix_pgm)*len(elem_wt_add), -1)
    else:
        elem_wt = np.array(list(itertools.product(elem_wt_pgm_np, elem_wt_add_np)))
        elem_wt = np.array(list(itertools.chain.from_iterable(
            list(itertools.chain.from_iterable(elem_wt))))).reshape(len(elem_wt_pgm)*len(elem_wt_add), -1)

    # (Essential_PGM_wt) + PGM_wt <= max_pgm_wt
    if 'Support' in cand_elem_cols:
        cand_elem_cols.remove('Support')
    else:
        pass
    elem_wt = pd.DataFrame(elem_wt, columns=cand_elem_cols)
    elem_wt_pgm_sum = elem_wt.loc[:, cand_elem_cols[:cand_pgm_num]].sum(axis=1)
    elem_wt_pgm_sum = pd.DataFrame(elem_wt_pgm_sum, columns=['PGM_sum'])
    elem_wt_sum = pd.concat([elem_wt, elem_wt_pgm_sum], axis=1)
    idx_PGM = elem_wt_sum.loc[:, 'PGM_sum'] <= max_pgm_wt
    elem_wt = elem_wt_sum[idx_PGM]
    elem_wt = elem_wt.drop('PGM_sum', axis=1)
    elem_wt = np.array(elem_wt)
    num_elem_wt = len(elem_wt)

    return num_elem_wt

def cand_separation(cand_data):
    condition = conditions.calc_condition()
    Reaction = condition['Reaction']
    cand_pgm_num, cand_add_num, cand_supp_num, CalT_num  = condition['cand_pgm_num'], condition['cand_add_num'], condition['supp_num'], condition['CalT_num']
    data_cols = conditions.data_columns(Reaction, condition)
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
    if cand_pgm_num != 0:
            for i in range(cand_pgm_num):
                    cand_sep_data_wt_pgm = cand_sep_data[f'PGM{i+1}_wt%'].str.strip(',')
                    sep_data_wt.append(cand_sep_data_wt_pgm)
    if cand_add_num != 0:
            for i in range(cand_add_num):
                    cand_sep_data_wt_add = cand_sep_data[f'Ad{i+1}_wt%'].str.strip(',')
                    sep_data_wt.append(cand_sep_data_wt_add)
    
    cand_wt = pd.DataFrame(sep_data_wt).transpose()
    cand_wt.columns = cand_wt_cols
    
    cand_elem = cand_sep_data.loc[:, cand_elem_cols]
    cand_ei = cand_sep_data.loc[:, 'ei']
    if CalT_num != 0:
            cand_CalT = cand_sep_data.loc[:, 'Calc. temp. (℃)']
            cand = pd.concat([cand_elem, cand_wt, cand_CalT, cand_ei], axis=1)
    elif CalT_num == 0:
            cand = pd.concat([cand_elem, cand_wt, cand_ei], axis=1)
    cand = cand.reindex(columns=cand_labels)
    cand[cand_wt_cols] = cand[cand_wt_cols].astype(float)
    return cand

def cand_str(cand, cand_pgm_num, cand_add_num, cand_supp_num, CalT_num, cand_wt_cols):
    if cand[cand_wt_cols[0]].dtypes == float:
        cand[cand_wt_cols] = cand[cand_wt_cols].astype(str)
    else:
        pass
    
    if cand_pgm_num == 0 and cand_add_num == 2 and cand_supp_num == 0:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 3 and cand_supp_num == 0:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 4 and cand_supp_num == 0:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3, cand_Ad4],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 5 and cand_supp_num == 0:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        cand_Ad5 = cand['Ad5'].str.cat([cand['Ad5_wt%']], sep=' ')
        cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 6 and cand_supp_num == 0:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        cand_Ad5 = cand['Ad5'].str.cat([cand['Ad5_wt%']], sep=' ')
        cand_Ad6 = cand['Ad6'].str.cat([cand['Ad6_wt%']], sep=' ')
        cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5, cand_Ad6],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 2 and cand_supp_num == 1:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 3 and cand_supp_num == 1:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 4 and cand_supp_num == 1:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3, cand_Ad4, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3, cand_Ad4, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 5 and cand_supp_num == 1:
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        cand_Ad5 = cand['Ad5'].str.cat([cand['Ad5_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_Ad1.str.cat([cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_Ad1, cand_Ad2, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_Ad1, cand_Ad2, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_Ad1, cand_Ad2, cand_Ad3, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_Ad1, cand_Ad2, cand_Ad3, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 1 and cand_add_num == 4 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 1 and cand_add_num == 5 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        cand_Ad5 = cand['Ad5'].str.cat([cand['Ad5_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_PGM2 = cand['PGM2'].str.cat([cand['PGM2_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_Ad1, cand_Ad2, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_Ad1, cand_Ad2, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_PGM2 = cand['PGM2'].str.cat([cand['PGM2_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_Ad1, cand_Ad2, cand_Ad3, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_Ad1, cand_Ad2, cand_Ad3, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 2 and cand_add_num == 4 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_PGM2 = cand['PGM2'].str.cat([cand['PGM2_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 2 and cand_add_num == 5 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_PGM2 = cand['PGM2'].str.cat([cand['PGM2_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        cand_Ad5 = cand['Ad5'].str.cat([cand['Ad5_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_PGM2 = cand['PGM2'].str.cat([cand['PGM2_wt%']], sep=' ')
        cand_PGM3 = cand['PGM3'].str.cat([cand['PGM3_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_PGM3, cand_Ad1, cand_Ad2, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_PGM3, cand_Ad1, cand_Ad2, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_PGM2 = cand['PGM2'].str.cat([cand['PGM2_wt%']], sep=' ')
        cand_PGM3 = cand['PGM3'].str.cat([cand['PGM3_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_PGM3, cand_Ad1, cand_Ad2, cand_Ad3, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_PGM3, cand_Ad1, cand_Ad2, cand_Ad3, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 3 and cand_add_num == 4 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_PGM2 = cand['PGM2'].str.cat([cand['PGM2_wt%']], sep=' ')
        cand_PGM3 = cand['PGM3'].str.cat([cand['PGM3_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_PGM3, cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_PGM3, cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_supp, cand_CalT],  sep=', ')

    elif cand_pgm_num == 3 and cand_add_num == 5 and cand_supp_num == 1:
        cand_PGM1 = cand['PGM1'].str.cat([cand['PGM1_wt%']], sep=' ')
        cand_PGM2 = cand['PGM2'].str.cat([cand['PGM2_wt%']], sep=' ')
        cand_PGM3 = cand['PGM3'].str.cat([cand['PGM3_wt%']], sep=' ')
        cand_Ad1 = cand['Ad1'].str.cat([cand['Ad1_wt%']], sep=' ')
        cand_Ad2 = cand['Ad2'].str.cat([cand['Ad2_wt%']], sep=' ')
        cand_Ad3 = cand['Ad3'].str.cat([cand['Ad3_wt%']], sep=' ')
        cand_Ad4 = cand['Ad4'].str.cat([cand['Ad4_wt%']], sep=' ')
        cand_Ad5 = cand['Ad5'].str.cat([cand['Ad5_wt%']], sep=' ')
        if CalT_num == 0:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_PGM3, cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5, cand_supp],  sep=', ')
        elif CalT_num == 1:
                cand_supp = cand.loc[:, 'Support'].astype(str)
                cand_supp = cand_supp.str.rstrip(',')
                cand_CalT = cand['Calc. temp. (℃)']
                cand['Top catal.'] = cand_PGM1.str.cat([cand_PGM2, cand_PGM3, cand_Ad1, cand_Ad2, cand_Ad3, cand_Ad4, cand_Ad5, cand_supp, cand_CalT],  sep=', ')
    cand[cand_wt_cols] = cand[cand_wt_cols].astype(float)
    return cand

def clus_high_str(clus_high, cand_pgm_num, cand_add_num, cand_supp_num, CalT_num, cand_wt_cols):
    if clus_high[cand_wt_cols[0]].dtypes == float:
        clus_high[cand_wt_cols] = clus_high[cand_wt_cols].astype(str)
    else:
        pass

    if cand_pgm_num == 0 and cand_add_num == 2 and cand_supp_num == 0:
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 3 and cand_supp_num == 0:
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_Ad3],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 4 and cand_supp_num == 0:
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high_Ad4 = clus_high['Ad4'].str.cat([clus_high['Ad4_wt%']], sep=' ')
        clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_Ad3, clus_high_Ad4],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 5 and cand_supp_num == 0:
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high_Ad4 = clus_high['Ad4'].str.cat([clus_high['Ad4_wt%']], sep=' ')
        clus_high_Ad5 = clus_high['Ad5'].str.cat([clus_high['Ad5_wt%']], sep=' ')
        clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_Ad5],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 6 and cand_supp_num == 0:
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high_Ad4 = clus_high['Ad4'].str.cat([clus_high['Ad4_wt%']], sep=' ')
        clus_high_Ad5 = clus_high['Ad5'].str.cat([clus_high['Ad5_wt%']], sep=' ')
        clus_high_Ad6 = clus_high['Ad6'].str.cat([clus_high['Ad6_wt%']], sep=' ')
        clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_Ad5, clus_high_Ad6],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 2 and cand_supp_num == 1:
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 3 and cand_supp_num == 1:
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_Ad3, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_Ad3, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 4 and cand_supp_num == 1:
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high_Ad4 = clus_high['Ad4'].str.cat([clus_high['Ad4_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 0 and cand_add_num == 5 and cand_supp_num == 1:
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high_Ad4 = clus_high['Ad4'].str.cat([clus_high['Ad4_wt%']], sep=' ')
        clus_high_Ad5 = clus_high['Ad5'].str.cat([clus_high['Ad5_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_Ad5, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_Ad1.str.cat([clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_Ad5, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_Ad1, clus_high_Ad2, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_Ad1, clus_high_Ad2, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 1 and cand_add_num == 4 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high_Ad4 = clus_high['Ad4'].str.cat([clus_high['Ad4_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 1 and cand_add_num == 5 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high_Ad4 = clus_high['Ad4'].str.cat([clus_high['Ad4_wt%']], sep=' ')
        clus_high_Ad5 = clus_high['Ad5'].str.cat([clus_high['Ad5_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_Ad5, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_Ad5, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_PGM2 = clus_high['PGM2'].str.cat([clus_high['PGM2_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_Ad1, clus_high_Ad2, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_Ad1, clus_high_Ad2, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_PGM2 = clus_high['PGM2'].str.cat([clus_high['PGM2_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 2 and cand_add_num == 4 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_PGM2 = clus_high['PGM2'].str.cat([clus_high['PGM2_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high_Ad4 = clus_high['Ad4'].str.cat([clus_high['Ad4_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 2 and cand_add_num == 5 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_PGM2 = clus_high['PGM2'].str.cat([clus_high['PGM2_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high_Ad4 = clus_high['Ad4'].str.cat([clus_high['Ad4_wt%']], sep=' ')
        clus_high_Ad5 = clus_high['Ad5'].str.cat([clus_high['Ad5_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_Ad5, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_Ad5, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_PGM2 = clus_high['PGM2'].str.cat([clus_high['PGM2_wt%']], sep=' ')
        clus_high_PGM3 = clus_high['PGM3'].str.cat([clus_high['PGM3_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_PGM3, clus_high_Ad1, clus_high_Ad2, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_PGM3, clus_high_Ad1, clus_high_Ad2, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_PGM2 = clus_high['PGM2'].str.cat([clus_high['PGM2_wt%']], sep=' ')
        clus_high_PGM3 = clus_high['PGM3'].str.cat([clus_high['PGM3_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_PGM3, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_PGM3, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 3 and cand_add_num == 4 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_PGM2 = clus_high['PGM2'].str.cat([clus_high['PGM2_wt%']], sep=' ')
        clus_high_PGM3 = clus_high['PGM3'].str.cat([clus_high['PGM3_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high_Ad4 = clus_high['Ad4'].str.cat([clus_high['Ad4_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_PGM3, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_PGM3, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_supp, clus_high_CalT],  sep=', ')

    elif cand_pgm_num == 3 and cand_add_num == 5 and cand_supp_num == 1:
        clus_high_PGM1 = clus_high['PGM1'].str.cat([clus_high['PGM1_wt%']], sep=' ')
        clus_high_PGM2 = clus_high['PGM2'].str.cat([clus_high['PGM2_wt%']], sep=' ')
        clus_high_PGM3 = clus_high['PGM3'].str.cat([clus_high['PGM3_wt%']], sep=' ')
        clus_high_Ad1 = clus_high['Ad1'].str.cat([clus_high['Ad1_wt%']], sep=' ')
        clus_high_Ad2 = clus_high['Ad2'].str.cat([clus_high['Ad2_wt%']], sep=' ')
        clus_high_Ad3 = clus_high['Ad3'].str.cat([clus_high['Ad3_wt%']], sep=' ')
        clus_high_Ad4 = clus_high['Ad4'].str.cat([clus_high['Ad4_wt%']], sep=' ')
        clus_high_Ad5 = clus_high['Ad5'].str.cat([clus_high['Ad5_wt%']], sep=' ')
        if CalT_num == 0:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_PGM3, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_Ad5, clus_high_supp],  sep=', ')
        elif CalT_num == 1:
                clus_high_supp = clus_high.loc[:, 'Support'].astype(str)
                clus_high_supp = clus_high_supp.str.rstrip(',')
                clus_high_CalT = clus_high['Calc. temp. (℃)']
                clus_high['Top catal.'] = clus_high_PGM1.str.cat([clus_high_PGM2, clus_high_PGM3, clus_high_Ad1, clus_high_Ad2, clus_high_Ad3, clus_high_Ad4, clus_high_Ad5, clus_high_supp, clus_high_CalT],  sep=', ')
    clus_high[cand_wt_cols] = clus_high[cand_wt_cols].astype(float)
    return clus_high
