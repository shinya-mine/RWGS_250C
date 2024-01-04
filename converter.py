#!/usr/bin/env python
import modules.conditions as conditions
import modules.setting as setting
import numpy as np
import pandas as pd
import itertools
import scipy.sparse as sparse
from scipy.stats import norm
import collections
import sys
sys.dont_write_bytecode = True


# 多重リストを平滑化する関数。
def flatten(l):
    for el in l:
        if isinstance(el, collections.abc.Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el


# 任意の整数の約数一覧(昇順)を作成する。
def divisors_list(num):
    divisors = []
    for i in range(1, num+1):
        if num % i == 0:
            divisors.append(i)
    return divisors


# ある整数(total_num)を、ある大きさ(size)以上で割り切れる最小の整数を計算する。
def num_jobs(total_num, divisors, size=int(1E6)):
    catal_comb_size = []
    for i in range(len(divisors)):
        catal_comb_size.append(total_num/divisors[i])
        if catal_comb_size[i] > size:
            slice_size = int(catal_comb_size[i])
            num_slice = int(total_num / slice_size)
        if total_num < size:
            slice_size = total_num
            num_slice = 1
    return slice_size, num_slice


# 計算に用いる元素リストを作成する。
def elem_list(condition, converted, desc_cols):
    pgm_num, prm_num = condition['pgm_num'], condition['prm_num']
    add_num, supp_num = condition['add_num'], condition['supp_num']
    drop_pgms, drop_prms, drop_adds = condition['drop_pgms'], condition['drop_prms'], condition['drop_adds']
    essential_supp, drop_supp_list = condition['essential_supp'], condition['drop_supp_list']
    CalT_list, ReacT_list, WHSV_list = condition['CalT_list'], condition['ReacT_list'], condition['WHSV_list']
    H2_ppm_list, CO_ppm_list = condition['H2_ppm_list'], condition['CO_ppm_list']
    NH3_ppm_list, C3H6_ppm_list = condition['NH3_ppm_list'], condition['C3H6_ppm_list']
    CH4_ppm_list = condition['CH4_ppm_list']

    pgm_desc, prm_desc = converted['pgm_desc'], converted['prm_desc']
    add_desc, supp_desc = converted['add_desc'], converted['supp_desc']
    supp_name_dict = desc_cols['supp_name_dict']

    if pgm_num != 0:
        all_pgm_list = pgm_desc.index[1:]
        if len(drop_pgms) == 0:
            use_pgm_list = all_pgm_list
        else:
            use_pgm_list = list(set(all_pgm_list) - set(drop_pgms))
    else:
        all_pgm_list, use_pgm_list = [], []

    if prm_num != 0:
        all_prm_list = prm_desc.index[1:]
        if len(drop_prms) == 0:
            use_prm_list = all_prm_list
        else:
            use_prm_list = list(set(all_prm_list) - set(drop_prms))
    else:
        all_prm_list, use_prm_list = [], []

    if add_num != 0:
        all_add_list = add_desc.index[1:]
        if len(drop_adds) == 0:
            use_add_list = all_add_list
        else:
            use_add_list = list(set(all_add_list) - set(drop_adds))
    else:
        all_add_list, use_add_list = [], []

    if supp_num != 0:
        all_supp_list = supp_desc.index
        if len(drop_supp_list) == 0:
            use_supp_list = all_supp_list
        else:
            use_supp_list = list(set(all_supp_list) - set(drop_supp_list))
        if len(essential_supp) != 0:
            for i in range(len(essential_supp)):
                essential_supp[i] = supp_name_dict[essential_supp[i]]
            use_supp_list = essential_supp
    else:
        all_supp_list, use_supp_list = [], []

    elem_lists_dict = {
        'all_pgm': all_pgm_list, 'use_pgm': use_pgm_list,
        'all_prm': all_prm_list, 'use_prm': use_prm_list,
        'all_add': all_add_list, 'use_add': use_add_list,
        'all_supp': all_supp_list, 'use_supp': use_supp_list,
        'CalT': CalT_list, 'ReacT': ReacT_list, 'WHSV': WHSV_list,
        'H2_ppm': H2_ppm_list, 'CO_ppm': CO_ppm_list,
        'NH3_ppm': NH3_ppm_list, 'C3H6_ppm': C3H6_ppm_list,
        'CH4_ppm': CH4_ppm_list,
    }
    return elem_lists_dict


# 探索する触媒組成に応じたdictを作成(後の関数内で呼び出す)
def cols_order(condition):
    cols_order = {}
    cand_pgm_num, cand_prm_num = condition['cand_pgm_num'], condition['cand_prm_num']
    cand_add_num, cand_supp_num = condition['cand_add_num'], condition['cand_supp_num']
    CalT_num, ReacT_num, WHSV_num = condition['CalT_num'], condition['ReacT_num'], condition['WHSV_num']
    H2_ppm_num, CO_ppm_num = condition['H2_ppm_num'], condition['CO_ppm_num']
    NH3_ppm_num, C3H6_ppm_num = condition['NH3_ppm_num'], condition['C3H6_ppm_num']
    CH4_ppm_num = condition['CH4_ppm_num']

    if cand_pgm_num != 0 and cand_prm_num != 0 and cand_add_num != 0 \
        and cand_supp_num != 0 and CalT_num != 0 and ReacT_num != 0 and WHSV_num != 0:
        cols_order = {'pgm': 0, 'prm': 1, 'add': 2,
                    'supp': 3, 'CalT': 4, 'ReacT': 5, 'WHSV': 6}
    elif cand_pgm_num != 0 and cand_prm_num != 0 and cand_add_num != 0 \
        and cand_supp_num != 0 and CalT_num != 0 and ReacT_num != 0 and WHSV_num == 0:
        cols_order = {'pgm': 0, 'prm': 1, 'add': 2,
                    'supp': 3, 'CalT': 4, 'ReacT': 5}
    elif cand_pgm_num != 0 and cand_prm_num != 0 and cand_add_num != 0 \
        and cand_supp_num != 0 and CalT_num != 0 and ReacT_num == 0 and WHSV_num != 0:
        cols_order = {'pgm': 0, 'prm': 1, 'add': 2,
                    'supp': 3, 'CalT': 4, 'WHSV': 5}
    elif cand_pgm_num != 0 and cand_prm_num == 0 and cand_add_num != 0 \
        and cand_supp_num != 0 and CalT_num != 0 and ReacT_num == 0 and WHSV_num == 0:
        cols_order = {'pgm': 0, 'add': 1, 'supp': 2, 'CalT': 3}
    elif cand_pgm_num != 0 and cand_prm_num == 0 and cand_add_num != 0 \
        and cand_supp_num != 0 and CalT_num == 0 and ReacT_num == 0 and WHSV_num == 0:
        if H2_ppm_num != 0 or CO_ppm_num != 0 or NH3_ppm_num != 0 or C3H6_ppm_num != 0 or CH4_ppm_num != 0:
            cols_order = {'pgm':  0, 'add': 1, 'supp': 2,
                        'H2_ppm': 3, 'CO_ppm': 4, 'NH3_ppm': 5, 'C3H6_ppm': 6, 'CH4_ppm': 7}
        else:
            cols_order = {'pgm': 0, 'add': 1, 'supp': 2}
    elif cand_pgm_num == 0 and cand_prm_num == 0 and cand_add_num != 0 \
        and cand_supp_num != 0 and CalT_num != 0 and ReacT_num == 0 and WHSV_num == 0:
        cols_order = {'add': 0, 'supp': 1, 'CalT': 2}
    elif cand_pgm_num == 0 and cand_prm_num == 0 and cand_add_num != 0 \
        and cand_supp_num == 0 and CalT_num == 0 and ReacT_num == 0 and WHSV_num == 0:
        cols_order = {'add': 0}

    return cols_order


# 実験データを読み込む関数
def read_data(condition, desc_cols, ReacT_range=[], read_sheet_name='data_sheet_name', idx=None):
    date, Reaction = condition['date'], condition['Reaction']
    data_sheet_name, skip_rows = condition[read_sheet_name], condition['skip_rows']
    supp_num = condition['supp_num']
    supp_name_dict = desc_cols['supp_name_dict']
    data_file_name = f'data/{date}_{Reaction}_dataset.xlsx'

    data = pd.read_excel(data_file_name, sheet_name=data_sheet_name, skiprows=skip_rows)
    if Reaction == 'EtOH_CO2' or Reaction == 'EtOH_CO':
        if len(ReacT_range) == 2:
            data = data[data.loc[:, 'Reaction T'] >= ReacT_range[0]].reset_index()
            data = data[data.loc[:, 'Reaction T'] <= ReacT_range[1]].reset_index()
        elif len(ReacT_range) == 1:
            data = data[data.loc[:, 'Reaction T'] >= ReacT_range[0]].reset_index()

    data_cols = setting.data_columns(condition)
    elem_cols, wt_cols = data_cols['elem'], data_cols['wt']

    data.loc[:, elem_cols] = data.loc[:, elem_cols].fillna('H')
    data.loc[:, wt_cols] = data.loc[:, wt_cols].fillna(0)

    if supp_num != 0:
        for i in range(len(data)):
            data.loc[i, 'Support_name'] = supp_name_dict[data.loc[i, 'Support_name']]

    if Reaction == 'N2O-SCR':
        gas_cols = data_cols['gas_cols']
        data.loc[:, gas_cols] = data.loc[:, gas_cols].fillna(0)

    if idx is not None:
        #print(Reaction, 'Iteration =', idx)
        data_idx = data[data['Iteration'] <= idx]
        data_idx = data_idx.reset_index()
        data = data_idx
    return data


#　　Descriptorファイルを読み込む
def read_desc(elem, condition, desc_cols, select_desc=True, local=False, standard=False):
    Reaction = condition['Reaction']
    pgm_num, prm_num = condition['pgm_num'], condition['prm_num']
    add_num, supp_num = condition['add_num'], condition['supp_num']
    elem_comp_list = [pgm_num, prm_num, add_num, supp_num]
    desc_file_name = 'data/Descriptors.xlsx'

    if elem == 'pgm' or elem == 'prm' or elem == 'add' or elem == 'all':
        desc_sheet_name = 'Descriptors_elem'
        use_desc = desc_cols['basic_desc_columns']
        elem_desc_dict = desc_cols['elem_desc_dict']
        for i, j in enumerate(use_desc):
            if j in list(elem_desc_dict.keys()):
                use_desc[i] = elem_desc_dict[j]
        desc = pd.read_excel(desc_file_name, sheet_name=desc_sheet_name, index_col='Symbol')
        desc = desc.rename(columns=elem_desc_dict)
    elif elem == 'supp':
        desc_sheet_name = 'Descriptors_supp'
        use_desc = desc_cols['all_supp_desc']
        supp_desc_dict = desc_cols['supp_desc_dict']
        for i, j in enumerate(use_desc):
            if j in list(supp_desc_dict.keys()):
                use_desc[i] = supp_desc_dict[j]
        desc = pd.read_excel(desc_file_name, sheet_name=desc_sheet_name, index_col='Support_name_python')
        desc = desc.rename(columns=supp_desc_dict)

    all_elems = list(desc.index)
    noble_gas = desc_cols['noble_gas']
    drop_elems = desc_cols['drop_elems']

    if elem == 'pgm':
        use_elem = desc_cols['pgm_plus_ReAu']
        if Reaction == 'N2O':
            desc = desc.loc[use_elem].drop('Os', axis=0)
        else:
            desc = desc.loc[use_elem].drop(['Re', 'Os'], axis=0)
        desc = desc[use_desc].fillna(desc.mean())
        if len(elem_comp_list) - elem_comp_list.count(0) != 1:
            desc = desc.rename(columns=lambda s: s+' [PGMs]')
            aw = desc.loc[:, 'AW [PGMs]']
        else:
            aw = desc.loc[:, 'AW']

    elif elem == 'prm':
        use_elem = ['H'] + desc_cols['prm_elems']
        desc = desc.loc[use_elem]
        desc = desc[use_desc].fillna(desc.mean())
        if len(elem_comp_list) - elem_comp_list.count(0) != 1:
            desc = desc.rename(columns=lambda s: s+' [Promoters]')
            aw = desc.loc[:, 'AW [Promoters]']
        else:
            aw = desc.loc[:, 'AW']

    elif elem == 'add':
        use_elem = [i for i in all_elems if i not in noble_gas+drop_elems]
        desc = desc.loc[use_elem]
        desc = desc[use_desc].fillna(desc.mean())
        if len(elem_comp_list) - elem_comp_list.count(0) != 1:
            desc = desc.rename(columns=lambda s: s+' [Additives]')
            aw = desc.loc[:, 'AW [Additives]']
        else:
            aw = desc.loc[:, 'AW']

    elif elem == 'supp':
        use_elem = desc.index[desc[Reaction] == 1]
        desc = desc.loc[use_elem]
        #desc = desc.set_index('Support_name_python', inplace=True)
        desc = desc[use_desc].replace('-', 0)
        desc = desc[use_desc].fillna(desc.mean())
        if len(elem_comp_list) - elem_comp_list.count(0) != 1:
            desc = desc.rename(columns=lambda s: s+' [Supports]')
            aw = desc.loc[:, 'MW [Supports]']
        else:
            aw = desc.loc[:, 'MW']

    elif elem == 'all':
        desc = desc[use_desc].fillna(desc.mean())
        aw = desc.loc[:, 'AW']

    if local == False:
        selec_desc = conditions.desc_select(Reaction)
    # else:
        #selec_desc = desc_select(Reaction)

    if elem == 'pgm' or elem == 'prm' or elem == 'add':
        for i, j in enumerate(selec_desc['use_'+elem]):
            if j in list(elem_desc_dict.keys()):
                selec_desc['use_'+elem][i] = elem_desc_dict[j]
        for i, j in enumerate(selec_desc['drop_'+elem]):
            if j in list(elem_desc_dict.keys()):
                selec_desc['drop_'+elem][i] = elem_desc_dict[j]
    elif elem == 'supp':
        for i, j in enumerate(selec_desc['use_'+elem]):
            if j in list(supp_desc_dict.keys()):
                selec_desc['use_'+elem][i] = supp_desc_dict[j]
        for i, j in enumerate(selec_desc['drop_'+elem]):
            if j in list(supp_desc_dict.keys()):
                selec_desc['drop_'+elem][i] = supp_desc_dict[j]

    if select_desc == True:
        if elem == 'pgm':
            selec_desc['use_'+elem] = [s +' [PGMs]' for s in selec_desc['use_'+elem]]
        if elem == 'prm':
            selec_desc['use_'+elem] = [s +' [Promoters]' for s in selec_desc['use_'+elem]]
        if elem == 'add':
            selec_desc['use_'+elem] = [s +' [Additives]' for s in selec_desc['use_'+elem]]
        if elem == 'supp':
            selec_desc['use_'+elem] = [s +' [Supports]' for s in selec_desc['use_'+elem]]
        if elem != 'all':
            desc = desc.loc[:, selec_desc['use_'+elem]]
    else:
        if elem == 'pgm':
            selec_desc['drop_'+elem] = [s +' [PGMs]' for s in selec_desc['drop_'+elem]]
        if elem == 'prm':
            selec_desc['drop_'+elem] = [s +' [Promoters]' for s in selec_desc['drop_'+elem]]
        if elem == 'add':
            selec_desc['drop_'+elem] = [s +' [Additives]' for s in selec_desc['drop_'+elem]]
        if elem == 'supp':
            selec_desc['drop_'+elem] = [s +' [Supports]' for s in selec_desc['drop_'+elem]]
        if elem != 'all':
            desc = desc.drop(selec_desc['drop_'+elem], axis=1)
    desc = desc.dropna(how='any', axis=1)

    if standard == True:
        desc = (desc - desc.mean()) / desc.std(ddof=1)
    return desc, aw


# 触媒候補データ(cand_sum)を読み込む関数
def read_cand(condition):
    date, Reaction = condition['date'], condition['Reaction']
    ML_model, Search_method = condition['ML_model'], condition['Search_method']
    pgm_num, prm_num = condition['pgm_num'], condition['prm_num']
    add_num, supp_num = condition['add_num'], condition['supp_num']
    pgm_model, prm_model = condition['pgm_model'], condition['prm_model']
    add_model, supp_model = condition['add_model'], condition['supp_model']

    # results_csv_file_all
    if pgm_num != 0 and prm_num != 0 and add_num != 0 and supp_num != 0:
        cand_file_name = f'results/{date}_{Reaction}_cand_sum_{ML_model}_prop{pgm_model}{prm_model}{add_model}{supp_model}_{Search_method}.csv'
    elif pgm_num != 0 and prm_num == 0 and add_num != 0 and supp_num != 0:
        cand_file_name = f'results/{date}_{Reaction}_cand_sum_{ML_model}_prop{pgm_model}{add_model}{supp_model}_{Search_method}.csv'
    elif pgm_num == 0 and prm_num == 0 and add_num != 0 and supp_num != 0:
        cand_file_name = f'results/{date}_{Reaction}_cand_sum_{ML_model}_prop{add_model}{supp_model}_{Search_method}.csv'
    elif pgm_num == 0 and prm_num == 0 and add_num != 0 and supp_num == 0:
        cand_file_name = f'results/{date}_{Reaction}_cand_sum_{ML_model}_prop{add_model}_{Search_method}.csv'

    data_cols = setting.data_columns(condition)
    cand_elem, cand_wt = data_cols['cand_elem'], data_cols['cand_wt']

    cand = cand_separation(cand_file_name, condition)
    cand = cand_str(cand, condition, data_cols)
    cand.drop_duplicates(subset=['Top catal.'], inplace=True)

    cand.loc[:, cand_elem] = cand.loc[:, cand_elem].fillna('H')
    cand.loc[:, cand_wt] = cand.loc[:, cand_wt].fillna(0)
    cand = cand.reset_index()
    return cand


# Feat.を定義する関数
def make_feat(elem, elem_num, model, data, desc, aw):
    if elem == 'PGM' or elem == 'Promoter' or elem == 'Ad':
        if elem_num != 0:
            feat_comp = pd.DataFrame(
                index=np.arange(len(data)), columns=desc.index)
            for i in range(len(data)):
                for j in range(elem_num):
                    feat_comp.loc[i, data.loc[i, f'{elem}{j+1}']] = data.loc[i, f'{elem}{j+1}_wt%']
            feat_comp = feat_comp.fillna(0)
            if 'H' in desc.index:
                feat_comp = feat_comp.drop('H', axis=1)
                feat_comp_cols = list(desc.index[1:])
            else:
                feat_comp_cols = list(desc.index)

            feat_desc = sum([
                np.multiply(np.array(desc.loc[data.loc[:, f"{elem}{i+1}"]]),
                            np.array(data[f"{elem}{i+1}_wt%"]).reshape(-1, 1)) /
                np.array(aw.loc[data.loc[:, f"{elem}{i+1}"]]).reshape(-1, 1) for i in range(elem_num)
            ])
            feat_desc_cols = list(desc.columns)

            if model == 0:
                feat = feat_comp
                feat_cols = feat_comp_cols
            elif model == 1:
                feat = np.hstack((feat_comp, feat_desc))
                feat_cols = list(feat_comp_cols) + list(feat_desc_cols)
            elif model == 2:
                feat = feat_desc
                feat_cols = feat_desc_cols
        else:
            feat, feat_cols = None, None

    elif elem == 'Supp':
        if elem_num != 0:
            feat_comp = pd.DataFrame(index=np.arange(len(data)), columns=desc.index)
            for i in range(len(data)):
                feat_comp.loc[i, data.loc[i, 'Support_name']] = 1
            feat_comp = feat_comp.fillna(0)
            feat_comp_cols = list(desc.index)

            feat_desc = desc.loc[data.loc[:, 'Support_name']]
            feat_desc = feat_desc.replace('-', 0)
            feat_desc_cols = list(desc.columns)

            if model == 0:
                feat = feat_comp
                feat_cols = feat_comp_cols
            elif model == 1:
                feat = np.hstack((feat_comp, feat_desc))
                feat_cols = list(feat_comp_cols) + list(feat_desc_cols)
            elif model == 2:
                feat = feat_desc
                feat_cols = feat_desc_cols
        else:
            feat, feat_cols = None, None

    elif elem == 'Calc T' or elem == 'Reaction T' or elem == 'WHSV (mL g-1 h-1)' or \
            elem == 'H2_ppm' or elem == 'CO_ppm' or elem == 'NH3_ppm' or elem == 'C3H6_ppm' or elem == 'CH4_ppm':
        if elem_num != 0:
            feat = np.array(data.loc[:, elem]).reshape(-1, 1)
            feat_cols = [elem]
        else:
            feat, feat_cols = None, None
    return feat, feat_cols


# Targetを定義する関数
def make_target(Reaction, data, name, temp):
    if Reaction == 'rwgs_250' or Reaction == 'rwgs_250_1wt' or \
            Reaction == 'rwgs_300' or Reaction == 'rwgs_TL' or Reaction == 'CH3OH' or \
            Reaction == 'H2SCR' or Reaction == 'EtOH_CO2' or Reaction == 'EtOH_CO':
        target = data.loc[:, name]
        target = target.fillna(target.mean())
        target = target.astype(float)

    elif Reaction == 'N2O' or Reaction == 'N2O-SCR' or Reaction == 'CH4':
        target = data.loc[:, f'T{temp} (℃)']
        target = target.fillna(900)  # sufficiently high temperature
        target = target.astype(float)

    elif Reaction == 'NH3SCR':
        target = data.loc[:, f'T{temp}_N2 yield']
        target = target.replace('-', 400)
        target = target.astype(float)
    return target


# SWEDを定義する関数
def make_swed(elem, elem_num, data, desc, aw, data_cols):
    if elem_num != 0:
        elem_swed_names = []
        for i in range(elem_num):
            for j in desc.columns:
                elem_swed_names.append(f'{j} [{i+1}]')

        swed_elem = pd.DataFrame(columns=elem_swed_names)
        if elem != 'supp':
            for i in range(len(data)):
                swed_elem.loc[i] = np.array(
                    list(
                        flatten(
                            [(
                                np.array(desc.loc[data.loc[i, f"{elem}{j+1}"]]) *
                                np.array(data.loc[i, f"{elem}{j+1}_wt%"]) /
                                np.array(aw.loc[data.loc[i, f"{elem}{j+1}"]])
                            ) for j in range(elem_num)]
                        )))
        else:
            elem_wt_sum = data.loc[:, data_cols['wt']].sum(axis=1)
            supp_wt = 100 - elem_wt_sum
            for i in range(len(data)):
                swed_elem.loc[i] = np.array(desc.loc[data.loc[i, 'Support_name']]) * \
                    np.array(supp_wt[i]) / \
                    np.array(aw.loc[data.loc[i, 'Support_name']])
    else:
        swed_elem = None
    return swed_elem


# datasetと DescriptorのExcelファイルを読み込み、feat, target等必要な形状に加工する。
def data_convert(condition, data_cols, desc_cols, sheet_name='data_sheet_name',
                select_desc=True, use_models=[], local=False, idx=None):
    converted = {}
    data = read_data(condition, desc_cols, ReacT_range=[condition['ReacT_min'], condition['ReacT_max']],
                    read_sheet_name=sheet_name, idx=idx)

    pgm_desc, pgm_aw = read_desc('pgm', condition, desc_cols, select_desc, local, standard=condition['desc_std'])
    prm_desc, prm_aw = read_desc('prm', condition, desc_cols, select_desc, local, standard=condition['desc_std'])
    add_desc, add_aw = read_desc('add', condition, desc_cols, select_desc, local, standard=condition['desc_std'])
    supp_desc, supp_mw = read_desc('supp', condition, desc_cols, select_desc, local, standard=condition['desc_std'])

    if len(use_models) == 0:
        pgm_model, prm_model = condition['pgm_model'], condition['prm_model']
        add_model, supp_model = condition['add_model'], condition['supp_model']
    elif len(use_models) == 1:
        add_model = use_models[0]
        pgm_model, prm_model, supp_model = 2, 2, 2
    elif len(use_models) == 2:
        add_model, supp_model = use_models[0], use_models[1]
        pgm_model, prm_model = 2, 2
    elif len(use_models) == 3:
        pgm_model, add_model, supp_model = use_models[0], use_models[1], use_models[2]
        prm_model = 2
    elif len(use_models) == 4:
        pgm_model, prm_model = use_models[0], use_models[1]
        add_model, supp_model = use_models[0], use_models[1]

    feat_pgm, feat_pgm_cols = make_feat('PGM', condition['pgm_num'], pgm_model, data, pgm_desc, pgm_aw)
    feat_prm, feat_prm_cols = make_feat('Promoter', condition['prm_num'], prm_model, data, prm_desc, prm_aw)
    feat_add, feat_add_cols = make_feat('Ad', condition['add_num'], add_model, data, add_desc, add_aw)
    feat_supp, feat_supp_cols = make_feat('Supp', condition['supp_num'], supp_model, data, supp_desc, supp_mw)

    feat_CalT, feat_CalT_cols = make_feat('Calc T', condition['CalT_num'], 0, data, add_desc, add_aw)
    feat_ReacT, feat_ReacT_cols = make_feat('Reaction T', condition['ReacT_num'], 0, data, add_desc, add_aw)
    feat_WHSV, feat_WHSV_cols = make_feat('WHSV (mL g-1 h-1)', condition['WHSV_num'], 0, data, add_desc, add_aw)

    feat_H2, feat_H2_cols = make_feat('H2_ppm', condition['H2_ppm_num'], 0, data, add_desc, add_aw)
    feat_CO, feat_CO_cols = make_feat('CO_ppm', condition['CO_ppm_num'], 0, data, add_desc, add_aw)
    feat_NH3, feat_NH3_cols = make_feat('NH3_ppm', condition['NH3_ppm_num'], 0, data, add_desc, add_aw)
    feat_C3H6, feat_C3H6_cols = make_feat('C3H6_ppm', condition['C3H6_ppm_num'], 0, data, add_desc, add_aw)
    feat_CH4, feat_CH4_cols = make_feat('CH4_ppm', condition['CH4_ppm_num'], 0, data, add_desc, add_aw)

    feat_list = [
        feat_pgm, feat_prm, feat_add,
        feat_supp, feat_CalT, feat_ReacT, feat_WHSV,
        feat_H2, feat_CO, feat_NH3, feat_C3H6, feat_CH4
    ]
    feat_cols_list = [
        feat_pgm_cols, feat_prm_cols, feat_add_cols, feat_supp_cols,
        feat_CalT_cols, feat_ReacT_cols, feat_WHSV_cols,
        feat_H2_cols, feat_CO_cols, feat_NH3_cols, feat_C3H6_cols, feat_CH4_cols
    ]
    feat_list = [i for i in feat_list if i is not None]
    feat = np.hstack(feat_list)
    feat_cols_list = [i for i in feat_cols_list if i is not None]
    feat_cols = list(itertools.chain.from_iterable(feat_cols_list))
    #feat_cols = [item for row in feat_cols_list for item in row]
    feat = pd.DataFrame(feat, columns=feat_cols)

    swed_pgm = make_swed('PGM', condition['pgm_num'], data, pgm_desc, pgm_aw, data_cols)
    swed_prm = make_swed('Promoter', condition['prm_num'], data, prm_desc, prm_aw, data_cols)
    swed_add = make_swed('Ad', condition['add_num'], data, add_desc, add_aw, data_cols)
    swed_supp = make_swed('supp', condition['supp_num'], data, supp_desc, supp_mw, data_cols)

    target = make_target(condition['Reaction'], data, condition['target_name'], condition['target_temp'])

    converted = {
        'data': data,
        'pgm_desc': pgm_desc, 'prm_desc': prm_desc, 'add_desc': add_desc, 'supp_desc': supp_desc,
        'pgm_aw': pgm_aw, 'prm_aw': prm_aw, 'add_aw': add_aw, 'supp_mw': supp_mw,
        'swed_pgm': swed_pgm, 'swed_prm': swed_prm, 'swed_add': swed_add, 'swed_supp': swed_supp,
        'feat': feat, 'target': target
    }
    return converted


# datasetと DescriptorのExcelファイルを読み込み、feat, target等必要な形状に加工する。
def cand_convert(condition, converted, data_cols, desc_cols, select_desc=True, use_models=[], local=False):
    cand_dict = {}
    
    data = converted['data']
    data = data_str(data, condition, data_cols)

    cand = read_cand(condition)
    dupl_cat = list(set(data['Catalyst']) & set(cand['Top catal.']))
    if len(dupl_cat) != 0:
        for i in range(len(dupl_cat)):
            cand.drop(list(cand[cand['Top catal.'] == dupl_cat[i]].index), inplace=True)
    cand.reset_index(drop=True, inplace=True)

    for i in list(cand.columns[:-1]): # drop 'Top catal.'
        if type(cand.loc[0, i]) == str and ',' in cand.loc[0, i]:
            cand[i] = cand[i].str.strip(',')
            if i != 'Support_name':
                cand[i] = cand[i].astype(float)

    pgm_desc, pgm_aw = read_desc('pgm', condition, desc_cols, select_desc, local, standard=condition['desc_std'])
    prm_desc, prm_aw = read_desc('prm', condition, desc_cols, select_desc, local, standard=condition['desc_std'])
    add_desc, add_aw = read_desc('add', condition, desc_cols, select_desc, local, standard=condition['desc_std'])
    supp_desc, supp_mw = read_desc('supp', condition, desc_cols, select_desc, local, standard=condition['desc_std'])

    if len(use_models) == 0:
        pgm_model, prm_model = condition['pgm_model'], condition['prm_model']
        add_model, supp_model = condition['add_model'], condition['supp_model']
    elif len(use_models) == 1:
        add_model = use_models[0]
        pgm_model, prm_model, supp_model = 2, 2, 2
    elif len(use_models) == 2:
        add_model, supp_model = use_models[0], use_models[1]
        pgm_model, prm_model = 2, 2
    elif len(use_models) == 3:
        pgm_model, add_model, supp_model = use_models[0], use_models[1], use_models[2]
        prm_model = 2
    elif len(use_models) == 4:
        pgm_model, prm_model = use_models[0], use_models[1]
        add_model, supp_model = use_models[0], use_models[1]

    feat_pgm, feat_pgm_cols = make_feat('PGM', condition['cand_pgm_num'], pgm_model, cand, pgm_desc, pgm_aw)
    feat_prm, feat_prm_cols = make_feat('Promoter', condition['cand_prm_num'], prm_model, cand, prm_desc, prm_aw)
    feat_add, feat_add_cols = make_feat('Ad', condition['cand_add_num'], add_model, cand, add_desc, add_aw)
    feat_supp, feat_supp_cols = make_feat('Supp', condition['cand_supp_num'], supp_model, cand, supp_desc, supp_mw)

    feat_CalT, feat_CalT_cols = make_feat('Calc T', condition['CalT_num'], 0, cand, add_desc, add_aw)
    feat_ReacT, feat_ReacT_cols = make_feat('Reaction T', condition['ReacT_num'], 0, cand, add_desc, add_aw)
    feat_WHSV, feat_WHSV_cols = make_feat('WHSV (mL g-1 h-1)', condition['WHSV_num'], 0, cand, add_desc, add_aw)

    feat_H2, feat_H2_cols = make_feat('H2_ppm', condition['H2_ppm_num'], 0, cand, add_desc, add_aw)
    feat_CO, feat_CO_cols = make_feat('CO_ppm', condition['CO_ppm_num'], 0, cand, add_desc, add_aw)
    feat_NH3, feat_NH3_cols = make_feat('NH3_ppm', condition['NH3_ppm_num'], 0, cand, add_desc, add_aw)
    feat_C3H6, feat_C3H6_cols = make_feat('C3H6_ppm', condition['C3H6_ppm_num'], 0, cand, add_desc, add_aw)
    feat_CH4, feat_CH4_cols = make_feat('CH4_ppm', condition['CH4_ppm_num'], 0, cand, add_desc, add_aw)

    feat_list = [
        feat_pgm, feat_prm, feat_add,
        feat_supp, feat_CalT, feat_ReacT, feat_WHSV,
        feat_H2, feat_CO, feat_NH3, feat_C3H6, feat_CH4
    ]
    feat_cols_list = [
        feat_pgm_cols, feat_prm_cols, feat_add_cols, feat_supp_cols,
        feat_CalT_cols, feat_ReacT_cols, feat_WHSV_cols,
        feat_H2_cols, feat_CO_cols, feat_NH3_cols, feat_C3H6_cols, feat_CH4_cols
    ]
    feat_list = [i for i in feat_list if i is not None]
    feat = np.hstack(feat_list)
    feat_cols_list = [i for i in feat_cols_list if i is not None]
    feat_cols = list(itertools.chain.from_iterable(feat_cols_list))
    #feat_cols = [item for row in feat_cols_list for item in row]
    feat = pd.DataFrame(feat, columns=feat_cols)
    target = cand.loc[:, 'ei']

    cand_dict = {
        'cand': cand,
        'pgm_desc': pgm_desc, 'prm_desc': prm_desc, 'add_desc': add_desc, 'supp_desc': supp_desc,
        'pgm_aw': pgm_aw, 'prm_aw': prm_aw, 'add_aw': add_aw, 'supp_mw': supp_mw,
        'feat': feat, 'target': target
        }
    return cand_dict


# 各元素成分ごとに、wtの組み合わせを作る関数
def elem_wt_comb(elem_wt, ess_elem_wt, cand_elem_num, fix_elem_num, max_elem_wt):
    if cand_elem_num != 0:
        elem_wt_comb_full = itertools.product(elem_wt, repeat=cand_elem_num-fix_elem_num)
        if fix_elem_num != 0:
            elem_wt_comb_fix = itertools.product(ess_elem_wt, repeat=fix_elem_num)
            elem_wt_comb = itertools.product(elem_wt_comb_fix, elem_wt_comb_full)
        else:
            elem_wt_comb = itertools.product(elem_wt, repeat=cand_elem_num)

        elem_wt_comb = np.array(list(flatten(elem_wt_comb))).reshape(-1, cand_elem_num)
        if max_elem_wt is not None:
            elem_wt_comb = elem_wt_comb[elem_wt_comb.sum(axis=1) < max_elem_wt]
        elem_wt_comb = iter(elem_wt_comb)
    else:
        elem_wt_comb = None
    return elem_wt_comb


# 上で定義したelem_wt_comb関数を成分の数だけ結合し、wtの全ての組み合わせを作成する。
def all_wt_comb(condition):
    pgm_wt, prm_wt, add_wt = condition['pgm_wt'], condition['prm_wt'], condition['add_wt']
    ess_pgm_wt, ess_prm_wt, ess_add_wt = condition['essential_pgm_wt'], condition['essential_prm_wt'], condition['essential_add_wt']
    cand_pgm_num, cand_prm_num = condition['cand_pgm_num'], condition['cand_prm_num'],
    cand_add_num, cand_supp_num = condition['cand_add_num'], condition['cand_supp_num']
    fix_pgm_num, fix_prm_num, fix_add_num = condition['fix_pgm_num'], condition['fix_prm_num'], condition['fix_add_num']
    max_pgm_wt, max_prm_wt, max_add_wt = condition['max_pgm_wt'], condition['max_prm_wt'], condition['max_add_wt']

    elem_wt_pgm = elem_wt_comb(pgm_wt, ess_pgm_wt, cand_pgm_num, fix_pgm_num, max_pgm_wt)
    elem_wt_prm = elem_wt_comb(prm_wt, ess_prm_wt, cand_prm_num, fix_prm_num, max_prm_wt)
    elem_wt_add = elem_wt_comb(add_wt, ess_add_wt, cand_add_num, fix_add_num, max_add_wt)

    _lis_wt = [elem_wt_pgm, elem_wt_prm, elem_wt_add]
    comb_wt = [z for z in _lis_wt if z]
    all_comb_wt = itertools.product(*comb_wt)
    contens = len(_lis_wt) - _lis_wt.count(None)
    if contens == 1:
        elem_wt = [x1[0] for x1 in all_comb_wt]
    elif contens == 2:
        elem_wt = [(*x1, *x2) for (x1, x2) in all_comb_wt]
    elif contens == 3:
        elem_wt = [(*x1, *x2, *x3) for (x1, x2, x3) in all_comb_wt]
    if cand_supp_num != 0:
        elem_wt = [elem_wt[i] + tuple([100-sum(elem_wt[i])]) for i in range(len(elem_wt))]
    num_elem_wt = len(elem_wt)
    return elem_wt, num_elem_wt


# 必須元素が登場するリストとその他の元素(未登場の必須元素も含まれ得る)が登場するリストをそれぞれ作成。
def gen_fix_and_res(use_elem_list, essential_elems, cand_elem_num, fix_elem_num, duplicate=True):
    # duplicate: True: 重複を含む, False: 重複を含まない
    residue_elems = list(set(use_elem_list) - set(essential_elems))
    elem_comb_fix = list(itertools.combinations(essential_elems, fix_elem_num))
    if duplicate == True:
        elem_comb_res = list(itertools.combinations(
            use_elem_list, cand_elem_num-fix_elem_num))
    else:
        elem_comb_res = list(itertools.combinations(
            residue_elems, cand_elem_num-fix_elem_num))
    return elem_comb_fix, elem_comb_res


# 上のgen_fix_and_res関数で作成した必須元素リストと、その他元素リストを結合する。
def fix_and_res_comb(fix_list, res_list):
    _lis = [fix_list, res_list]
    comb = [z for z in _lis if z is not None]
    comb_ = itertools.product(*comb)
    comb_list = [(*x1, *x2) for (x1, x2) in comb_]
    return comb_list


# 上のfix_and_res_comb関数で作成したリストに含まれる重複(e.g. ['Pt', 'Pt'])成分を落とす。
def drop_duplicate(elem_comb, fix_elem_num):
    elem_comb_dup = []
    for i in range(len(elem_comb)):
        if len(list(set(elem_comb[i][:fix_elem_num]) & set(elem_comb[i][fix_elem_num:]))) == 0:
            elem_comb_dup.append(elem_comb[i])
    return elem_comb_dup


# drop_duplicate関数で作成した重複を含まないリストを要素ごとに分解する。
def make_elem_comp(elem_comb, cand_elem_num):
    elem_comp_dict = {}
    for i in range(cand_elem_num):
        elem_comp_dict[i] = [x[i] for x in elem_comb]
    return elem_comp_dict


# 触媒組成 (+ 実験条件) の全組み合わせを計算する。
def catal_comb(condition, elem_lists):
    cand_pgm_num, cand_prm_num = condition['cand_pgm_num'], condition['cand_prm_num']
    cand_add_num, cand_supp_num = condition['cand_add_num'], condition['cand_supp_num']
    CalT_num, ReacT_num, WHSV_num = condition['CalT_num'], condition['ReacT_num'], condition['WHSV_num']
    H2_ppm_num, CO_ppm_num = condition['H2_ppm_num'], condition['CO_ppm_num']
    NH3_ppm_num, C3H6_ppm_num = condition['NH3_ppm_num'], condition['C3H6_ppm_num']
    CH4_ppm_num = condition['CH4_ppm_num']
    fix_pgm_num, fix_prm_num, fix_add_num = condition['fix_pgm_num'], condition['fix_prm_num'], condition['fix_add_num']
    essential_pgms, essential_prms, essential_adds = condition['essential_pgms'], condition['essential_prms'], condition['essential_adds']
    use_pgm_list, use_prm_list = elem_lists['use_pgm'], elem_lists['use_prm']
    use_add_list, use_supp_list = elem_lists['use_add'], elem_lists['use_supp']
    CalT_list, ReacT_list, WHSV_list = elem_lists['CalT'], elem_lists['ReacT'], elem_lists['WHSV']
    H2_ppm_list, CO_ppm_list = elem_lists['H2_ppm'], elem_lists['CO_ppm']
    NH3_ppm_list, C3H6_ppm_list = elem_lists['NH3_ppm'], elem_lists['C3H6_ppm']
    CH4_ppm_list = elem_lists['CH4_ppm']

    if cand_pgm_num != 0:
        if fix_pgm_num != 0:
            elem_comb_pgm_fix, elem_comb_pgm_res = \
                gen_fix_and_res(use_pgm_list, essential_pgms, cand_pgm_num, fix_pgm_num, duplicate=True)
            elem_comb_pgm = fix_and_res_comb(elem_comb_pgm_fix, elem_comb_pgm_res)
            del elem_comb_pgm_fix
            del elem_comb_pgm_res
            elem_comb_pgm = drop_duplicate(elem_comb_pgm, fix_pgm_num)
        else:
            elem_comb_pgm = itertools.combinations(use_pgm_list, cand_pgm_num)
            elem_comb_pgm = list(elem_comb_pgm)
        comb_pgm_num = len(elem_comb_pgm)
    else:
        elem_comb_pgm = None
        comb_pgm_num = 1

    if cand_prm_num != 0:
        if fix_prm_num != 0:
            elem_comb_prm_fix, elem_comb_prm_res = \
                gen_fix_and_res(use_prm_list, essential_prms, cand_prm_num, fix_prm_num, duplicate=True)
            elem_comb_prm = fix_and_res_comb(elem_comb_prm_fix, elem_comb_prm_res)
            del elem_comb_prm_fix
            del elem_comb_prm_res
            elem_comb_prm = drop_duplicate(elem_comb_prm, fix_prm_num)
        else:
            elem_comb_prm = itertools.combinations(use_prm_list, cand_prm_num)
            elem_comb_prm = list(elem_comb_prm)
        comb_prm_num = len(elem_comb_prm)
    else:
        elem_comb_prm = None
        comb_prm_num = 1

    if cand_add_num != 0:
        if fix_add_num != 0:
            elem_comb_add_fix, elem_comb_add_res \
                = gen_fix_and_res(use_add_list, essential_adds, cand_add_num, fix_add_num, duplicate=True)
            elem_comb_add = fix_and_res_comb(elem_comb_add_fix, elem_comb_add_res)
            del elem_comb_add_fix
            del elem_comb_add_res
            elem_comb_add = drop_duplicate(elem_comb_add, fix_add_num)
        else:
            elem_comb_add = itertools.combinations(use_add_list, cand_add_num)
            elem_comb_add = list(elem_comb_add)
        comb_add_num = len(elem_comb_add)
    else:
        elem_comb_add = None
        comb_add_num = 1

    if cand_supp_num != 0:
        elem_comb_supp = use_supp_list
        comb_supp_num = len(elem_comb_supp)
    else:
        elem_comb_supp = None
        comb_supp_num = 1

    if CalT_num != 0:
        elem_comb_CalT = CalT_list
        comb_CalT_num = len(elem_comb_CalT)
    else:
        elem_comb_CalT = None
        comb_CalT_num = 1

    if ReacT_num != 0:
        elem_comb_ReacT = ReacT_list
        comb_ReacT_num = len(elem_comb_ReacT)
    else:
        elem_comb_ReacT = None
        comb_ReacT_num = 1

    if WHSV_num != 0:
        elem_comb_WHSV = WHSV_list
        comb_WHSV_num = len(elem_comb_WHSV)
    else:
        elem_comb_WHSV = None
        comb_WHSV_num = 1

    if H2_ppm_num != 0:
        elem_comb_H2 = H2_ppm_list
        comb_H2_num = len(elem_comb_H2)
    else:
        elem_comb_H2 = None
        comb_H2_num = 1

    if CO_ppm_num != 0:
        elem_comb_CO = CO_ppm_list
        comb_CO_num = len(elem_comb_CO)
    else:
        elem_comb_CO = None
        comb_CO_num = 1

    if NH3_ppm_num != 0:
        elem_comb_NH3 = NH3_ppm_list
        comb_NH3_num = len(elem_comb_NH3)
    else:
        elem_comb_NH3 = None
        comb_NH3_num = 1

    if C3H6_ppm_num != 0:
        elem_comb_C3H6 = C3H6_ppm_list
        comb_C3H6_num = len(elem_comb_C3H6)
    else:
        elem_comb_C3H6 = None
        comb_C3H6_num = 1

    if CH4_ppm_num != 0:
        elem_comb_CH4 = CH4_ppm_list
        comb_CH4_num = len(elem_comb_CH4)
    else:
        elem_comb_CH4 = None
        comb_CH4_num = 1

    total_comb_num = comb_pgm_num * comb_prm_num * comb_add_num \
        * comb_supp_num * comb_CalT_num * comb_ReacT_num * comb_WHSV_num \
        * comb_H2_num * comb_CO_num * comb_NH3_num * comb_C3H6_num * comb_CH4_num

    _lis = [
        elem_comb_pgm, elem_comb_prm, elem_comb_add,
        elem_comb_supp, elem_comb_CalT, elem_comb_ReacT, elem_comb_WHSV,
        elem_comb_H2, elem_comb_CO, elem_comb_NH3, elem_comb_C3H6, elem_comb_CH4
    ]
    comb = [z for z in _lis if z is not None]
    all_comb = itertools.product(*comb)
    return all_comb, total_comb_num


# 触媒組成を one-hot表現に変換する関数
def one_hot(read_col_name, elem_list, cand_elem_num, cand_cols_order, condition, elem_lists, slice_size, start, end):
    elem_comp_dict = {}
    elem_cols = {}
    ones = np.ones(slice_size)
    rows = np.arange(0, slice_size, 1)
    int_code = {x: j for (j, x) in enumerate(elem_list)}
    read_col = cand_cols_order[read_col_name]

    if cand_elem_num == 1:
        all_comb, total_comb_num = catal_comb(condition, elem_lists)
        all_comb_slice = itertools.islice(all_comb, start, end)
        elem_comp_dict = [u[read_col] for u in all_comb_slice]
        elem_cols = np.asarray([int_code[v] for v in elem_comp_dict])
        elem_oh = [sparse.csr_matrix(
            (ones, (rows, elem_cols)), shape=(slice_size, len(elem_list)))]
    else:
        for n in range(cand_elem_num):
            all_comb, total_comb_num = catal_comb(condition, elem_lists)
            all_comb_slice = itertools.islice(all_comb, start, end)
            elem_comp_dict[n] = [u[read_col][n] for u in all_comb_slice]
            elem_cols[n] = np.asarray([int_code[v] for v in elem_comp_dict[n]])
        elem_oh = [sparse.csr_matrix((ones, (rows, elem_cols[n])), shape=(
            slice_size, len(elem_list))) for n in range(cand_elem_num)]
    return elem_oh


# 実験条件を入力表現に変換する関数
def exp_rep(read_col_name, cand_cols_order, condition, elem_lists, start, end):
    read_col = cand_cols_order[read_col_name]
    all_comb, total_comb_num = catal_comb(condition, elem_lists)
    all_comb_slice = itertools.islice(all_comb, start, end)
    exp_dict = [u[read_col] for u in all_comb_slice]
    exp_col = np.array(exp_dict).reshape(-1, 1)
    exp_rep = sparse.csr_matrix(exp_col)
    return exp_rep


# 元素記述子を入力表現に変換する関数
def elem_rep(read_col_name, desc, aw, cand_elem_num, cand_cols_order, condition, elem_lists, start, end):
    elem_dict = {}
    read_col = cand_cols_order[read_col_name]
    all_comb, total_comb_num = catal_comb(condition, elem_lists)
    all_comb_slice = itertools.islice(all_comb, start, end)
    elem_comp_dict = [u[read_col] for u in all_comb_slice]

    if cand_elem_num == 1:
        elem_rep = [np.array(desc.loc[elem_comp_dict]) /
                    np.array(aw.loc[elem_comp_dict]).reshape(-1, 1)]
        elem_rep = [sparse.csr_matrix(m) for m in elem_rep]

    else:
        elem_dict = [np.asarray([i[n] for i in elem_comp_dict]) for n in range(cand_elem_num)]
        elem_rep = [np.array(desc.loc[elem_dict[i]]) /
                    np.array(aw.loc[elem_dict[i]]).reshape(-1, 1) for i in range(cand_elem_num)]
        elem_rep = [sparse.csr_matrix(m) for m in elem_rep]
    return elem_rep


# 上で定義した3種の関数(one_hot, exp_rep, elem_rep)を使って、EI値計算のための最終的な入力表現を作成する関数
def input_comps(condition, converted, elem_lists, cand_cols_order, slice_size, start, end):
    pgm_model, prm_model = condition['pgm_model'], condition['prm_model']
    add_model, supp_model = condition['add_model'], condition['supp_model']
    cand_pgm_num, cand_prm_num = condition['cand_pgm_num'], condition['cand_prm_num']
    cand_add_num, cand_supp_num = condition['cand_add_num'], condition['supp_num']
    CalT_num, ReacT_num, WHSV_num = condition['CalT_num'], condition['ReacT_num'], condition['WHSV_num']
    H2_ppm_num, CO_ppm_num = condition['H2_ppm_num'], condition['CO_ppm_num'],
    NH3_ppm_num, C3H6_ppm_num = condition['NH3_ppm_num'], condition['C3H6_ppm_num']
    CH4_ppm_num = condition['CH4_ppm_num']

    pgm_desc, prm_desc = converted['pgm_desc'], converted['prm_desc'],
    add_desc, supp_desc = converted['add_desc'], converted['supp_desc']
    pgm_aw, prm_aw = converted['pgm_aw'], converted['prm_aw']
    add_aw, supp_mw = converted['add_aw'], converted['supp_mw']
    all_pgm_list, all_prm_list = elem_lists['all_pgm'], elem_lists['all_prm'],
    all_add_list, all_supp_list = elem_lists['all_add'], elem_lists['all_supp']

    all_comb, total_comb_num = catal_comb(condition, elem_lists)
    comb_slice = itertools.islice(all_comb, start, end)
    comb_slice_list = list(comb_slice)

    if pgm_model == 0 and cand_pgm_num != 0:
        pgm_rep = None
        pgm_oh = one_hot('pgm', all_pgm_list, cand_pgm_num, cand_cols_order,
                        condition, elem_lists, slice_size, start, end)
    elif pgm_model == 1 and cand_pgm_num != 0:
        pgm_rep = elem_rep('pgm', pgm_desc, pgm_aw, cand_pgm_num,
                        cand_cols_order, condition, elem_lists, start, end)
        pgm_oh = one_hot('pgm', all_pgm_list, cand_pgm_num, cand_cols_order,
                        condition, elem_lists, slice_size, start, end)
    elif pgm_model == 2 and cand_pgm_num != 0:
        pgm_rep = elem_rep('pgm', pgm_desc, pgm_aw, cand_pgm_num,
                        cand_cols_order, condition, elem_lists, start, end)
        pgm_oh = None
    else:
        pgm_rep, pgm_oh = None, None

    if prm_model == 0 and cand_prm_num != 0:
        prm_rep = None
        prm_oh = one_hot('prm', all_prm_list, cand_prm_num, cand_cols_order,
                        condition, elem_lists, slice_size, start, end)
    elif prm_model == 1 and cand_prm_num != 0:
        prm_rep = elem_rep('prm', prm_desc, prm_aw, cand_prm_num,
                        cand_cols_order, condition, elem_lists, start, end)
        prm_oh = one_hot('prm', all_prm_list, cand_prm_num, cand_cols_order,
                        condition, elem_lists, slice_size, start, end)
    elif prm_model == 2 and cand_prm_num != 0:
        prm_rep = elem_rep('prm', prm_desc, prm_aw, cand_prm_num,
                        cand_cols_order, condition, elem_lists, start, end)
        prm_oh = None
    else:
        prm_rep, prm_oh = None, None

    if add_model == 0 and cand_add_num != 0:
        add_rep = None
        add_oh = one_hot('add', all_add_list, cand_add_num, cand_cols_order,
                        condition, elem_lists, slice_size, start, end)
    elif add_model == 1 and cand_add_num != 0:
        add_rep = elem_rep('add', add_desc, add_aw, cand_add_num,
                        cand_cols_order, condition, elem_lists, start, end)
        add_oh = one_hot('add', all_add_list, cand_add_num, cand_cols_order,
                        condition, elem_lists, slice_size, start, end)
    elif add_model == 2 and cand_add_num != 0:
        add_rep = elem_rep('add', add_desc, add_aw, cand_add_num,
                        cand_cols_order, condition, elem_lists, start, end)
        add_oh = None
    else:
        add_rep, add_oh = None, None

    if supp_model == 0 and cand_supp_num != 0:
        supp_rep = None
        supp_oh = one_hot('supp', all_supp_list, cand_supp_num,
                        cand_cols_order, condition, elem_lists, slice_size, start, end)
    elif supp_model == 1 and cand_supp_num != 0:
        supp_rep = elem_rep('supp', supp_desc, supp_mw, cand_supp_num,
                        cand_cols_order, condition, elem_lists, start, end)
        supp_oh = one_hot('supp', all_supp_list, cand_supp_num,
                        cand_cols_order, condition, elem_lists, slice_size, start, end)
    elif supp_model == 2 and cand_supp_num != 0:
        supp_rep = elem_rep('supp', supp_desc, supp_mw, cand_supp_num,
                        cand_cols_order, condition, elem_lists, start, end)
        supp_oh = None
    else:
        supp_rep, supp_oh = None, None

    if CalT_num != 0:
        CalT_rep = exp_rep('CalT', cand_cols_order,
                        condition, elem_lists, start, end)
    else:
        CalT_rep = None

    if ReacT_num != 0:
        ReacT_rep = exp_rep('ReacT', cand_cols_order,
                            condition, elem_lists, start, end)
    else:
        ReacT_rep = None

    if WHSV_num != 0:
        WHSV_rep = exp_rep('WHSV', cand_cols_order,
                        condition, elem_lists, start, end)
    else:
        WHSV_rep = None

    if H2_ppm_num != 0:
        H2_ppm_rep = exp_rep('H2_ppm', cand_cols_order,
                        condition, elem_lists, start, end)
    else:
        H2_ppm_rep = None

    if CO_ppm_num != 0:
        CO_ppm_rep = exp_rep('CO_ppm', cand_cols_order,
                        condition, elem_lists, start, end)
    else:
        CO_ppm_rep = None

    if NH3_ppm_num != 0:
        NH3_ppm_rep = exp_rep('NH3_ppm', cand_cols_order,
                        condition, elem_lists, start, end)
    else:
        NH3_ppm_rep = None

    if C3H6_ppm_num != 0:
        C3H6_ppm_rep = exp_rep('C3H6_ppm', cand_cols_order,
                        condition, elem_lists, start, end)
    else:
        C3H6_ppm_rep = None

    if CH4_ppm_num != 0:
        CH4_ppm_rep = exp_rep('CH4_ppm', cand_cols_order,
                        condition, elem_lists, start, end)
    else:
        CH4_ppm_rep = None

    return comb_slice_list, pgm_rep, pgm_oh, prm_rep, prm_oh, add_rep, add_oh, supp_rep, supp_oh, \
        CalT_rep, ReacT_rep, WHSV_rep, H2_ppm_rep, CO_ppm_rep, NH3_ppm_rep, C3H6_ppm_rep, CH4_ppm_rep


# 標準偏差を計算する関数　(鈴木氏作成, 峯改良)
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


# EI値を計算する関数　(鈴木氏作成, 峯改良)
def EI(Reaction, mu, sigma, cur_min, cur_max):
    if Reaction == 'N2O' or Reaction == 'N2O-SCR' or Reaction == 'NH3SCR' or Reaction == 'CH4':
        Z = (cur_min - mu) / sigma
        ei = (cur_min - mu) * norm.cdf(Z) + sigma*norm.pdf(Z)
    else:
        Z = (mu - cur_max) / sigma
        ei = (mu - cur_max) * norm.cdf(Z) + sigma*norm.pdf(Z)
    return ei


# EI値を計算した結果、上位に残った触媒組成+実験条件のリストの形を整形する関数(encodeしbarhグラフをプロットするために使用)
def make_cand_comp(cand_elem, condition):
    cand_pgm_num, cand_prm_num = condition['cand_pgm_num'], condition['cand_prm_num']
    cand_add_num, cand_supp_num = condition['cand_add_num'], condition['cand_supp_num']
    CalT_num, ReacT_num, WHSV_num = condition['CalT_num'], condition['ReacT_num'], condition['WHSV_num']
    H2_ppm_num, CO_ppm_num = condition['H2_ppm_num'], condition['CO_ppm_num']
    NH3_ppm_num, C3H6_ppm_num = condition['NH3_ppm_num'], condition['C3H6_ppm_num']
    CH4_ppm_num = condition['CH4_ppm_num']

    total_elem_num = cand_pgm_num + cand_prm_num + cand_add_num + cand_supp_num
    total_elem_exp_num = total_elem_num + CalT_num + ReacT_num + WHSV_num
    cand_cols_order = cols_order(condition)
    elem_comp_all = {}
    elem_comp_dict = {}

    elem_comp_pgm = {}
    if cand_pgm_num != 0:
        for i in range(0, cand_pgm_num):
            elem_comp_dict[i] = [x[cand_cols_order['pgm']][i] for x in cand_elem]
            elem_comp_pgm[i] = np.array([x for x in elem_comp_dict[i]]).reshape(-1, 1)
        elem_comp_all['pgm'] = np.hstack([elem_comp_pgm[i] for i in range(cand_pgm_num)]).reshape(-1, cand_pgm_num)

    elem_comp_prm = {}
    if cand_prm_num != 0:
        for i in range(0, cand_prm_num):
            elem_comp_dict[i+cand_pgm_num] = [x[cand_cols_order['prm']][i] for x in cand_elem]
            elem_comp_prm[i] = np.array([x for x in elem_comp_dict[i+cand_pgm_num]]).reshape(-1, 1)
        elem_comp_all['prm'] = np.hstack([elem_comp_prm[i] for i in range(cand_prm_num)]).reshape(-1, cand_prm_num)

    elem_comp_add = {}
    if cand_add_num != 0:
        for i in range(0, cand_add_num):
            elem_comp_dict[i+cand_pgm_num+cand_prm_num] = [x[cand_cols_order['add']][i] for x in cand_elem]
            elem_comp_add[i] = np.array([x for x in elem_comp_dict[i+cand_pgm_num+cand_prm_num]]).reshape(-1, 1)
        elem_comp_all['add'] = np.hstack([elem_comp_add[i] for i in range(cand_add_num)]).reshape(-1, cand_add_num)

    elem_comp_supp = {}
    if cand_supp_num != 0:
        for i in range(0, cand_supp_num):
            elem_comp_dict[i+cand_pgm_num+cand_prm_num+cand_add_num] = [x[cand_cols_order['supp']] for x in cand_elem]
            elem_comp_supp[i] = np.array([x for x in elem_comp_dict[i+cand_pgm_num+cand_prm_num+cand_add_num]]).reshape(-1, 1)
        elem_comp_all['supp'] = np.hstack([elem_comp_supp[i] for i in range(cand_supp_num)]).reshape(-1, cand_supp_num)

    elem_comp_CalT = {}
    if CalT_num != 0:
        for i in range(0, CalT_num):
            elem_comp_dict[i+total_elem_num] = [x[cand_cols_order['CalT']] for x in cand_elem]
            elem_comp_CalT[i] = np.array([x for x in elem_comp_dict[i+total_elem_num]]).reshape(-1, 1)
        elem_comp_all['CalT'] = np.hstack([elem_comp_CalT[i] for i in range(CalT_num)]).reshape(-1, CalT_num)

    elem_comp_ReacT = {}
    if ReacT_num != 0:
        for i in range(0, ReacT_num):
            elem_comp_dict[i+total_elem_num+CalT_num] = [x[cand_cols_order['ReacT']] for x in cand_elem]
            elem_comp_ReacT[i] = np.array([x for x in elem_comp_dict[i+total_elem_num+CalT_num]]).reshape(-1, 1)
        elem_comp_all['ReacT'] = np.hstack([elem_comp_ReacT[i] for i in range(ReacT_num)]).reshape(-1, ReacT_num)

    elem_comp_WHSV = {}
    if WHSV_num != 0:
        for i in range(0, WHSV_num):
            elem_comp_dict[i+total_elem_num+CalT_num+ReacT_num] = [x[cand_cols_order['WHSV']] for x in cand_elem]
            elem_comp_WHSV[i] = np.array([x for x in elem_comp_dict[i+total_elem_num+CalT_num+ReacT_num]]).reshape(-1, 1)
        elem_comp_all['WHSV'] = np.hstack([elem_comp_WHSV[i] for i in range(WHSV_num)]).reshape(-1, WHSV_num)

    elem_comp_H2_ppm = {}
    if H2_ppm_num != 0:
        for i in range(0, H2_ppm_num):
            elem_comp_dict[i+total_elem_exp_num] = [x[cand_cols_order['H2_ppm']] for x in cand_elem]
            elem_comp_H2_ppm[i] = np.array([x for x in elem_comp_dict[i+total_elem_exp_num]]).reshape(-1, 1)
        elem_comp_all['H2_ppm'] = np.hstack([elem_comp_H2_ppm[i] for i in range(H2_ppm_num)]).reshape(-1, H2_ppm_num)

    elem_comp_CO_ppm = {}
    if CO_ppm_num != 0:
        for i in range(0, CO_ppm_num):
            elem_comp_dict[i+total_elem_exp_num+H2_ppm_num] = [x[cand_cols_order['CO_ppm']] for x in cand_elem]
            elem_comp_CO_ppm[i] = np.array([x for x in elem_comp_dict[i+total_elem_exp_num+H2_ppm_num]]).reshape(-1, 1)
        elem_comp_all['CO_ppm'] = np.hstack([elem_comp_CO_ppm[i] for i in range(CO_ppm_num)]).reshape(-1, CO_ppm_num)

    elem_comp_NH3_ppm = {}
    if NH3_ppm_num != 0:
        for i in range(0, NH3_ppm_num):
            elem_comp_dict[i+total_elem_exp_num+H2_ppm_num+CO_ppm_num] = [x[cand_cols_order['NH3_ppm']] for x in cand_elem]
            elem_comp_NH3_ppm[i] = np.array([x for x in elem_comp_dict[i+total_elem_exp_num+H2_ppm_num+CO_ppm_num]]).reshape(-1, 1)
        elem_comp_all['NH3_ppm'] = np.hstack([elem_comp_NH3_ppm[i] for i in range(NH3_ppm_num)]).reshape(-1, NH3_ppm_num)

    elem_comp_C3H6_ppm = {}
    if C3H6_ppm_num != 0:
        for i in range(0, C3H6_ppm_num):
            elem_comp_dict[i+total_elem_exp_num+H2_ppm_num+CO_ppm_num+NH3_ppm_num] = [x[cand_cols_order['C3H6_ppm']] for x in cand_elem]
            elem_comp_C3H6_ppm[i] = np.array([x for x in elem_comp_dict[i+total_elem_exp_num+H2_ppm_num+CO_ppm_num+NH3_ppm_num]]).reshape(-1, 1)
        elem_comp_all['C3H6_ppm'] = np.hstack([elem_comp_C3H6_ppm[i] for i in range(C3H6_ppm_num)]).reshape(-1, C3H6_ppm_num)

    elem_comp_CH4_ppm = {}
    if CH4_ppm_num != 0:
        for i in range(0, CH4_ppm_num):
            elem_comp_dict[i+total_elem_exp_num+H2_ppm_num+CO_ppm_num+NH3_ppm_num+C3H6_ppm_num] = [x[cand_cols_order['CH4_ppm']] for x in cand_elem]
            elem_comp_CH4_ppm[i] = np.array([x for x in elem_comp_dict[i+total_elem_exp_num+H2_ppm_num+CO_ppm_num+NH3_ppm_num+C3H6_ppm_num]]).reshape(-1, 1)
        elem_comp_all['CH4_ppm'] = np.hstack([elem_comp_CH4_ppm[i] for i in range(CH4_ppm_num)]).reshape(-1, CH4_ppm_num)

    cand_comp = np.hstack([elem_comp_all[i] for i in list(cand_cols_order.keys())])
    return cand_comp


# 触媒組成＋実験条件を文字列(str)として整形し、棒(barh)グラフの軸要素(yticks)を作成する関数
def encode(elem_wt, elem, condition):
    cand_pgm_num, cand_prm_num = condition['cand_pgm_num'], condition['cand_prm_num']
    cand_add_num, cand_supp_num = condition['cand_add_num'], condition['cand_supp_num']
    CalT_num, ReacT_num, WHSV_num = condition['CalT_num'], condition['ReacT_num'], condition['WHSV_num']
    H2_ppm_num, CO_ppm_num = condition['H2_ppm_num'], condition['CO_ppm_num']
    NH3_ppm_num, C3H6_ppm_num = condition['NH3_ppm_num'], condition['C3H6_ppm_num']
    CH4_ppm_num = condition['CH4_ppm_num']

    total_elem_num = cand_pgm_num + cand_prm_num + cand_add_num + cand_supp_num
    total_cand_num = total_elem_num + CalT_num + ReacT_num + WHSV_num \
        + H2_ppm_num + CO_ppm_num + NH3_ppm_num + C3H6_ppm_num + CH4_ppm_num
    encode_name = ""
    for i in range(0, cand_pgm_num):
        if i == total_cand_num - 1:
            encode_name += f"{elem[i]} {elem_wt[i]}"
        else:
            encode_name += f"{elem[i]} {elem_wt[i]}, "

    for i in range(cand_pgm_num, cand_pgm_num+cand_prm_num):
        if i == total_cand_num - 1:
            encode_name += f"{elem[i]} {elem_wt[i]}"
        else:
            encode_name += f"{elem[i]} {elem_wt[i]}, "

    for i in range(cand_pgm_num+cand_prm_num, cand_pgm_num+cand_prm_num+cand_add_num):
        if i == total_cand_num - 1:
            encode_name += f"{elem[i]} {elem_wt[i]}"
        else:
            encode_name += f"{elem[i]} {elem_wt[i]}, "

    for i in range(cand_pgm_num+cand_prm_num+cand_add_num, total_cand_num):
        if i == total_cand_num - 1:
            encode_name += f"{elem[i]}"
        else:
            encode_name += f"{elem[i]}, "
    return encode_name


# 上で定義した2つの関数(make_cand_comp, encode)を順番に呼び出して処理を実行する関数。
def get_encode(idx, elem_wt, all_comb_slice, condition):
    catal_comb_idx = [all_comb_slice[i] for i in idx]
    catal_comp = make_cand_comp(catal_comb_idx, condition)
    catal_name = [encode(elem_wt, elem, condition) for elem in catal_comp]
    return catal_name


# 上で定義した get_encode関数で文字列(str)に変換した触媒組成を、再び要素ごとに分解してpd.DataFrame形式で出力する関数
def cand_separation(cand_data, condition):
    cand_pgm_num, cand_prm_num, cand_add_num = condition['cand_pgm_num'], condition['cand_prm_num'], condition['cand_add_num']
    CalT_num, ReacT_num, WHSV_num = condition['CalT_num'], condition['ReacT_num'], condition['WHSV_num']
    H2_ppm_num, CO_ppm_num = condition['H2_ppm_num'], condition['CO_ppm_num']
    NH3_ppm_num, C3H6_ppm_num = condition['NH3_ppm_num'], condition['C3H6_ppm_num']
    CH4_ppm_num = condition['CH4_ppm_num']

    data_cols = setting.data_columns(condition)
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
    if cand_prm_num != 0:
        for i in range(cand_prm_num):
            cand_sep_data_wt_prm = cand_sep_data[f'Promoter{i+1}_wt%'].str.strip(',')
            sep_data_wt.append(cand_sep_data_wt_prm)
    if cand_add_num != 0:
        for i in range(cand_add_num):
            cand_sep_data_wt_add = cand_sep_data[f'Ad{i+1}_wt%'].str.strip(',')
            sep_data_wt.append(cand_sep_data_wt_add)

    cand_wt = pd.DataFrame(sep_data_wt).transpose()
    cand_wt.columns = cand_wt_cols

    cand_elem = cand_sep_data.loc[:, cand_elem_cols]
    cand_ei = cand_sep_data.loc[:, 'ei']
    if CalT_num != 0 and ReacT_num == 0 and WHSV_num == 0:
        cand_CalT = cand_sep_data.loc[:, 'Calc T']
        cand = pd.concat([cand_elem, cand_wt, cand_CalT, cand_ei], axis=1)
    elif CalT_num != 0 and ReacT_num != 0 and WHSV_num != 0:
        cand_CalT = cand_sep_data.loc[:, 'Calc T']
        cand_RecT = cand_sep_data.loc[:, 'Reaction T']
        cand_WHSV = cand_sep_data.loc[:, 'WHSV (mL g-1 h-1)']
        cand = pd.concat([cand_elem, cand_wt, cand_CalT,
                        cand_RecT, cand_WHSV, cand_ei], axis=1)
    elif CalT_num != 0 and ReacT_num != 0 and WHSV_num == 0:
        cand_CalT = cand_sep_data.loc[:, 'Calc T']
        cand_RecT = cand_sep_data.loc[:, 'Reaction T']
        cand = pd.concat([cand_elem, cand_wt, cand_CalT,
                        cand_RecT, cand_ei], axis=1)
        cand_labels = cand_labels.remove('WHSV (mL g-1 h-1)')
    elif CalT_num != 0 and ReacT_num == 0 and WHSV_num != 0:
        cand_CalT = cand_sep_data.loc[:, 'Calc T']
        cand_WHSV = cand_sep_data.loc[:, 'WHSV (mL g-1 h-1)']
        cand = pd.concat([cand_elem, cand_wt, cand_CalT,
                        cand_WHSV, cand_ei], axis=1)
        cand_labels = cand_labels.remove('Reaction T')
    elif CalT_num == 0 and ReacT_num == 0 and WHSV_num == 0:
        if H2_ppm_num != 0 or CO_ppm_num != 0 or NH3_ppm_num != 0 or C3H6_ppm_num != 0 or CH4_ppm_num != 0:
            cand_H2_ppm = cand_sep_data.loc[:, 'H2_ppm']
            cand_CO_ppm = cand_sep_data.loc[:, 'CO_ppm']
            cand_NH3_ppm = cand_sep_data.loc[:, 'NH3_ppm']
            cand_C3H6_ppm = cand_sep_data.loc[:, 'C3H6_ppm']
            cand_CH4_ppm = cand_sep_data.loc[:, 'CH4_ppm']
            cand = pd.concat([cand_elem, cand_wt, cand_H2_ppm, cand_CO_ppm,
                            cand_NH3_ppm, cand_C3H6_ppm, cand_CH4_ppm, cand_ei], axis=1)
        else:
            cand = pd.concat([cand_elem, cand_wt, cand_ei], axis=1)
    cand = cand.reindex(columns=cand_labels)
    cand[cand_wt_cols] = cand[cand_wt_cols].astype(float)
    
    return cand


def data_str(data, condition, data_cols):
    pgm_num, prm_num = condition['pgm_num'], condition['prm_num']
    add_num, supp_num = condition['add_num'], condition['supp_num']
    CalT_num, ReacT_num, WHSV_num = condition['CalT_num'], condition['ReacT_num'], condition['WHSV_num']
    H2_ppm_num, CO_ppm_num = condition['H2_ppm_num'], condition['CO_ppm_num']
    NH3_ppm_num, C3H6_ppm_num = condition['NH3_ppm_num'], condition['C3H6_ppm_num']
    CH4_ppm_num = condition['CH4_ppm_num']
    elem_cols, wt_cols = data_cols['elem'], data_cols['wt']

    data['Catalyst'] = ''
    for i in range(len(data)):
        for j in range(pgm_num+prm_num+add_num):
            if data[wt_cols[j]][i] > 0.00:
                data['Catalyst'][i] += data[elem_cols[j]][i] + \
                    '(' + data[wt_cols[j]][i].astype(str) + '), '
        if supp_num != 0:
            data['Catalyst'][i] += data['Support_name'][i]
        if CalT_num == 1:
            data['Catalyst'][i] += ', ' + \
                '(' + data['Calc T'][i].astype(str) + ')'
        if ReacT_num == 1:
            data['Catalyst'][i] += ', ' + \
                '(' + data['Reaction T'][i].astype(str) + ')'
        if WHSV_num == 1:
            data['Catalyst'][i] += ', ' + \
                '(' + data['WHSV (mL g-1 h-1)'][i].astype(str) + ')'
        if H2_ppm_num == 1:
            data['Catalyst'][i] += ', ' + \
                '(' + data['H2_ppm'][i].astype(str) + ')'
        if CO_ppm_num == 1:
            data['Catalyst'][i] += ', ' + \
                '(' + data['CO_ppm'][i].astype(str) + ')'
        if NH3_ppm_num == 1:
            data['Catalyst'][i] += ', ' + \
                '(' + data['NH3_ppm'][i].astype(str) + ')'
        if C3H6_ppm_num == 1:
            data['Catalyst'][i] += ', ' + \
                '(' + data['C3H6_ppm'][i].astype(str) + ')'
        if CH4_ppm_num == 1:
            data['Catalyst'][i] += ', ' + \
                '(' + data['CH4_ppm'][i].astype(str) + ')'
    return data


# 上で定期した cand_separation関数で要素ごとに分解した触媒組成を、再び文字列(str)として結合し直す関数(K-cluster後plotのために使用)
def cand_str(cand, condition, data_cols):
    cand_pgm_num, cand_prm_num = condition['cand_pgm_num'], condition['cand_prm_num']
    cand_add_num, cand_supp_num = condition['cand_add_num'], condition['cand_supp_num']
    CalT_num, ReacT_num, WHSV_num = condition['CalT_num'], condition['ReacT_num'], condition['WHSV_num']
    H2_ppm_num, CO_ppm_num = condition['H2_ppm_num'], condition['CO_ppm_num']
    NH3_ppm_num, C3H6_ppm_num = condition['NH3_ppm_num'], condition['C3H6_ppm_num']
    CH4_ppm_num = condition['CH4_ppm_num']
    cand_elem_cols, cand_wt_cols = data_cols['cand_elem'], data_cols['cand_wt']

    cand['Top catal.'] = ''
    for i in range(len(cand)):
        for j in range(cand_pgm_num+cand_prm_num+cand_add_num):
            if cand[cand_wt_cols[j]][i] > 0.00:
                cand['Top catal.'][i] += cand[cand_elem_cols[j]][i] + \
                    '(' + cand[cand_wt_cols[j]][i].astype(str) + '), '
        if cand_supp_num != 0:
            cand['Top catal.'][i] += cand['Support_name'][i]
        if CalT_num == 1:
            cand['Top catal.'][i] += ' ' + '(' + cand['Calc T'][i] + ')'
        if ReacT_num == 1:
            cand['Top catal.'][i] += ', ' + '(' + cand['Reaction T'][i] + ')'
        if WHSV_num == 1:
            cand['Top catal.'][i] += ', ' + \
                '(' + cand['WHSV (mL g-1 h-1)'][i] + ')'
        if H2_ppm_num == 1:
            cand['Top catal.'][i] += ', ' + \
                '(' + cand['H2_ppm'][i] + ')'
        if CO_ppm_num == 1:
            cand['Top catal.'][i] += ', ' + \
                '(' + cand['CO_ppm'][i] + ')'
        if NH3_ppm_num == 1:
            cand['Top catal.'][i] += ', ' + \
                '(' + cand['NH3_ppm'][i] + ')'
        if C3H6_ppm_num == 1:
            cand['Top catal.'][i] += ', ' + \
                '(' + cand['C3H6_ppm'][i] + ')'
        if CH4_ppm_num == 1:
            cand['Top catal.'][i] += ', ' + \
                '(' + cand['CH4_ppm'][i] + ')'
    return cand