#!/usr/bin/env python
import numpy as np
import pandas as pd
import sys, os
sys.dont_write_bytecode = True
import warnings
warnings.filterwarnings('ignore')

def calc_condition():
    condition = {}

    ### Basic calculation conditions ###
    date = 20220318
    computer = 'I' # Select 'A' or 'B' or 'I'
    core = 36 # default: A==68, B==36, I==36
    hyper_threading = False # default: False
    Reaction = 'rwgs_250'
    test_cal = False # default: False
    data_sheet_name = 'test_data'
    train_data_sheet_name = 'train_data'
    ML_model = 'ETR' # Now, only support ETR
    
    ### Define ML model (0: conv, 1: prop1, 2: prop2) ###
    add_model = 2
    Search_method = 'all' # or 'local' (in the future...)
    
    ### Input catalyst shapes (Check your dataset.xlsx files) ###
    add_num = 5                         # 3, 4, 5

    ### Output condtions ###
    # 0. Basics
    save_depth = 100                   # if save_depth == -1: full depth
    catal_number = 100
    K_cluster = 100
    
    ### Define the inpot catalytic compositions ###

    # 1. Additives
    cand_add_num = 5
    add_wt = [0, 0.5, 1, 2, 3, 5, 7, 10]
    fix_add_num = 0                     # 0, 1, 2, ..., max_add_num-1
    essential_adds = [] # for example ['V', 'Ce', 'Nb', 'Rb', 'Mo', 'Ba', 'Re'] <- if fix_add_num == 0: this list has no meaning.

    ### Define target name ###
    target_name = 'CO formation rate_mmol min-1 gcat-1'
    
    if test_cal == True:
        add_wt = [1]
    else:
        pass

    condition = {
        'date': date, 'computer': computer, 'core': core, 'Reaction': Reaction, 'hyper_threading': hyper_threading,
        'data_sheet_name': data_sheet_name, 'train_data_sheet_name': train_data_sheet_name,
        'ML_model': ML_model, 'add_model': add_model, 'Search_method': Search_method,
        'cand_add_num': cand_add_num, 'add_num': add_num, 'add_wt': add_wt, 'fix_add_num': fix_add_num, 'essential_adds': essential_adds,
        'save_depth': save_depth, 'catal_number': catal_number, 'K_cluster': K_cluster, 'target_name': target_name
        }
    return condition

def use_desc():
    desc = {}
    add_use_desc = [
        # 'AN',
        'Group',
        # 'Period',
        'AW', ### Do not comment out !!! ###
        # 'Atom rad',
        # 'Ion rad',
        # 'Coval rad',
        # 'vdW rad',
        # 'Crys rad',
        'EN_Allred',
        # 'EN_Pauling',
        # 'EN_Mulliken',
        # 'EN_Allen',
        # 'EN_M&B',
        # 'EN_T&O',
        'm.p.',
        # 'b.p.',
        'd_fus_H',
        'Density',
        # 'a x 106 ',
        # 'Heat cap',
        # 'Therm Cond',
        # 'Ion E',
        # 'EA',
        # 'VE',
        # 'Surf E',
        # 'E_form_M',
        # 'WF_M',
        'Eads_CO2',
        # 'OCO_angle',
        # 'Eads_O-atom',
        # 'MW_Ox',
        # 'Ox Num_Oxide',
        # 'E_form_Ox',
        # 'Density_Ox',
        'BG_Oxide',
        # 'Ox Num_Precursor',
        # 'MW_Precursor',
        # 'Solubility',
        # 'Valency'
        ]
    desc = {'add': add_use_desc}
    
    return desc

def data_columns(condition):
    cand_add_num = condition['cand_add_num']
    data_cols = {}
    
    elem_cols = ['Ad1', 'Ad2', 'Ad3', 'Ad4', 'Ad5']
    wt_cols = [
        f'{elem_cols[0]}_wt%', f'{elem_cols[1]}_wt%', f'{elem_cols[2]}_wt%',
        f'{elem_cols[3]}_wt%', f'{elem_cols[4]}_wt%'
        ]
    data_labels = [
        'No.',
        'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%','Ad5', 'Ad5_wt%',
        'CO Yield_%','CO formation rate_mmol min-1 gcat-1',
        'Iteration',
        ]
    use_add_desc = [ # 'AW' is not used as Descriptor.
        'Group', 'AW', 'EN_Allred', 'm.p.', 'd_fus_H',
        'Density', 'BG_Oxide', 'Ox Num_Oxide',
        'Eads_CO2',
        ]

    if cand_add_num == 2:
        cand_elem_cols = ['Ad1', 'Ad2']
        cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%']
        cand_labels = ['Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'ei']
    elif cand_add_num == 3:
        cand_elem_cols = ['Ad1', 'Ad2', 'Ad3']
        cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%']
        cand_labels = ['Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'ei']
    elif cand_add_num == 4:
        cand_elem_cols = ['Ad1', 'Ad2', 'Ad3', 'Ad4']
        cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
        cand_labels = ['Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%', 'ei']
    elif cand_add_num == 5:
        cand_elem_cols = ['Ad1', 'Ad2', 'Ad3', 'Ad4', 'Ad5']
        cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
        cand_labels = ['Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%', 'Ad5', 'Ad5_wt%', 'ei']

    data_cols = {
        'elem': elem_cols, 'wt': wt_cols, 'labels': data_labels,
        'cand_elem': cand_elem_cols, 'cand_wt': cand_wt_cols, 'cand_labels': cand_labels,
        'add_desc': use_add_desc,
        }
    return data_cols

def desc_columns():
    desc_cols = {}
    basic_desc_columns = [
        'AN', 'Group', 'Period', 'AW',
        'Atom rad', 'Ion rad', 'Coval rad', 'vdW rad', 'Crys rad',
        'EN_Allred', 'EN_Pauling', 'EN_Mulliken', 'EN_Allen', 'EN_M&B', 'EN_T&O',
        'm.p.', 'b.p.', 'd_fus_H', 'Density', 'a x 106 ', 'Heat cap', 'Therm Cond',
        'Ion E', 'EA', 'VE', 'Surf E', 'E_form_M', 'WF_M',
        'Eads_CO2', 'OCO_angle', 'Eads_O-atom', 'MW_Ox',
        'Ox Num_Oxide', 'E_form_Ox', 'Density_Ox', 'BG_Oxide',
        'Ox Num_Precursor', 'MW_Precursor', 'Solubility', 'Valency'
        ]
    
    noble_gas = ['He', 'Ne', 'Ar', 'Kr', 'Xe']
    pgm_plus_ReAu = ['H', 'Ru', 'Rh', 'Pd', 'Re', 'Os', 'Ir', 'Pt', 'Au']
    
    # 'Re' is treated as Additives.
    drop_elems = [
        'Be', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'As', 'Se',
        'Br', 'Tc', 'I', 'Pm', 'Hg', 'Tl',
        'Ru', 'Rh', 'Pd','Os', 'Ir', 'Pt', 'Au'
        ]

    desc_cols = {
        'noble_gas': noble_gas, 'pgm_plus_ReAu': pgm_plus_ReAu, 'drop_elems': drop_elems,
        'basic_desc_columns': basic_desc_columns
        }

    return desc_cols
