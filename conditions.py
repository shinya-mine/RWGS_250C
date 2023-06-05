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
    computer = 'A' # Select 'A' or 'B' or 'I'
    core = 36 # default: A==68, B==36, I==36
    hyper_threading = False # default: False
    Reaction = 'rwgs_250'
    test_cal = True # default: False
    data_sheet_name = 'test_data'
    train_data_sheet_name = 'train_data'
    ML_model = 'ETR' # or 'RFR'
    
    ### Define ML model (0: conv, 1: prop1, 2: prop2) ###
    pgm_model = 2
    add_model = 2
    supp_model = 2
    Search_method = 'all' # or 'local' (in the future...)
    
    ### Input catalyst shapes (Check your dataset.xlsx files) ###
    pgm_num = 0                         # 0, 1, 2, 3
    add_num = 5                         # 3, 4, 5
    supp_num = 0                        # 0 or 1

    ### Output condtions ###
    # 0. Basics
    save_depth = 100                   # if save_depth == -1: full depth
    catal_number = 100
    K_cluster = 100
    # 0-1. Option of extractions.
    extractions = []
    
    ### Define the inpot catalytic compositions ###
    # 1. PGMs
    cand_pgm_num = 0
    pgm_wt = []
    fix_pgm_num = 0                     # # 0, 1, 2, ..., max_pgm_num-1
    max_pgm_wt = 0
    essential_pgms = []                # if fix_pgm_num == 0: this list has no meaning.
    essential_pgm_wt = []               # 'CH3OH':'Pt', 'N2O':'essential_pgm'
    specific_pgms = []
    specific_pgm_wt = []

    # 2. Additives
    cand_add_num = 5
    add_wt = [0, 0.2, 0.5, 0.7, 1.0, 1.2, 1.5, 2]
    fix_add_num = 0                     # 0, 1, 2, ..., max_add_num-1
    essential_adds = [] # ['V', 'Ce', 'Nb', 'Rb', 'Mo', 'Ba', 'Re']                 # if fix_add_num == 0: this list has no meaning.
    essential_add_wt = []
    specific_adds = []
    specific_add_wt = []

    # 3. Supports
    cand_supp_num = 0
    fix_supp_num = 0                    # 0, 1, 2, ..., max_supp_num
    essential_supp = []                 # if fix_supp_num == 0: this list has no meaning.(ex. 'Al2O3_gamma', ...)

    # 4. Calcination temp.
    if Reaction == 'CH3OH' or Reaction == 'NH3SCR':
        CalT_num = 1                    # 0(False) or 1(True)
        CalT_list = [400, 500, 600]     # if CalT_num == 0: this list has no meaning.
    else:
        CalT_num = 0
        CalT_list = []
    
    ### Define target name ###
    if Reaction == 'rwgs_250' or Reaction == 'rwgs_300' or Reaction == 'rwgs_250_1wt':
        target_name = 'CO formation rate_mmol min-1 gcat-1'
    elif Reaction == 'CH3OH':
        target_name = 'CH3OH rate (mmol g-1 h-1)'
    else:
        target_name = 'target'
        
    if Reaction == 'N2O' or Reaction == 'H2SCR' or Reaction == 'NH3SCR' or Reaction == 'CH4' or Reaction == 'EtOH':
        target_temp = 20
    else:
        target_temp = -1
    
    ### Other option (Skip rows) ###
    if Reaction == 'H2SCR' or Reaction == 'NH3SCR':
        skip_rows = 1
    else:
        skip_rows = 0
    
    if test_cal == True:
        pgm_wt = [0.1]
        essential_pgm_wt = [0.1]
        specific_pgm_wt = [0.1]
        add_wt = [1]
        essential_add_wt = [1]
        specific_add_wt = [1]
    else:
        pass

    condition = {
        'date': date, 'computer': computer, 'core': core, 'Reaction': Reaction, 'hyper_threading': hyper_threading,
        'data_sheet_name': data_sheet_name, 'train_data_sheet_name': train_data_sheet_name,
        'ML_model': ML_model, 'pgm_model': pgm_model, 'add_model': add_model, 'supp_model': supp_model, 'Search_method': Search_method,
        'cand_pgm_num': cand_pgm_num, 'cand_add_num': cand_add_num, 'cand_supp_num': cand_supp_num,
        'fix_pgm_num': fix_pgm_num, 'pgm_num': pgm_num, 'add_num': add_num, 'supp_num': supp_num,
        'max_pgm_wt': max_pgm_wt, 'essential_pgms': essential_pgms, 'essential_pgm_wt': essential_pgm_wt,
        'pgm_wt': pgm_wt, 'add_wt': add_wt, 'fix_add_num': fix_add_num, 'essential_adds': essential_adds, 'essential_add_wt': essential_add_wt,
        'specific_pgms': specific_pgms, 'specific_pgm_wt': specific_pgm_wt, 'specific_adds': specific_adds, 'specific_add_wt': specific_add_wt,
        'fix_supp_num': fix_supp_num, 'essential_supp': essential_supp, 'CalT_num': CalT_num, 'CalT_list': CalT_list,
        'save_depth': save_depth, 'catal_number': catal_number, 'K_cluster': K_cluster, 'extractions': extractions,
        'target_name': target_name, 'target_temp': target_temp, 'skip_rows': skip_rows,
        }
    return condition

def use_desc():
    desc = {}
    pgm_use_desc = [
        #'AN (PGM)',
        #'Group (PGM)',
        #'Period (PGM)',
        'AW (PGM)',                                ### Do not comment out !!!
        'Atom rad (PGM)',
        # 'Ion rad (PGM)',
        #'Coval rad (PGM)',
        #'vdW rad (PGM)',
        #'Crys rad (PGM)',
        #'EN_Allred (PGM)',
        #'EN_Pauling (PGM)',
        #'EN_Mulliken (PGM)',
        #'EN_Allen (PGM)',
        #'EN_M&B (PGM)',
        #'EN_T&O (PGM)',
        #'m.p. (PGM)',
        #'b.p. (PGM)',
        #'d_fus_H (PGM)',
        #'Density (PGM)',
        #'a x 106  (PGM)',
        #'Heat cap (PGM)',
        #'Therm Cond (PGM)',
        #'Ion E (PGM)',
        #'EA (PGM)',
        #'VE (PGM)',
        #'Surf E (PGM)',
        #'d_H_des_homo (PGM)',
        #'d_H_des_hetero (PGM)',
        #'FMM M (PGM)',
        #'E_form_M (PGM)',
        #'BG_M (PGM)',
        #'WF_M (PGM)',
        # 'Eads_CO2 (PGM)',
        #'OCO_angle (PGM)',
        #'Eads_O-atom (PGM)',
        #'Eads_C-atom (PGM)',
        #'Eads_H-atom (PGM)',
        #'Eads_N2O (PGM)',
        #'NNO_angle (PGM)',
        #'Eads_N-atom (PGM)',
        #'MW_Ox (PGM)',
        #'Ox Num_Oxide (PGM)',
        #'FMM_Ox (PGM)',
        #'E_form_Ox (PGM)',
        #'Density_Ox (PGM)',
        #'BG_Oxide (PGM)',
        #'Surf E Ox (PGM)',
        #'WF_Ox (PGM)',
        #'Ox Num_Precursor (PGM)',
        #'Valency (PGM)'
        ]

    add_use_desc = [
        #'AN (Additive)',
        'Group (Additive)',
        #'Period (Additive)',
        'AW (Additive)',                                   ### Do not comment out !!!
        #'Atom rad (Additive)',
        #'Ion rad (Additive)',
        #'Coval rad (Additive)',
        #'vdW rad (Additive)',
        #'Crys rad (Additive)',
        'EN_Allred (Additive)',
        #'EN_Pauling (Additive)',
        #'EN_Mulliken (Additive)',
        #'EN_Allen (Additive)',
        #'EN_M&B (Additive)',
        #'EN_T&O (Additive)',
        'm.p. (Additive)',
        #'b.p. (Additive)',
        'd_fus_H (Additive)',
        'Density (Additive)',
        #'a x 106  (Additive)',
        #'Heat cap (Additive)',
        #'Therm Cond (Additive)',
        #'Ion E (Additive)',
        #'EA (Additive)',
        #'VE (Additive)',
        #'Surf E (Additive)',
        #'d_H_des_homo (Additive)',
        #'d_H_des_hetero (Additive)',
        #'FMM M (Additive)',
        #'E_form_M (Additive)',
        #'BG_M (Additive)',
        #'WF_M (Additive)',
        #'Eads_CO2 (Additive)',
        #'OCO_angle (Additive)',
        #'Eads_O-atom (Additive)',
        #'Eads_C-atom (Additive)',
        #'Eads_H-atom (Additive)',
        #'Eads_N2O (Additive)',
        #'NNO_angle (Additive)',
        #'Eads_N-atom (Additive)',
        #'MW_Ox (Additive)',
        #'Ox Num_Oxide (Additive)',
        #'FMM_Ox (Additive)',
        #'E_form_Ox (Additive)',
        #'Density_Ox (Additive)',
        #'BG_Oxide (Additive)',
        #'Surf E Ox (Additive)',
        #'WF_Ox (Additive)',
        #'Ox Num_Precursor (Additive)',
        #'Valency (Additives)'
        ]

    supp_use_desc = [ # index_col = 'Support_name'
        'MW (Support)', ### Do not comment out !!!
        #'Evac (Support)',
        'Bulk_EN (Support)',
        #'E_form (Support)',
        #'Density (Support)',
        #'Dielect Prop (Support)',
        #'Elasticity (Support)',
        #'BG (Support)',
        #'Surf E (Support)',
        #'Ion Pot (Support)',
        #'EA (Support)',
        #'WF (Support)',
        # 'SSA (Support)',
        #'SSA cal at 700C (Support)',
        #'SiO2/Al2O3 ratio (Support)',
        #'Pore size (Support)',
        #'Cation (Support)'
    ]

    desc = {'pgm': pgm_use_desc, 'add': add_use_desc, 'supp': supp_use_desc}
    
    return desc

def data_columns(Reaction, condition):
    cand_pgm_num, cand_add_num, cand_supp_num, CalT_num = condition['cand_pgm_num'], condition['cand_add_num'], condition['cand_supp_num'], condition['CalT_num']
    data_cols = {}
    
    if Reaction == 'rwgs_250':
        elem_cols = ['Ad1', 'Ad2', 'Ad3', 'Ad4', 'Ad5']
        wt_cols = [
            f'{elem_cols[0]}_wt%', f'{elem_cols[1]}_wt%', f'{elem_cols[2]}_wt%',
            f'{elem_cols[3]}_wt%', f'{elem_cols[4]}_wt%'
            ]
        data_labels = [
            'No.',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%','Ad5', 'Ad5_wt%', #'Support_name',
            'CO Yield_%','CO formation rate_mmol min-1 gcat-1',
            'Iteration', 'Catal prep', 'Reaction', 'Note', 'Experimantal condition'
            ]
        use_pgm_desc = ['AW'] # 'AW' is not used as Descriptor.
        use_add_desc = [ # 'AW' is not used as Descriptor.
            'Group', 'AW', 'EN_Allred', 'm.p.', 'd_fus_H',
            'Density', 'BG_Oxide', 'Ox Num_Oxide',
            'Eads_CO2',
            ]
        use_supp_desc = ['MW (Support)'] # 'MW (Support)' is not used as Descriptor.

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
    
    elif Reaction == 'rwgs_300':
        elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Ad3', 'Ad4']
        wt_cols = [
            f'{elem_cols[0]}_wt%', f'{elem_cols[1]}_wt%', f'{elem_cols[2]}_wt%',
            f'{elem_cols[3]}_wt%', f'{elem_cols[4]}_wt%', f'{elem_cols[5]}_wt%'
            ]
        data_labels = [
            'No.',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2','Ad2_wt%', 'Ad3', 'Ad3_wt%',
            'Support_name', 'Calc. temp. (℃)',
            'CO Yield_%', 'CO formation rate_mmol min-1 gcat-1','CO selec',
            'Iteration', 'Catal prep', 'Reaction', 'Note', 'Experimantal condition'
            ]
        use_pgm_desc = ['AW (PGM)']
        use_add_desc = ['AW (Additive)']
        use_supp_desc = ['MW (Support)']
        
        if cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]

    elif Reaction == 'rwgs_250_1wt':
        elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Ad4', 'Ad5']
        wt_cols = [
            f'{elem_cols[0]}_wt%', f'{elem_cols[1]}_wt%', f'{elem_cols[2]}_wt%',
            f'{elem_cols[3]}_wt%', f'{elem_cols[4]}_wt%', f'{elem_cols[5]}_wt%', f'{elem_cols[6]}_wt%', f'{elem_cols[7]}_wt%'
            ]
        data_labels = [
            'No.', 'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
            'Ad5', 'Ad5_wt%', 'Support_name', 'Calc. temp. (℃)', 'CO Yield_%',
            'CO selec', 'CO formation rate_mmol min-1 gcat-1', 'Iteration',
            'Catal prep', 'Reaction', 'Date', 'Note'
            ]
        use_pgm_desc = ['AW (PGM)']
        use_add_desc = ['AW (Additive)']
        use_supp_desc = ['MW (Support)']

        if cand_pgm_num == 3 and cand_add_num == 5 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Ad4', 'Ad5', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%', f'{cand_elem_cols[6]}_wt%', f'{cand_elem_cols[7]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%', 'Ad5', 'Ad5_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 4 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Ad4', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%', f'{cand_elem_cols[6]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 5 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Ad3', 'Ad4', 'Ad5', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%', f'{cand_elem_cols[6]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%', 'Ad5', 'Ad5_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 4 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Ad3', 'Ad4', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 5 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Ad3', 'Ad4', 'Ad5', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%', 'Ad5', 'Ad5_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 4 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Ad3', 'Ad4', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        
    elif Reaction == 'CH3OH':
        elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Ad4']
        wt_cols = [
            f'{elem_cols[0]}_wt%', f'{elem_cols[1]}_wt%', f'{elem_cols[2]}_wt%',
            f'{elem_cols[3]}_wt%', f'{elem_cols[4]}_wt%', f'{elem_cols[5]}_wt%',
            f'{elem_cols[6]}_wt%'
            ]
        data_labels = [
            'No.', 'Cat amount/g',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
            'Support_name', 'Calc. temp. (℃)',
            'CH3OH rate (mmol g-1 h-1)', 'CH3OH yield (%)',
            'CH3OH selec (%)', 'CH3OH rate x selec',
            'CO rate (mmol g-1 h-1)', 'CO yield (%)', 'CO selec (%)',
            'CH4 rate (mmol g-1 h-1)', 'CH4 yield (%)', 'CH4 selec (%)',
            'CH3OH-CH4', # rate (mmol g-1 h-1)
            'Iteration', 'Experimenter',
            'CH3OH (mmol)', 'CO (mmol)', 'CH4 (mmol)',
            'Date', 'Preparation', 'Reaction', 'Note'
            ]
        use_pgm_desc = ['AW (PGM)']
        use_add_desc = ['AW (Additive)']
        use_supp_desc = ['MW (Support)']

        if cand_pgm_num == 3 and cand_add_num == 4 and cand_supp_num == 1 and CalT_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Ad4', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%', f'{cand_elem_cols[6]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
                'Support', 'Calc. temp. (℃)', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 1 and CalT_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'Calc. temp. (℃)', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 1 and CalT_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'Calc. temp. (℃)', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 1 and cand_supp_num == 1 and CalT_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'Calc. temp. (℃)', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 4 and cand_supp_num == 1 and CalT_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Ad3', 'Ad4', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
                'Support', 'Calc. temp. (℃)', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 4 and cand_supp_num == 1 and CalT_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Ad3', 'Ad4', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
                'Support', 'Calc. temp. (℃)', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 1 and CalT_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'Calc. temp. (℃)', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 1 and CalT_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'Calc. temp. (℃)', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 1 and CalT_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'Calc. temp. (℃)', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 1 and CalT_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'Calc. temp. (℃)', 'ei'
                ]
        
    elif Reaction == 'N2O':
        elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3']
        wt_cols = [
            f'{elem_cols[0]}_wt%', f'{elem_cols[1]}_wt%', f'{elem_cols[2]}_wt%',
            f'{elem_cols[3]}_wt%', f'{elem_cols[4]}_wt%', f'{elem_cols[5]}_wt%'
            ]
        data_labels = [
            'No.',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
            'Support_name', 'Calc. temp. (℃)',
            'T20 (℃)', 'T30 (℃)', 'T50 (℃)', 'T80 (℃)',
            'Conv at 500℃ (during pretreatment)', 'Cat. amount (mg)', 'Iteration',
            'Prep.', 'Reaction', 'Comment', 'Date'
            ]
        use_pgm_desc = ['AW (PGM)']
        use_add_desc = ['AW (Additive)']
        use_supp_desc = ['MW (Support)']
        
        if cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]

    elif Reaction == 'H2SCR':
        elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3']
        wt_cols = [
            f'{elem_cols[0]}_wt%', f'{elem_cols[1]}_wt%', f'{elem_cols[2]}_wt%',
            f'{elem_cols[3]}_wt%', f'{elem_cols[4]}_wt%', f'{elem_cols[5]}_wt%'
            ]
        data_labels = [
            'No.',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2','Ad2_wt%', 'Ad3', 'Ad3_wt%',
            'Support_name', 'Calc. temp. (℃)',
            'N2 yield_T10', 'N2 yield_T15', 'N2 yield_T20', 'N2 yield_T30', 'N2 yield_T50', 'N2 yield_T80', # target_cols
            'NO conv_T50', 'H2 conv_T50', 'NO conv_T80', 'H2 conv_T80'
            'NO conv [300℃] (%)','H2 conv [300℃] (%)', 'N2O yield [300℃] (%)',
            'N2 yield [300℃] (%)','NH3 yield [300℃] (%)', 'NO2 yield [300℃] (%)', 'N2 selec [300℃] (%)',
            'Cat. amount (mg)', 'Iteration', 'Preparation', 'Reaction', 'Date'
            ]
        use_pgm_desc = ['AW (PGM)']
        use_add_desc = ['AW (Additive)']
        use_supp_desc = ['MW (Support)']

        if cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]

    elif Reaction == 'NH3SCR':
        elem_cols = ['Ad1', 'Ad2', 'Ad3', 'Ad4']
        wt_cols = [
            f'{elem_cols[0]}_wt%', f'{elem_cols[1]}_wt%', f'{elem_cols[2]}_wt%',
            f'{elem_cols[3]}_wt%',
            ]
        data_labels = [
            'No.',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4','Ad4_wt%',
            'Support_name', 'Calc. temp. (℃)',
            'T20_NO conv', 'T30_NO conv', 'T50_NO conv',
            'T20_N2 yield', 'T30_N2 yield', 'T50_N2 yield',
            'T20_NH3 conv', 'T30_NH3 conv', 'T50_NH3 conv',
            'T20_N2 yield_2','T30_N2 yield_2', 'T50_N2 yield_2',
            'Cat. amount (mg)', 'Iteration',
            'B acid 200 C', 'L acid 200 C', 'Preparation', 'Reaction', 'Date'
            ]
        use_pgm_desc = []
        use_add_desc = ['AW']
        use_supp_desc = ['MW (Support)']

        if cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['Ad1', 'Ad2']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%']
            cand_labels = ['Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Support', 'Calc. temp. (℃)', 'ei']
        elif cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['Ad1', 'Ad2', 'Ad3']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%']
            cand_labels = ['Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Support', 'Calc. temp. (℃)', 'ei']
        elif cand_add_num == 4 and cand_supp_num == 1:
            cand_elem_cols = ['Ad1', 'Ad2', 'Ad3', 'Ad4']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = ['Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%', 'Support', 'Calc. temp. (℃)', 'ei']
        elif cand_add_num == 5 and cand_supp_num == 1:
            cand_elem_cols = ['Ad1', 'Ad2', 'Ad3', 'Ad4', 'Ad5']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = ['Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%', 'Ad5', 'Ad5_wt%', 'Support', 'Calc. temp. (℃)', 'ei']

    elif Reaction == 'CH4':
        elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3']
        wt_cols = [
            f'{elem_cols[0]}_wt%', f'{elem_cols[1]}_wt%', f'{elem_cols[2]}_wt%',
            f'{elem_cols[3]}_wt%', f'{elem_cols[4]}_wt%', f'{elem_cols[5]}_wt%'
            ]
        data_labels = [
            'No.',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
            'Support_name', 'Calc. temp. (℃)',
            'T20 (℃)', 'T30  (℃)', 'T50  (℃)',
            'Cat. amount (mg)', 'Iteration', 'Preparation','Reaction', 'Date'
            ]
        use_pgm_desc = ['AW (PGM)']
        use_add_desc = ['AW (Additive)']
        use_supp_desc = ['MW (Support)']

        if cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]
    
    elif Reaction == 'EtOH_CO2' and Reaction == 'EtOH_CO':
        elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Ad4']
        wt_cols = [
            f'{elem_cols[0]}_wt%', f'{elem_cols[1]}_wt%', f'{elem_cols[2]}_wt%',
            f'{elem_cols[3]}_wt%', f'{elem_cols[4]}_wt%', f'{elem_cols[5]}_wt%', f'{elem_cols[6]}_wt%'
            ]
        data_labels =[
            'Exp No.', 'Cat No',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
            'Support', 'Calc. temp. (℃)', 'Reac. temp. (℃)', 'GHSV', 'Pressure',
            'C3+OH rate (mmol g-1 h-1)', 'C3+OH yield (%)',
            'C2H5OH rate (mmol g-1 h-1)', 'C2H5OH yield (%)',
            'CH3OH rate (mmol g-1 h-1)', 'CH3OH yield (%)',
            'CO rate (mmol g-1 h-1)', 'CO yield (%)', 'C2-C4 rate (mmol g-1 h-1)',
            'C2-C4 yield (%)', 'CH4 rate (mmol g-1 h-1)', 'CH4 yield (%)',
            'CO2 conv (%)', 'Iteration',
            'Experimenter', 'Date', 'Note'
            ]
        use_pgm_desc = ['AW (PGM)']
        use_add_desc = ['AW (Additive)']
        use_supp_desc = ['MW (Support)']

        if cand_pgm_num == 3 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%', f'{cand_elem_cols[5]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 3 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Ad3', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%', f'{cand_elem_cols[4]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%', f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 3 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'PGM3', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%',
                            f'{cand_elem_cols[3]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 1 and cand_add_num == 2 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'Ad1', 'Ad2', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%',
                            f'{cand_elem_cols[1]}_wt%', f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%',
                'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%',
                'Support', 'ei'
                ]
        elif cand_pgm_num == 2 and cand_add_num == 1 and cand_supp_num == 1:
            cand_elem_cols = ['PGM1', 'PGM2', 'Ad1', 'Support']
            cand_wt_cols = [f'{cand_elem_cols[0]}_wt%', f'{cand_elem_cols[1]}_wt%',
                            f'{cand_elem_cols[2]}_wt%']
            cand_labels = [
                'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
                'Ad1', 'Ad1_wt%',
                'Support', 'ei'
                ]

    data_cols = {
        'elem': elem_cols, 'wt': wt_cols, 'labels': data_labels,
        'cand_elem': cand_elem_cols, 'cand_wt': cand_wt_cols, 'cand_labels': cand_labels,
        'pgm_desc': use_pgm_desc, 'add_desc': use_add_desc, 'supp_desc': use_supp_desc
        }
    return data_cols

def desc_columns(Reaction):
    desc_cols = {}
    all_desc_columns = [
        'AN', 'Name', 'Pseudopotential', 'Precursor', 'Group', 'Period', 'AW', 'Atom rad', 'Ion rad', 'Coval rad', 'vdW rad', 'Crys rad',
        'EN_Allred', 'EN_Pauling', 'EN_Mulliken', 'EN_Allen', 'EN_M&B', 'EN_T&O', 'm.p.', 'b.p.', 'd_fus_H', 'Density', 'a x 106 ', 'Heat cap',
        'Therm Cond', 'Ion E', 'EA', 'VE', 'Surf E', 'Homonuclear bonds', 'd_H_des_homo', 'Heteronuclear bonds', 'd_H_des_hetero', 'ID_M', 'Space group M',
        'Im-3m_M', 'P63/mmc_M', 'Fm-3m_M', 'Fd-3m_M', 'R-3m_M', 'Cmce_M', 'P2/c_M', 'P3121_M', 'Crys Sys M', 'Mag Ord M', 'Surf M', 'FMM M', 'E_form_M',
        'BG_M', 'WF_M', 'Eads_CO2', 'OCO_angle', 'Eads_O-atom', 'Eads_C-atom', 'Eads_N2O', 'NNO_angle', 'Eads_N-atom', 'ID_Ox', 'Composition_Ox', 'MW_Ox',
        'Ox Num_Oxide', 'Space group Ox', 'Group No', 'Crys Sys Ox', 'Mag Order Ox', 'FMM_Ox', 'E_form_Ox', 'Density_Ox', 'BG_Oxide', 'Stable surf', 'Surf E Ox',
        'WF_Ox', 'Ox Num_Precursor', 'Counter ion', 'H^+', 'NH4^+', 'OH^-', 'NO2^-', 'NO3^-', 'O^2-', 'CO3^2-', 'F^-', 'Cl^-', 'HC2O4^-', '(OCH2CH3)^-',
        '(OOCCH3)^-', '(C2O4NH4)^-', 'CAS No', 'MW_Precursor', 'Solubility', 'Supplier', 'Valency',
        'K 1s', 'L 2s', 'L 2p', 'M 3s', 'M 3p', 'M 3d', 'N 4s', 'N 4p', 'N 4d', 'N 4f', 'O 5s', 'O 5p', 'O 5d', 'O 5f', 'P 6s', 'P 6p', 'P 6d', 'K 1s pot',
        'L 2s pot', 'L 2p pot', 'M 3s pot', 'M 3p pot', 'M 3d pot', 'N 4s pot', 'N 4p pot', 'N 4d pot', 'N 4f pot', 'O 5s pot', 'O 5p pot',
        'O 5d pot', 'O 5f pot', 'P 6s pot', 'P 6p pot', 'P 6d pot'
        ]

    supp_desc_columns = [ # index_col = 'Support_name'
        'Compound', 'Structure', 'Surface (Support)', 'MW (Support)', 'Evac (Support)', 'Bulk_EN (Support)', 'E_form (Support)', 'Density (Support)',
        'Dielect Prop (Support)', 'Elasticity (Support)', 'BG (Support)', 'Surf E (Support)', 'Ion Pot (Support)', 'EA (Support)', 'WF (Support)',
        'SSA (Support)', 'SSA cal at 700C (Support)', 'SiO2/Al2O3 ratio (Support)', 'Pore size (Support)', 'Cation (Support)',
        'rwgs_250', 'rwgs_250_1st', 'rwgs_300', 'CH3OH', 'EtOH_CO2', 'EtOH_CO', 'H2SCR', 'NH3SCR', 'N2O', 'CH4',
        'Source','Sample name', 'Note', 'Ref'
        ]

    basic_desc_columns = [
        'AN', 'Group', 'Period', 'AW',
        'Atom rad', 'Ion rad', 'Coval rad', 'vdW rad', 'Crys rad',
        'EN_Allred', 'EN_Pauling', 'EN_Mulliken', 'EN_Allen', 'EN_M&B', 'EN_T&O',
        'm.p.', 'b.p.', 'd_fus_H', 'Density', 'a x 106 ', 'Heat cap', 'Therm Cond', 'Ion E', 'EA', 'VE', 'Surf E',
        'd_H_des_homo', 'd_H_des_hetero',
        'FMM M', 'E_form_M', 'BG_M', 'WF_M', 'Eads_CO2', 'OCO_angle', 'Eads_O-atom', 'Eads_C-atom', 'Eads_N2O', 'NNO_angle', 'Eads_N-atom',
        'MW_Ox', 'Ox Num_Oxide','FMM_Ox', 'E_form_Ox', 'Density_Ox', 'BG_Oxide',
        'Surf E Ox', 'WF_Ox', 'Ox Num_Precursor','Valency'
    ]

    Metal_space_group = [
        'Im-3m_M', 'P63/mmc_M', 'Fm-3m_M', 'Fd-3m_M', 'R-3m_M', 'Cmce_M', 'P2/c_M', 'P3121_M'
    ]

    Counter_ions = [
        'H^+', 'NH4^+', 'OH^-', 'NO2^-', 'NO3^-', 'O^2-', 'CO3^2-', 'F^-', 'Cl^-',
        'HC2O4^-', '(OCH2CH3)^-', '(OOCCH3)^-', '(C2O4NH4)^-'
    ]

    electron_config = [
        'K 1s', 'L 2s', 'L 2p', 'M 3s', 'M 3p', 'M 3d', 'N 4s', 'N 4p', 'N 4d', 'N 4f',
        'O 5s', 'O 5p', 'O 5d', 'O 5f', 'P 6s', 'P 6p', 'P 6d'
    ]

    potcar_config = [
        'K 1s pot', 'L 2s pot', 'L 2p pot', 'M 3s pot', 'M 3p pot', 'M 3d pot', 'N 4s pot','N 4p pot',
        'N 4d pot', 'N 4f pot', 'O 5s pot', 'O 5p pot', 'O 5d pot', 'O 5f pot', 'P 6s pot', 'P 6p pot', 'P 6d pot'
    ]
    
    noble_gas = ['He', 'Ne', 'Ar', 'Kr', 'Xe']
    pgm_plus_ReAu = ['H', 'Ru', 'Rh', 'Pd', 'Re', 'Os', 'Ir', 'Pt', 'Au']
    
    if Reaction == 'rwgs_250' or Reaction == 'rwgs_300' or Reaction == 'rwgs_250_1wt':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'As', 'Se',
            'Br', 'Tc', 'I', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd','Os', 'Ir', 'Pt', 'Au'
            ]

    elif Reaction == 'CH3OH':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au',
            ]

    elif Reaction == 'N2O':
        # 'Re' is treated as PGMs.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Re', 'Os', 'Ir', 'Pt', 'Au'
            ]

    elif Reaction == 'H2SCR':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au'
            ]

    elif Reaction == 'NH3SCR':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au'
            ]

    elif Reaction == 'CH4':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'As', 'Se',
            'Br', 'Tc', 'I', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au'
            ]

    elif Reaction == 'EtOH':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au'
            ]

    desc_cols = {
        'noble_gas': noble_gas, 'pgm_plus_ReAu': pgm_plus_ReAu, 'drop_elems': drop_elems, 'all_desc_columns': all_desc_columns,
        'basic_desc_columns': basic_desc_columns, 'supp_desc_columns': supp_desc_columns, 'Metal_space_group': Metal_space_group,
        'Counter_ions': Counter_ions, 'electron_config': electron_config, 'potcar_config': potcar_config,
        }

    return desc_cols
