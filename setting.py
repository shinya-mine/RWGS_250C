#!/usr/bin/env python
import pandas as pd
import sys
sys.dont_write_bytecode = True


def make_cols(elem_name='elem', elem_num=0):
    elem_cols = []
    for i in range(elem_num):
        elem_cols.append(f'{elem_name}{i+1}')
    wt_cols = list(map(lambda s: s + '_wt%', elem_cols))
    return elem_cols, wt_cols


def make_cand_labels(elem_name='cand_elem', cand_elem_num=0):
    cand_labels = []
    for i in range(cand_elem_num):
        cand_labels.append(f'{elem_name}{i+1}')
        cand_labels.append(f'{elem_name}{i+1}_wt%')
    return cand_labels


def make_list(condition):
    pgm_cols, pgm_wt_cols = make_cols(
        elem_name='PGM', elem_num=condition['pgm_num'])
    cand_pgm_cols, cand_pgm_wt_cols = make_cols(
        elem_name='PGM', elem_num=condition['cand_pgm_num'])

    prm_cols, prm_wt_cols = make_cols(
        elem_name='Promoter', elem_num=condition['prm_num'])
    cand_prm_cols, cand_prm_wt_cols = make_cols(
        elem_name='Promoter', elem_num=condition['cand_prm_num'])

    add_cols, add_wt_cols = make_cols(
        elem_name='Ad', elem_num=condition['add_num'])
    cand_add_cols, cand_add_wt_cols = make_cols(
        elem_name='Ad', elem_num=condition['cand_add_num'])

    elem_cols = pgm_cols + prm_cols + add_cols
    if condition['supp_num'] == 1:
        cand_elem_cols = cand_pgm_cols + \
            cand_prm_cols + cand_add_cols + ['Support_name']
    else:
        cand_elem_cols = cand_pgm_cols + cand_prm_cols + cand_add_cols
    wt_cols = pgm_wt_cols + prm_wt_cols + add_wt_cols
    cand_wt_cols = cand_pgm_wt_cols + cand_prm_wt_cols + cand_add_wt_cols
    gas_cols = ['H2_ppm', 'CO_ppm', 'NH3_ppm', 'C3H6_ppm', 'CH4_ppm']

    return elem_cols, cand_elem_cols, wt_cols, cand_wt_cols, gas_cols


def make_labels(condition):
    date, Reaction = condition['date'], condition['Reaction']
    skip_rows = condition['skip_rows']
    data_file_name = f'data/{date}_{Reaction}_dataset.xlsx'

    data = pd.read_excel(data_file_name, skiprows=skip_rows)
    data_labels = data.columns
    data_labels = [s for s in data_labels if 'Unnamed:' not in s]

    pgm_cand_labels = make_cand_labels(elem_name='PGM', cand_elem_num=condition['cand_pgm_num'])
    prm_cand_labels = make_cand_labels(elem_name='Promoter', cand_elem_num=condition['cand_prm_num'])
    add_cand_labels = make_cand_labels(elem_name='Ad', cand_elem_num=condition['cand_add_num'])

    cand_labels = pgm_cand_labels + prm_cand_labels + add_cand_labels

    if condition['cand_supp_num'] != 0:
        cand_labels += ['Support_name']
    if condition['CalT_num'] == 1:
        cand_labels += ['Calc T']
    if condition['ReacT_num'] == 1:
        cand_labels += ['Reaction T']
    if condition['WHSV_num'] == 1:
        cand_labels += ['WHSV (mL g-1 h-1)']
    if condition['H2_ppm_num'] == 1:
        cand_labels += ['H2_ppm']
    if condition['CO_ppm_num'] == 1:
        cand_labels += ['CO_ppm']
    if condition['NH3_ppm_num'] == 1:
        cand_labels += ['NH3_ppm']
    if condition['C3H6_ppm_num'] == 1:
        cand_labels += ['C3H6_ppm']
    if condition['CH4_ppm_num'] == 1:
        cand_labels += ['CH4_ppm']
    cand_labels += ['ei']
    return data_labels, cand_labels


def data_columns(condition):
    data_cols = {}

    elem_cols, cand_elem_cols, wt_cols, cand_wt_cols, gas_cols = make_list(condition)
    data_labels, cand_labels = make_labels(condition)

    data_cols = {
        'elem': elem_cols, 'wt': wt_cols, 'labels': data_labels,
        'cand_elem': cand_elem_cols, 'cand_wt': cand_wt_cols, 'cand_labels': cand_labels,
        'gas_cols': gas_cols,
    }
    return data_cols


def desc_columns(Reaction):
    desc_cols = {}

    elem_desc_columns = [
        'AN', 'Symbol', 'Name', 'Group', 'Period', 'AW', 'Atom rad',
        'Atom rad emp', 'Ion rad', 'Coval rad', 'Coval rad 3', 'vdW rad',
        'Crys rad', 'Metal bond rad', 'EN_Allred', 'EN_Pauling', 'EN_Mulliken',
        'EN_Allen', 'EN_M&B', 'EN_T&O', 'm.p.', 'b.p.', 'd_fus_H', 'Density',
        'Thermal expansion', 'Heat capacity', 'Thermal conductivity',
        'Ionization energy', 'Electron affinity', 'Valence electron',
        'Valence electron VASP', 'E_surf', 'Homonuclear bonds',
        'd_desoc_H_homo', 'Heteronuclear bonds', 'd_desoc_H_hetero', 'FMM Metal',
        'E_form Metal', 'BG Metal', 'Work function', 'Eads_CO2', 'OCO_angle',
        'Eads_H', 'Eads_H2', 'Eads_C', 'Eads_N', 'Eads_N2', 'Eads_O', 'Eads_O2',
        'Composition_Ox', 'MW_Oxide', 'Stable oxidation number', 'FMM Oxide',
        'E_form Oxide', 'Density Oxide', 'BG Oxide', 'Precursor',
        'Precursor oxidation number', 'MW_Precursor', 'Solubility', 'CPK color code'
    ]

    supp_desc_columns = [
        'Support_name', 'Compound', 'Structure', 'Surface', 'ID', 'Space_Group',
        'MW', 'EN_bulk', 'E_form', 'BG_PBEsol', 'E_surf', 'Ionization potential',
        'Electron affinity', 'Work function', 'Evac',
        'BG MP', 'E_form MP', 'Magnetization', 'Density', 'Poisson Ratio',
        'Elastic_Prop_Bulk_Voigt', 'm.p.', 'b.p.', 'd_form_H_298K', 'std_mol_S',
        'std_mol_Cp', 'Eads_H2', 'Eads_CO2', 'Eads_CO', 'Eads_NO', 'Eads_N2O',
        'Eads_CH4', 'Eads_CH3OH', 'Eads_C2H5OH', 'Edisads_HCOOH',
        'Edisads_CH3COOH', 'SSA specification', 'SSA w/o cal',
        'SSA cal at 700℃', 'N2 uptake at P/P0=0.015',
        'N2 uptake at P/P0=0.015 at 700℃', 'IR L acid intensity',
        'IR L acid cm-1', 'IR B acid intensity', 'IR B acid cm-1',
        'IR L acid intensity cal at 700℃', 'IR L acid cm-1 cal at 700℃',
        'IR B acid intensity cal at 700℃', 'IR B acid cm-1 cal at 700℃',
        'Si/Al ratio', 'Pore size', 'Cation', 'rwgs_250', 'rwgs_250_1wt',
        'rwgs_300', 'rwgs_TL', 'CH3OH', 'CH3OH_flow', 'EtOH_CO2', 'EtOH_CO',
        'H2SCR', 'NH3SCR', 'N2O', 'N2O-SCR', 'CH4', 'Support_name_python',
        'Source', 'Sample name', 'Note', 'Note (detail)'
    ]
    
    # https://matplotlib.org/stable/tutorials/text/mathtext.html
    elem_desc_dict = {
        'd_fus_H': '$\Delta_{fus}$$\it{H}$$^{\ominus}$',
        'E_surf': '$\it{E}$$_{surf}$',
        'd_desoc_H_homo': '$\Delta_{Dissoc}$$\it{H}$$_{Homo}$',
        'd_desoc_H_hetero': '$\Delta_{Dissoc}$$\it{H}$$_{Hetero}$',
        'E_form Metal': '$\it{E}$$_{form}$ Metal',
        'BG Metal': '$\it{E}$$_{gap}$ Metal',
        'Eads_CO2': '$\it{E}$$_{ads}$ C$O_{2}$',
        'Eads_H': '$\it{E}$$_{ads}$ H',
        'Eads_H2': '$\it{E}$$_{ads}$ $H_{2}$',
        'Eads_C': '$\it{E}$$_{ads}$ C',
        'Eads_N': '$\it{E}$$_{ads}$ N',
        'Eads_N2': '$\it{E}$$_{ads}$ $N_{2}$',
        'Eads_O': '$\it{E}$$_{ads}$ O',
        'Eads_O2': '$\it{E}$$_{ads}$ $O_{2}$',
        'E_form Oxide': '$\it{E}$$_{form}$ Oxide',
        'BG Oxide': '$\it{E}$$_{gap}$ Oxide',
        }
    
    supp_desc_dict = {
        'E_form': '$\it{E}$$_{form}$',
        'BG_PBEsol': '$\it{E}$$_{gap}$ PBEsol',
        'E_surf': '$\it{E}$$_{surf}$',
        'Evac': '$\it{E}$$_{vac}$',
        'BG MP': '$\it{E}$$_{gap}$ MP',
        'E_form MP': '$\it{E}$$_{form}$ MP',
        'd_form_H_298K':'$\Delta_{f}$$\it{H}$$^{\ominus}_{298K}$',
        'std_mol_S': 'S$^{\ominus}_{298K}$',
        'std_mol_Cp': 'C$^{\ominus}_{p}$',
        'Eads_H2': '$\it{E}$$_{ads}$ $H_{2}$',
        'Eads_CO2': '$\it{E}$$_{ads}$ C$O_{2}$',
        'Eads_CO': '$\it{E}$$_{ads}$ CO',
        'Eads_NO': '$\it{E}$$_{ads}$ NO',
        'Eads_N2O': '$\it{E}$$_{ads}$ $N_{2}$O',
        'Eads_CH4': '$\it{E}$$_{ads}$ C$H_{4}$',
        'Eads_CH3OH': '$\it{E}$$_{ads}$ C$H_{3}$OH',
        'Eads_C3H5OH': '$\it{E}$$_{ads}$ $C_{2}H_{5}$OH',
        'Edisads_HCOOH': '$\it{E}$$_{disads}$ HCOOH',
        'Edisads_CH3COOH': '$\it{E}$$_{disads}$ $CH_{3}$COOH',
        'N2 uptake at P/P0=0.015': '$N_{2}$ uptake at P/$P_{0}$=0.015',
        'N2 uptake at P/P0=0.015 at 700℃': '$N_{2}$ uptake at P/$P_{0}$=0.015 at 700℃',
        'IR L acid cm-1': 'IR L acid $cm^{-1}$',
        'IR B acid cm-1': 'IR B acid $cm^{-1}$',
        'IR L acid cm-1 cal at 700℃': 'IR L acid $cm^{-1}$ cal at 700℃',
        'IR B acid cm-1 cal at 700℃': 'IR B acid $cm^{-1}$ cal at 700℃',
        }

    supp_name_dict = {
        'MgO': 'MgO',
        'Al2O3_gamma': '$\\gamma$-$Al_{2}O_{3}$',
        'Al2O3_800Ccal': '$Al_{2}O_{3}$_800℃cal',
        'Al2O3_1000Ccal': '$Al_{2}O_{3}$_1000℃cal',
        'Al2O3_1100Ccal': '$Al_{2}O_{3}$_1100℃cal',
        'Al2O3_AKP-G07': '$Al_{2}O_{3}$_AKP-G07',
        'Al2O3_AKP-50': '$Al_{2}O_{3}$_AKP-50',
        'TiO2_STR-100N': '$TiO_{2}$_STR100N',
        'TiO2_ST01': '$TiO_{2}$_ST01',
        'TiO2_P25': '$TiO_{2}$_P25',
        'ZrO2_JRC5_NND': '$ZrO_{2}$_JRC5_NND',
        'ZrO2_monocliric': '$ZrO_{2}$_monocliric',
        'ZrO2_JRC4_EP': '$ZrO_{2}$_JRC4_EP',
        'ZrO2_JRC3_RC-100': '$ZrO_{2}$_JRC3_RC-100',
        'ZrO2_JRC2_EP-D': '$ZrO_{2}$_JRC2_EP-D',
        'Nb2O5': '$Nb_{2}O_{5}$',
        'CeO2_JRC3': '$CeO_{2}$_JRC3',
        'CeO2_JRC2': '$CeO_{2}$_JRC2',
        'CeO2(25%)-ZrO2': '$CeO_{2}$(25%)-$ZrO_{2}$',
        'CeO2(50%)-ZrO2': '$CeO_{2}$(50%)-$ZrO_{2}$',
        'Carbon': 'Carbon',
        'Mo2C': '$Mo_{2}$C',
        'WC': 'WC',
        'TaN': 'TaN',
        'ZrC': 'ZrC',
        'TaC': 'TaC',
        'SiO2': '$SiO_{2}$',
        'SiO2-Al2O3_JRC-SAL-6': '$SiO_{2}$-$Al_{2}O_{3}$_JRC-SAL-6',
        'SiO2-Al2O3_JRC-SAH-7': '$SiO_{2}$-$Al_{2}O_{3}$_JRC-SAH-7',
        'Co3O4': '$Co_{3}O_{4}$',
        'In2O3': '$In_{2}O_{3}$',
        'SnO2': 'Sn$O_{2}$',
        'BaTiO3': 'BaTi$O_{3}$',
        'SrTiO3': 'SrTi$O_{3}$',
        'CaTiO3': 'CaTi$O_{3}$',
        'SrFe12O19': 'Sr$Fe_{12}O_{19}$',
        'H-MOR_9': 'H-MOR_9',
        'H-MOR_10': 'H-MOR_10',
        'H-MOR_45': 'H-MOR_45',
        'H-MOR_110': 'H-MOR_110',
        'H-ZSM-5_11': 'H-ZSM-5_11',
        'H-ZSM-5_45': 'H-ZSM-5_45',
        'Na-ZSM-5_45': 'Na-ZSM-5_45',
        'H-Y_2.8': 'H-Y_2.8',
        'Na-Y_2.8': 'Na-Y_2.8',
        'H-Y_15': 'H-Y_15',
        'H-Y_50': 'H-Y_50',
        'H-BEA_9': 'H-BEA_9',
        'H-BEA_20': 'H-BEA_20',
        'H-BEA_75': 'H-BEA_75',
        'H-BEA_255': 'H-BEA_255',
        'H-CHA_11': 'H-CHA_11'
        }

    basic_desc_columns = elem_desc_columns[elem_desc_columns.index('AN'): elem_desc_columns.index('Solubility')+1]
    drop_desc_list = ['Symbol', 'Name', 'Homonuclear bonds', 'Heteronuclear bonds', 'Composition_Ox', 'Precursor']
    basic_desc_columns = [s for s in basic_desc_columns if s not in drop_desc_list]

    all_supp_desc = supp_desc_columns[supp_desc_columns.index('MW'): supp_desc_columns.index('Pore size')+1]


    noble_gas = ['He', 'Ne', 'Ar', 'Kr', 'Xe']
    pgm_plus_ReAu = ['H', 'Ru', 'Rh', 'Pd', 'Re', 'Os', 'Ir', 'Pt', 'Au']
    heavy_elems = ['Po', 'At', 'Rn', 'Fr',
                    'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu']
    alkali_metal = ['Li', 'Na', 'K', 'Rb', 'Cs']
    alkaline_earth_metal = ['Mg', 'Ca', 'Sr', 'Ba']  # drop 'Be'
    prm_elems = alkali_metal + alkaline_earth_metal

    if Reaction == 'rwgs_250' or Reaction == 'rwgs_300' or Reaction == 'rwgs_250_1wt' or Reaction == 'rwgs_TL':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'Cl', 'As',
            'Br', 'Tc', 'I', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au',
            'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu'
        ]

    elif Reaction == 'CH3OH' or Reaction == 'CH3OH_flow':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au',
            'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu'
        ]

    elif Reaction == 'N2O':
        # 'Re' is treated as PGMs.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Re', 'Os', 'Ir', 'Pt', 'Au',
            'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu'
        ]

    elif Reaction == 'N2O-SCR':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au',
            'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu'
            ]

    elif Reaction == 'H2SCR':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au',
            'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu'
        ]

    elif Reaction == 'NH3SCR':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au',
            'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu'
        ]

    elif Reaction == 'CH4':
        # 'Re' is treated as Additives.
        drop_elems = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'I', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au',
            'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu'
        ]
        # 2022/12/26 remove 'P', 'S', 'Se'

    elif Reaction == 'EtOH_CO2' or Reaction == 'EtOH_CO':
        # 'Re' is treated as Additives.
        drop_elems_1st = [
            'Be', 'B', 'C', 'N', 'O', 'F', 'Cl', 'As',
            'Br', 'Tc', 'Pm', 'Hg', 'Tl',
            'Ru', 'Rh', 'Pd', 'Os', 'Ir', 'Pt', 'Au',
            'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu'
        ]
        drop_elems = drop_elems_1st + prm_elems
    
    if Reaction == 'rwgs_250' or Reaction == 'rwgs_300' or Reaction == 'rwgs_250_1wt':
        target_name = 'CO formation rate_mmol min-1 gcat-1'
    elif Reaction == 'CH3OH':
        target_name = 'CH3OH rate (mmol g-1 h-1)'
    elif Reaction == 'EtOH_CO2' or Reaction == 'EtOH_CO':
        target_name = 'C2H5OH rate x selec'
    elif Reaction == 'H2SCR':
        target_name = 'Avg N2 yield_50-150'
    else:
        target_name = 'target'
    
    if Reaction == 'rwgs_250':
        eda_use_cols = [
            'No.', 'Catalyst',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
            'Ad4', 'Ad4_wt%', 'Ad5', 'Ad5_wt%',
            target_name, 'Iteration'
            ]
    
    elif Reaction == 'rwgs_250_1wt':
        eda_use_cols = [
            'No.', 'Catalyst',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4',
            'Ad4_wt%', 'Ad5', 'Ad5_wt%',
            'Support_name', 'Calc T',
            target_name, 'Iteration'
        ]
    
    elif Reaction == 'CH3OH':
        eda_use_cols = [
            'No.', 'Catalyst',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
            'Support_name', 'Calc T',
            target_name, 'Iteration'
        ]
    
    elif Reaction == 'N2O':
        eda_use_cols = [
            'No.', 'Catalyst',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%', 'PGM4', 'PGM4_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
            'Support_name', 'Calc T',
            target_name, 'Iteration'
        ]
    
    elif Reaction == 'N2O-SCR':
        eda_use_cols = [
            'Exp No.', 'Cat No.', 'Catalyst',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%', 'PGM4', 'PGM4_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
            'Support_name', 'Calc T', 'H2_ppm', 'CO_ppm', 'NH3_ppm', 'C3H6_ppm', 'CH4_ppm',
            target_name, 'Iteration'
        ]
    
    elif Reaction == 'CH4':
        eda_use_cols = [
            'No.', 'Catalyst',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
            'PGM4', 'PGM4_wt%', 'PGM5', 'PGM5_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
            'Support_name', 'Calc T',
            target_name, 'Iteration'
        ]
    
    elif Reaction == 'H2SCR':
        eda_use_cols = [
            'No.', 'Catalyst',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%', 'PGM3', 'PGM3_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
            'Support_name', 'Calc T',
            target_name, 'Iteration'
        ]
    
    elif Reaction == 'NH3SCR':
        eda_use_cols = [
            'No.', 'Catalyst',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%', 'Ad4', 'Ad4_wt%',
            'Support_name', 'Calc T',
            target_name, 'Iteration'
            ]
    
    elif Reaction == 'EtOH_CO2' or Reaction == 'EtOH_CO':
        eda_use_cols = [
            'Exp No.', 'Cat No.', 'Catalyst',
            'PGM1', 'PGM1_wt%', 'PGM2', 'PGM2_wt%',
            'Promoter1', 'Promoter1_wt%', 'Promoter2', 'Promoter2_wt%',
            'Ad1', 'Ad1_wt%', 'Ad2', 'Ad2_wt%', 'Ad3', 'Ad3_wt%',
            'Ad4', 'Ad4_wt%', 'Ad5', 'Ad5_wt%',
            'Support_name', 'Calc T', 'Reaction T', 'WHSV (mL g-1 h-1)',
            target_name, 'Iteration'
        ]

    desc_cols = {
        'noble_gas': noble_gas, 'heavy_elems': heavy_elems, 'pgm_plus_ReAu': pgm_plus_ReAu,
        'prm_elems': prm_elems, 'drop_elems': drop_elems,
        'elem_desc_columns': elem_desc_columns, 'elem_desc_dict': elem_desc_dict,
        'supp_desc_dict': supp_desc_dict, 'supp_name_dict': supp_name_dict,
        'basic_desc_columns': basic_desc_columns, 'supp_desc_columns': supp_desc_columns,
        'all_supp_desc': all_supp_desc, 'eda_use_cols': eda_use_cols
    }
    return desc_cols
