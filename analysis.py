#!/usr/bin/env python
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt

import sklearn
from sklearn.ensemble import ExtraTreesRegressor, RandomForestRegressor
from sklearn.model_selection import ShuffleSplit, KFold, GridSearchCV, train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from sklearn.model_selection import train_test_split, LeaveOneOut
from adjustText import adjust_text
import shap

import sys, os
sys.dont_write_bytecode = True
import conditions

import random
random.seed(1107)
np.random.seed(1107)

import warnings
warnings.filterwarnings('ignore')

PATH = 'Figures'
os.makedirs(PATH, exist_ok = True)

plt.rcParams['font.size'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['font.family'] = 'Hiragino sans'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['xtick.major.width'] = 1.0
plt.rcParams['ytick.major.width'] = 1.0
plt.rcParams['axes.linewidth'] = 1.0
plt.rcParams['axes.grid'] = False
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.linewidth'] = 2
plt.rcParams["legend.markerscale"] = 2
plt.rcParams["legend.fancybox"] = False
plt.rcParams["legend.framealpha"] = 1
plt.rcParams["legend.edgecolor"] = 'black'
plt.rcParams.update({'mathtext.default':  'regular' })

def analysis_data_convert(condition, data_sheet, use_models=[], idx=None):
    """_summary_

    Args:
        condition (_type_): _description_
        use_models (list, optional): _description_. Defaults to [2,2,2].
        idx (_type_, optional): _description_. Defaults to None.
    """
    converted = {}
    condition = conditions.calc_condition()
    date = condition['date']
    Reaction = condition['Reaction']

    data_file_name = f'data/{date}_{Reaction}_dataset.xlsx'
    desc_file_name = 'data/Descriptors.xlsx'
    elem_desc_sheet = 'Descriptors_elem'

    data_cols = conditions.data_columns(condition)
    elem_cols, wt_cols = data_cols['elem'], data_cols['wt']
    target_name = condition['target_name']

    if len(use_models) == 0:
        add_model = condition['add_model']
    else:
        add_model = use_models[0]
    add_num = condition['add_num']

    desc_cols = conditions.desc_columns()
    basic_desc_cols = desc_cols['basic_desc_columns']
    noble_gas = desc_cols['noble_gas']
    drop_elems = desc_cols['drop_elems']

    if Reaction == 'rwgs_250':
        add_use_desc = data_cols['add_desc']
    else:
        desc = conditions.use_desc()
        add_use_desc = desc['add']

	### Read dataset and descriptors ###
    data = pd.read_excel(data_file_name, sheet_name=data_sheet)
    data.loc[:, elem_cols] = data.loc[:, elem_cols].fillna('H')
    data.loc[:, wt_cols] = data.loc[:, wt_cols].fillna(0)

    if idx == None:
        print(date, Reaction, 'all data')
    else:
        print(Reaction, 'Iteration =', idx)
        data_idx = data[data['Iteration'] <= idx]
        data_idx = data_idx.reset_index()
        data = data_idx
    
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
        feat_cols = list(add_desc.index[1:])

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
        feat_cols = list(add_desc.index[1:]) + list(add_desc.columns)

    elif add_model == 2 and add_num != 0:
        feat = sum([
            np.multiply(np.array(add_desc.loc[data.loc[:, f"Ad{i+1}"]]),
            np.array(data[f"Ad{i+1}_wt%"]).reshape(-1, 1)) /
            np.array(add_aw.loc[data.loc[:, f"Ad{i+1}"]]).reshape(-1, 1) for i in range(add_num)
            ])
        feat_cols = list(add_desc.columns)
    feat = pd.DataFrame(feat, columns=feat_cols)

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

def grid_search(feat, target, split_type='KFold', split_num=10, use_model='ETR'):
    print(use_model)
    if split_type == 'KFold':
        cvf = KFold(n_splits=split_num, shuffle=True, random_state=1107)
    elif split_type == 'ShuffleSplit':
        cvf = ShuffleSplit(n_splits=split_num, random_state=1107, test_size=0.2)
    if 'ETR' == use_model:
        cvmodel = GridSearchCV(ExtraTreesRegressor(n_jobs=1, random_state=1107),
                               param_grid={"n_estimators": [100, 250, 500, 1000, 1500]},
                               n_jobs=4,cv=5)
        error = crossvalid(feat, target, cvmodel, cvf)

        model = ExtraTreesRegressor(n_estimators=cvmodel.best_params_['n_estimators'],
                        n_jobs=4, random_state=1107)

    if 'RFR' == use_model:
        cvmodel = GridSearchCV(RandomForestRegressor(n_jobs=1, random_state=1107),
                               param_grid={"n_estimators": [100, 250, 500, 1000, 1500]},
                               n_jobs=4,cv=5)
        error = crossvalid(feat, target, cvmodel, cvf)

        model = RandomForestRegressor(n_estimators=cvmodel.best_params_['n_estimators'],
                        n_jobs=4, random_state=1107)
    print('best_params:', cvmodel.best_params_['n_estimators'])
    return model, error

def crossvalid(xx, yy, model, cvf):
    mse_trn = []
    mse_tes = []
    mae_trn = []
    mae_tes = []
    r_2_tes = []
    r_2_trn = []
    for train_index, test_index in cvf.split(xx):
        x_trn = np.array(xx)[train_index]
        x_tes = np.array(xx)[test_index]
        y_trn = np.array(yy)[train_index]
        y_tes = np.array(yy)[test_index]
        model.fit(x_trn, y_trn)
        x_trn_pred = model.predict(x_trn)
        x_tes_pred = model.predict(x_tes)

        mse_tes.append(mean_squared_error(x_tes_pred, y_tes))
        mse_trn.append(mean_squared_error(x_trn_pred, y_trn))
        mae_tes.append(mean_absolute_error(x_tes_pred, y_tes))
        mae_trn.append(mean_absolute_error(x_trn_pred, y_trn))
        r_2_tes.append(r2_score(y_tes, x_tes_pred))
        r_2_trn.append(r2_score(y_trn, x_trn_pred))

    rmse_tes = np.sqrt(np.array(mse_tes))
    rmse_trn = np.sqrt(np.array(mse_trn))
    print("Train ... RMSE: %1.3f, MSE: %1.3f, MAE: %1.3f, R2: %1.3f, RMSE detail (sd: %1.3f, min:%1.3f, max:%1.3f)" 
    % (rmse_trn.mean(), np.array(mse_trn).mean(), np.array(mae_trn).mean(), np.array(r_2_trn).mean(), rmse_trn.std(), rmse_trn.min(), rmse_trn.max()))
    print("Test ... RMSE: %1.3f, MSE: %1.3f, MAE: %1.3f, R2: %1.3f, RMSE detail (sd: %1.3f, min:%1.3f, max:%1.3f)" 
    % (rmse_tes.mean(), np.array(mse_tes).mean(), np.array(mae_tes).mean(), np.array(r_2_tes).mean(), rmse_tes.std(), rmse_tes.min(), rmse_tes.max()))
    ret = {}
    ret['trn_rmse'] = rmse_trn.mean()
    ret['trn_mse'] = np.array(mse_trn).mean()
    ret['trn_mae'] = np.array(mae_trn).mean()
    ret['trn_r2'] = np.array(r_2_trn).mean()
    ret['tes_rmse'] = rmse_tes.mean()
    ret['tes_mse'] = np.array(mse_tes).mean()
    ret['tes_mae'] = np.array(mae_tes).mean()
    ret['tes_r2'] = np.array(r_2_tes).mean()
    return ret

def one_shot_plot(condition, feat, target, model, Reaction, test_size=0.1, random_state=1107, save=False):
    date, ML_model = condition['date'], condition['ML_model']
    add_model = condition['add_model']

    plt.figure(facecolor='white', figsize=(6,6))
    plt.subplot().set_aspect('equal')
    x_train, x_test, y_train, y_test = train_test_split(feat, target, test_size=test_size, random_state=random_state)
    model.fit(x_train, y_train)
    y_train_pred, y_test_pred = model.predict(x_train), model.predict(x_test)
    
    plt.plot(y_test, y_test_pred, 'o', c='red', markersize=5, alpha=0.6, label='Test')
    plt.plot(y_train, y_train_pred, 'o', c='blue', markersize=5, alpha=0.6, label='Train')
    
    plt.plot([0, int(target.max())+0.5], [0, int(target.max())+0.5], c='0', ls='-', lw=1.0)
    plt.legend()
    plt.xlabel('Experimental CO form rate (mmol $min^{{-}1}$ $gcat^{{-}1}$)')
    plt.ylabel('Predicted CO form rate (mmol $min^{{-}1}$ $gcat^{{-}1}$)')

    if save is not False:
        plt.savefig(f'{PATH}/Fig1_{date}_one-shot_{ML_model}_prop{add_model}.png', bbox_inches="tight", dpi=600)
    else:
        pass
    print(f'{int(1/test_size)}-fold One-shot plot {ML_model} prop{add_model} was finished.')

def loo_plot(condition, feat, target, feat_cols, model, plot_texts=True, save=False):
    loo = LeaveOneOut()
    date, ML_model, Reaction = condition['date'], condition['ML_model'], condition['Reaction']
    add_model = condition['add_model']
    
    for train_index, test_index in loo.split(feat):
        X_train, X_test = feat.loc[train_index, feat_cols], feat.loc[test_index, feat_cols]
        y_train, y_test = target.iloc[train_index], target.iloc[test_index]
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        feat.loc[test_index, 'original'] = y_test
        feat.loc[test_index, 'pred'] = y_pred
        feat.loc[test_index, 'RMSE'] = math.sqrt(mean_squared_error(y_test, y_pred))
        
    plt.figure(facecolor='white', figsize=(6, 6))
    f_max = feat.loc[:, 'original'].max()*1.1
    plt.plot([0, f_max], [0, f_max], c = 'black')
    plt.scatter(feat.loc[:, 'original'], feat.loc[:, 'pred'], color='red', s=25, alpha=0.5)
    plt.xlim([0, f_max])
    plt.ylim([0, f_max])

    plt.xlabel('Experimental CO form rate (mmol $min^{{-}1}$ $gcat^{{-}1}$)')
    plt.ylabel('Predicted CO form rate (mmol $min^{{-}1}$ $gcat^{{-}1}$)')

    if plot_texts == True:
        texts = [plt.text(feat.loc[i, 'original'], feat.loc[i, 'pred'], f"{i+1}", fontsize=7) for i in range(len(feat))]
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', linewidth=0.4, alpha=0.8))
    else:
        pass

    if save is not False:
        plt.savefig(f'{PATH}/Fig2_{date}_LOO_{ML_model}_prop{add_model}.png', bbox_inches="tight", dpi=1200)
    else:
        pass
    print(f'LeaveOneOut plot {ML_model} prop{add_model} was finished.')

def plot_importance(condition, model, labels, topk, save=False):
    date, ML_model = condition['date'], condition['ML_model']
    add_model = condition['add_model']

    plt.figure(facecolor='white', figsize =(6,6))
    importances = model.feature_importances_
    indices = np.argsort(importances)
    topk_idx = indices[-topk:]
    plt.barh(range(len(topk_idx)), importances[topk_idx], color='blue', edgecolor='blue', alpha=0.6, align='center')
    plt.yticks(range(len(topk_idx)), labels[topk_idx])
    plt.ylim([-1, len(topk_idx)])
    plt.xlabel("Feature Importance")

    if save is not False:
        plt.savefig(f'{PATH}/Fig3_{date}_Feature_importance_{ML_model}_prop{add_model}.png', bbox_inches="tight", dpi=600)
    else:
        pass
    print(f'Feature_importance plot {ML_model} prop{add_model} was finished.')

def shap_summary_plot(condition, model, feat, target, save=False):
    date, ML_model = condition['date'], condition['ML_model']
    add_model = condition['add_model']
    model.fit(feat, target)
    explainer = shap.TreeExplainer(model=model)
    shap_values = explainer.shap_values(feat)
    shap.summary_plot(
        shap_values,  features=feat, feature_names=feat.columns, max_display=15,
        plot_type='dot', color=None, axis_color='k', title=None, alpha=0.8,
        show=False, sort=True, color_bar=True, plot_size='auto', layered_violin_max_num_bins=20,
        class_names=None, class_inds=None, color_bar_label='Feature value',
        #cmap='cool',
        #cmap='cool_r',
        auto_size_plot=None, use_log_scale=False)
    if save is not False:
        plt.savefig(f'{PATH}/Fig4_{date}_SHAP_sumplot_{ML_model}_prop{add_model}.png', bbox_inches='tight', dpi=600)
    else:
        pass
    print(f'SHAP_summary_plot {ML_model} prop{add_model} was finished.')