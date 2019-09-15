#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to run survival regression on features extracted from TCGA data.
"""

import pandas as pd
import numpy as np
from util_functions import safe_mkdir

import matplotlib.pyplot as plt
from sklearn.model_selection import *
from lifelines.utils import concordance_index
from statsmodels.stats.multitest import multipletests
from sklearn import preprocessing

from lifelines import CoxPHFitter
from sksurv.svm import FastSurvivalSVM

pd.set_option('display.max_columns', 200)
part_number = 2

plt.style.use('seaborn-white')

cancer_locs = ['BRCA', 'HNSC', 'BLCA']

colors_protected = {'baseline': '#A1E4C5', 'SciClone': '#AE88FF',
                    'MATH': '#C9AB3C', 'PhyloWGS': '#FF7DA4', 'CSR': '#9B9B9B',
                    'PyClone': '#84B8FF', 'expands': '#FFC89B',
                    'Combination': '#FFFFFF'}

tmp_df = pd.read_csv('tmp/20190903_ith_method_metrics_final_runtime_immunity.csv', sep='\t')

tmp_df = tmp_df[tmp_df.survival_days>0]
base_cols = ['nb_clones', 'clonal_prop', 'smallest_vaf', 'shannon_index', 'most_populated_clone_vaf']


# new approach, get only features that increase with clonality. So a few easy transformations
norm_base_cols = ['norm_{}'.format(c) for c in base_cols]

tmp_df = tmp_df.assign(bool_vital_status=tmp_df.binary_vital_status.astype(bool))
nb_clones_col = [i for i in tmp_df.columns if 'nb_clones' in i]


def identity(x):
    return x


def one_minus_x(x):
    return 1-x


def hundred_minus_x(x):
    return 100-x

function_dict = {'nb_clones': identity, 'clonal_prop': hundred_minus_x,
                 'smallest_vaf': one_minus_x, 'shannon_index': identity,
                 'most_populated_clone_vaf': one_minus_x}
for c in base_cols:
    concerned_cols = [col for col in tmp_df.columns if c in col]
    print(concerned_cols)
    for real_c in concerned_cols:
        tmp_df = tmp_df.assign(**{'norm_{}'.format(real_c): lambda x: function_dict[c](x[real_c])})
        tmp_df = tmp_df.assign(**{'norm_{}_squared'.format(real_c): lambda x: x['norm_{}'.format(real_c)]**2})
tmp_df = tmp_df.assign(**{'math_public_squared': lambda x: x['math_cn_public']**2})
tmp_df = tmp_df.assign(**{'math_protected_squared': lambda x: x['math_cn_protected']**2})
tmp_df = tmp_df.assign(**{'math_protected_abs_squared': lambda x: x['math_cn_protected_abs']**2})

cols_restricted = [c for c in tmp_df.columns if ('norm' in c)]

extra_cols = ['age_at_diagnosis',
              'purity',
              'purity_absolute',
              'ploidy',
              'ploidy_absolute',
              'perc_non_diploid_ascat',
              'perc_non_diploid_absolute',
              'mutation_count_protected',
              'mutation_count_public',
              'math_protected',
              'math_cn_protected',
              'math_public',
              'math_cn_public',
              'math_cn_protected_abs']

nb_clones_restricted = [c for c in cols_restricted if 'nb_clones' in c]
cols_restricted += extra_cols
nb_clones_restricted += extra_cols

tmp_df_c = tmp_df.dropna(subset=nb_clones_col).copy()
output_path = 'results/20190901_results_rerun_newsksurv'
safe_mkdir(output_path)
restricted_df_dict = dict()
restricted_df_dict['no_purity1'] = tmp_df_c

if __name__ == '__main__':
    clinical_data_dict = dict()
    for cancer_loc in cancer_locs:
        # load clinical data
        clinical_data = pd.read_csv('data/{}/clinical/data_bcr_clinical_data_patient.txt'.format(cancer_loc), sep='\t', skiprows=4)
        clinical_data_sample = pd.read_csv('data/{}/clinical/data_bcr_clinical_data_sample.txt'.format(cancer_loc), sep='\t', skiprows=4)
        clinical_data = clinical_data.assign(binary_vital_status=clinical_data.apply(lambda x: 1 if x.OS_STATUS == 'DECEASED' else (0 if x.OS_STATUS == 'LIVING' else np.nan), axis=1))
        clinical_data = clinical_data[clinical_data.OS_MONTHS != '[Not Available]']
        clinical_data = clinical_data[clinical_data.DAYS_TO_BIRTH != '[Not Available]']
        clinical_data = clinical_data.assign(survival_days=clinical_data.apply(lambda x: 30.5*float(x['OS_MONTHS']), axis=1))
        clinical_data = clinical_data.assign(age_at_diagnosis=clinical_data.apply(lambda x: int(np.round(-float(x.DAYS_TO_BIRTH)/365)), axis=1))
        clinical_data = clinical_data.assign(cancer_loc=cancer_loc)
        # Th objective is to select relevant clinical variables
        # eliminate variables that are constant for all patients
        clinical_data = clinical_data.loc[:,clinical_data.apply(pd.Series.nunique) != 1]

        # eliminate features with na values if the feature is numerical
        # if the feature is categorical, one hot encoding allows to encode
        # missing values easily, so we keep them.
        clinical_data_bkp = clinical_data.copy()
        clinical_data = clinical_data.replace('[Not Available]', np.nan)
        clinical_data = clinical_data.replace('[Not Applicable]', np.nan)

        clinical_data = clinical_data.apply(pd.to_numeric, errors='ignore')
        uu = clinical_data.dtypes
        object_cols = uu[uu==object].index.tolist()
        clinical_data.loc[:,object_cols] = clinical_data_bkp.loc[:, object_cols]

        clinical_data = clinical_data.dropna(axis=1)

        clinical_data = clinical_data.assign(AJCC_TUMOR_PATHOLOGIC_PT_simple = clinical_data.AJCC_TUMOR_PATHOLOGIC_PT.str[:2])
        clinical_data = clinical_data.assign(AJCC_NODES_PATHOLOGIC_PN_simple = clinical_data.AJCC_NODES_PATHOLOGIC_PN.str[:2])
        clinical_data = clinical_data.assign(AJCC_METASTASIS_PATHOLOGIC_PM_simple = clinical_data.AJCC_METASTASIS_PATHOLOGIC_PM.str[:2])
        cols_to_exclude = ['OTHER_PATIENT_ID', 'PATIENT_ID', 'FORM_COMPLETION_DATE', 'DAYS_TO_BIRTH', 'OS_STATUS',
                           'OS_MONTHS', 'AGE', 'AJCC_TUMOR_PATHOLOGIC_PT', 'AJCC_NODES_PATHOLOGIC_PN',
                           'AJCC_METASTASIS_PATHOLOGIC_PM', 'TUMOR_STATUS', 'VITAL_STATUS', 'DFS_STATUS',
                           'NEW_TUMOR_EVENT_AFTER_INITIAL_TREATMENT', 'TREATMENT_OUTCOME_FIRST_COURSE']
        clinical_data = clinical_data.drop(columns=set(cols_to_exclude).intersection(set(clinical_data.columns)))

        data_x_numeric = pd.get_dummies(clinical_data)
        uu = data_x_numeric.astype(object).describe(include='all')
        # avoid almost constant values
        data_x_numeric_filter = data_x_numeric.loc[:,(data_x_numeric.shape[0]-50>uu.loc['freq'])]

        all_cols = data_x_numeric_filter.columns.tolist()
        all_cols.remove('survival_days')
        all_cols.remove('binary_vital_status')
        summary_list = list()
        for i, cc in enumerate(all_cols):
            cph = CoxPHFitter()
            cph.fit(data_x_numeric_filter[['survival_days', 'binary_vital_status', cc]],
                    duration_col='survival_days', event_col='binary_vital_status',
                    show_progress=False)
            summary_list.append(cph.summary)

        big_summary = pd.concat(summary_list, axis=0, sort=True)

        big_summary = big_summary.assign(corrected_pvalue=multipletests(big_summary['p'].values.ravel(), method='fdr_bh')[1])

        final_selected_features = data_x_numeric[big_summary[big_summary.corrected_pvalue <0.05].index]
        big_summary[['exp(coef)', 'p', 'corrected_pvalue']].to_csv('tmp/{}_clinical_variable_significance.csv'.format(cancer_loc), sep='\t', float_format='%g')
        clinical_data_dict[cancer_loc] = pd.concat([clinical_data_bkp[['PATIENT_ID']], final_selected_features], axis=1, sort=True)



    class SurvivalStratifiedKFold(StratifiedKFold):
        def _make_test_folds(self, X, y=None):
            try:
                y = y[y.dtype.names[0]]
            except TypeError:
                y = y
            return super(SurvivalStratifiedKFold, self)._make_test_folds(X, y)

    def my_ci_pvalue(event_times, predicted_event_times, event_observed, N=10000):
        myCI = concordance_index(event_times, predicted_event_times, event_observed)
        n = 0
        for i in range(N):
            random_pred = np.random.permutation(predicted_event_times)
            if concordance_index(event_times, random_pred, event_observed) > myCI:
                n += 1
        try:
            pval = n/N
        except ZeroDivisionError:
            pval = 1
        return myCI, pval

    def apply_survival_regression(col_list, sub_restricted_tmp_df, regression_df, folder, ith_method, additional_cols, loc, tcga_y, test_res_df, extra_pred_train=None, extra_pred_test=None, repeats=0):
        X = pd.get_dummies(sub_restricted_tmp_df[col_list])

        X_train, X_test = X.loc[train_index, X.columns[2:]], X.loc[test_index, X.columns[2:]]
        y_train, y_test = tcga_y[train_index], tcga_y[test_index]

        X_train = X_train.loc[:,X_train.nunique()!=1]
        X_test = X_test.loc[:, X_train.columns]
        scaler = preprocessing.StandardScaler().fit(X_train)
        X_train_s = scaler.transform(X_train)
        X_test_s = scaler.transform(X_test)  

        XX_train = X.loc[train_index]
        XX_test = X.loc[test_index]
        XX_train = XX_train.loc[:,XX_train.nunique()!=1]
        XX_test = XX_test.loc[:, XX_train.columns]

        # survival svm
        lin_svm = FastSurvivalSVM(rank_ratio=0.8, fit_intercept=True, max_iter=200)
        lin_svm.fit(X_train_s, y_train)
        T_pred_train = lin_svm.predict(X_train_s)
        ci_train, pval_train = my_ci_pvalue(XX_train['survival_days'], T_pred_train, XX_train['binary_vital_status'], repeats)
        T_pred_test = lin_svm.predict(X_test_s)
        ci_test, pval_test = my_ci_pvalue(XX_test['survival_days'], T_pred_test, XX_test['binary_vital_status'], repeats)
        variable_base_name = [r.split('_{}_'.format(folder))[1] if ith_method in r else r for r in X_train.columns]
        coef_dict = dict(zip(['{}_coef'.format(c) for c in variable_base_name], lin_svm.coef_))
        t = pd.DataFrame({'cancer_loc':loc, 'ith_method': ith_method, 'folder': folder, 'additional_cols': additional_cols[0], 'regression_method': 'linear_survival_svm', 'train_score': ci_train, 'pval_train': pval_train, 'test_score': ci_test, 'pval_test': pval_test, **coef_dict}, index=[0])
        regression_df = pd.concat((regression_df, t), sort=True)
        t = pd.DataFrame({'ith_method': ith_method, 'folder': folder, 'additional_cols': additional_cols[0], 'train_score': ci_train, 'pval_train': pval_train, 'test_score': ci_test, 'pval_test': pval_test, **{obs_idx: T_pred_test[obs_idx] for obs_idx in range(len(T_pred_test))}}, index=[0])
        test_res_df = pd.concat((test_res_df, t), sort=True)
        if extra_pred_test is not None:
            for i, c in enumerate(extra_pred_test):
                ci_test, pval_test = my_ci_pvalue(XX_test['survival_days'], T_pred_test+extra_pred_test[i], XX_test['binary_vital_status'], 0)
                ci_train, pval_train = my_ci_pvalue(XX_train['survival_days'], T_pred_train+extra_pred_train[i], XX_train['binary_vital_status'], 0)
                t = pd.DataFrame({'cancer_loc':loc, 'ith_method': ith_method, 'folder': folder, 'additional_cols': additional_cols[i+1], 'regression_method': 'linear_survival_svm', 'train_score': ci_train, 'pval_train': pval_train, 'test_score': ci_test, 'pval_test': pval_test}, index=[0])
                regression_df = pd.concat((regression_df, t), sort=True)
                t = pd.DataFrame({'ith_method': ith_method, 'folder': folder, 'additional_cols': additional_cols[i+1], 'train_score': ci_train, 'pval_train': pval_train, 'test_score': ci_test, 'pval_test': pval_test, **{obs_idx: T_pred_test+extra_pred_test[i][obs_idx] for obs_idx in range(len(T_pred_test+extra_pred_test[i]))}}, index=[0])
                test_res_df = pd.concat((test_res_df, t), sort=True)
        return regression_df, test_res_df, T_pred_test, T_pred_train




    NGS_no_ith_cols = ['ploidy', 'purity', 'perc_non_diploid_ascat',
                       'mutation_count_protected', 'mutation_count_public']
    pval_repeats = 1
    for k in restricted_df_dict.keys():
        tmp_res_path = 'tmp/20190901_{}_test_results'.format(k)
        safe_mkdir(tmp_res_path)
        nb_repeats = 1
        np.random.seed(567)
        restricted_tmp_df = restricted_df_dict[k]
        regression_df = pd.DataFrame(columns=['cancer_loc', 'ith_method', 'folder', 'regression_method', 'additional_cols', 'train_score', 'test_score'] + ['{}_coef'.format(c) for c in base_cols] + ['{}_pvalue'.format(c) for c in base_cols])
        for loc in cancer_locs:
            clinical_cols = clinical_data_dict[loc].columns.tolist()
            sub_restricted_tmp_df = restricted_tmp_df[restricted_tmp_df.cancer_loc == loc]
            sub_restricted_tmp_df = pd.merge(sub_restricted_tmp_df, clinical_data_dict[loc], left_on='PATIENT_ID', right_on='PATIENT_ID', how='left', suffixes=['', '_y'])
            sub_restricted_tmp_df.dropna(subset=clinical_cols+NGS_no_ith_cols, inplace=True)
            sub_restricted_tmp_df.index = range(len(sub_restricted_tmp_df))
            tcga_y = sub_restricted_tmp_df[['bool_vital_status', 'survival_days']].to_records(index=False)
            tcga_y.dtype.names=['vital_status', 'survival_days']
            for i in range(nb_repeats):
                cv = SurvivalStratifiedKFold(n_splits=5, shuffle=True)
                nb_cv = 0
                for train_index, test_index in cv.split(np.zeros(len(tcga_y)), tcga_y):
                    nb_cv+=1
                    y_train, y_test = tcga_y[train_index], tcga_y[test_index]
                    test_res_df = pd.DataFrame(columns=['ith_method', 'folder', 'additional_cols', 'train_score', 'test_score', 'pval_train', 'pval_test'] + [c for c in range(len(test_index))])
                    test_res_df = pd.concat((test_res_df, pd.DataFrame({'ith_method': 'groundTruth', 'folder': '', 'additional_cols': 'survival_days', **{obs_idx: y_test['survival_days'][obs_idx] for obs_idx in range(len(y_test))}}, index=[0])), sort=True)
                    test_res_df = pd.concat((test_res_df, pd.DataFrame({'ith_method': 'groundTruth', 'folder': '', 'additional_cols': 'vital_status', **{obs_idx: y_test['vital_status'][obs_idx] for obs_idx in range(len(y_test))}}, index=[0])), sort=True)
                    total_col_list = list()
                    total_col_list_squared = list()
                    for ith_method, folder in [('pyclone', 'protected_hg38_vcf'), ('pyclone', 'public_hg38_vcf'),
                                               ('pyclone', 'protected_hg38_vcf_absolute'),
                                               ('sciclone', 'protected_hg38_vcf'), ('sciclone', 'public_hg38_vcf'),
                                               ('sciclone', 'protected_hg38_vcf_absolute'),
                                               ('PhyloWGS', 'protected_hg38_vcf'), ('PhyloWGS', 'public_hg38_vcf'),
                                               ('PhyloWGS', 'protected_hg38_vcf_absolute'),
                                               ('CSR', 'protected_hg38_vcf'), ('CSR', 'public_hg38_vcf'),
                                               ('CSR', 'protected_hg38_vcf_absolute'),
                                               ('expands', 'protected_hg38_vcf'), ('expands', 'public_hg38_vcf'),
                                               ('expands', 'protected_hg38_vcf_absolute'),]:
                        col_list = ['survival_days', 'binary_vital_status'] + ['norm_{}_{}_{}'.format(ith_method, folder, c) for c in base_cols]
                        squared_cols = ['norm_{}_{}_{}_squared'.format(ith_method, folder, c) for c in base_cols]
                        if sum(pd.isnull(sub_restricted_tmp_df[col_list[2]])):
                            continue
                        total_col_list += col_list[2:]
                        total_col_list_squared += squared_cols
                        regression_df, test_res_df, _, _ = apply_survival_regression(col_list, sub_restricted_tmp_df, regression_df, folder, ith_method, ['clonality'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                        regression_df, test_res_df, _, _ = apply_survival_regression(col_list[0:3], sub_restricted_tmp_df, regression_df, folder, ith_method, ['nb_clones'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                        regression_df, test_res_df, _, _ = apply_survival_regression(col_list+squared_cols, sub_restricted_tmp_df, regression_df, folder, ith_method, ['clonality_squared'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                        regression_df, test_res_df, _, _ = apply_survival_regression(col_list[0:3]+squared_cols[0:1], sub_restricted_tmp_df, regression_df, folder, ith_method, ['nb_clones_squared'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                        regression_df, test_res_df, _, _ = apply_survival_regression(col_list+clinical_cols, sub_restricted_tmp_df, regression_df, folder, ith_method, ['clonality+clinical'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                        regression_df, test_res_df, _, _ = apply_survival_regression(col_list+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, folder, ith_method, ['clonality+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                        regression_df, test_res_df, _, _ = apply_survival_regression(col_list+clinical_cols+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, folder, ith_method, ['clonality+clinical+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_public'], sub_restricted_tmp_df, regression_df, 'public_hg38_vcf', 'math_score', ['clonality'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_public', 'math_public_squared'], sub_restricted_tmp_df, regression_df, 'public_hg38_vcf', 'math_score', ['clonality_squared'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_public']+clinical_cols, sub_restricted_tmp_df, regression_df, 'public_hg38_vcf', 'math_score', ['clonality+clinical'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_public']+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'public_hg38_vcf', 'math_score', ['clonality+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_public']+clinical_cols+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'public_hg38_vcf', 'math_score', ['clonality+clinical+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_protected'], sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf', 'math_score', ['clonality'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_protected', 'math_protected_squared'], sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf', 'math_score', ['clonality_squared'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_protected']+clinical_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf', 'math_score', ['clonality+clinical'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_protected']+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf', 'math_score', ['clonality+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_protected']+clinical_cols+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf', 'math_score', ['clonality+clinical+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_protected_abs'], sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf_absolute', 'math_score', ['clonality'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_protected_abs', 'math_protected_abs_squared'], sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf_absolute', 'math_score', ['clonality_squared'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_protected_abs']+clinical_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf_absolute', 'math_score', ['clonality+clinical'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_protected_abs']+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf_absolute', 'math_score', ['clonality+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status', 'math_cn_protected_abs']+clinical_cols+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf_absolute', 'math_score', ['clonality+clinical+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    public_cols = ['survival_days', 'binary_vital_status', 'math_cn_public'] + [c for c in total_col_list if 'public' in c]
                    public_cols_squared = public_cols + ['math_public_squared'] + [c for c in total_col_list_squared if 'public' in c]
                    protected_cols = ['survival_days', 'binary_vital_status', 'math_cn_protected'] + [c for c in total_col_list if (('protected' in c) and ('absolute' not in c))]
                    protected_cols_squared = protected_cols + ['math_protected_squared'] + [c for c in total_col_list_squared if (('protected' in c) and ('absolute' not in c))]
                    protected_abs_cols = ['survival_days', 'binary_vital_status', 'math_cn_protected_abs'] + [c for c in total_col_list if (('protected' in c) and ('absolute' in c))]
                    protected_abs_cols_squared = protected_abs_cols + ['math_protected_abs_squared'] + [c for c in total_col_list_squared if (('protected' in c) and ('absolute' in c))]
                    regression_df, test_res_df, _, _ = apply_survival_regression(public_cols, sub_restricted_tmp_df, regression_df, 'public_hg38_vcf', 'combination', ['clonality'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(public_cols_squared, sub_restricted_tmp_df, regression_df, 'public_hg38_vcf', 'combination', ['clonality_squared'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(public_cols+clinical_cols, sub_restricted_tmp_df, regression_df, 'public_hg38_vcf', 'combination', ['clonality+clinical'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(public_cols+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'public_hg38_vcf', 'combination', ['clonality+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(public_cols+clinical_cols+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'public_hg38_vcf', 'combination', ['clonality+clinical+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(protected_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf', 'combination', ['clonality'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(protected_cols_squared, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf', 'combination', ['clonality_squared'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(protected_cols+clinical_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf', 'combination', ['clonality+clinical'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(protected_cols+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf', 'combination', ['clonality+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(protected_cols+clinical_cols+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf', 'combination', ['clonality+clinical+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)

                    regression_df, test_res_df, _, _ = apply_survival_regression(protected_abs_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf_absolute', 'combination', ['clonality'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(protected_abs_cols_squared, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf_absolute', 'combination', ['clonality_squared'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(protected_abs_cols+clinical_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf_absolute', 'combination', ['clonality+clinical'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(protected_abs_cols+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf_absolute', 'combination', ['clonality+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(protected_abs_cols+clinical_cols+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'protected_hg38_vcf_absolute', 'combination', ['clonality+clinical+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status'] + clinical_cols, sub_restricted_tmp_df, regression_df, 'none', 'none', ['clonality+clinical'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status'] +NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'none', 'none', ['clonality+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    regression_df, test_res_df, _, _ = apply_survival_regression(['survival_days', 'binary_vital_status'] + clinical_cols+NGS_no_ith_cols, sub_restricted_tmp_df, regression_df, 'none', 'none', ['clonality+clinical+NGS'], loc, tcga_y, test_res_df, repeats=pval_repeats)
                    test_res_df.to_csv('{}/{}_fold{}_test_res.csv'.format(tmp_res_path, loc, nb_cv), index=False, sep='\t')

        regression_df = regression_df.assign(method_folder=regression_df[['ith_method', 'folder']].astype(str).sum(axis=1))
        regression_df.to_csv('tmp/20190901_regression_df_{}_scale_mean_model_squared.csv'.format(k), sep='\t', index=False)



