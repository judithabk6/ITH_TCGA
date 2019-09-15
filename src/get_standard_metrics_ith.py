#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to take as input a patient, and to compute all metrics (1B, 1C, 2A)
between all the methods to get a complete overview
"""

import pandas as pd
import numpy as np
from util_functions import safe_mkdir
import sys
import scipy as sp
import itertools
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from sklearn.metrics.cluster import v_measure_score
from scipy.spatial.distance import cosine, euclidean


"""
# to execute just once to save the patient list to run jobs
tmp_df = pd.read_csv('tmp/20190903_ith_method_metrics_final_runtime_immunity.csv', sep='\t')
nb_clones_col = [i for i in tmp_df.columns if 'nb_clones' in i]
tmp_df_c = tmp_df.dropna(subset=nb_clones_col).copy()
patient_df = tmp_df_c[['PATIENT_ID']]
patient_df.index = list(range(1, len(patient_df)+1))
patient_df.to_csv('intersection_patient_list.csv', sep=',', header=False, index=True)
"""

patient_id = sys.argv[1]

# patient_id = 'TCGA-A7-A0CG'

def score1B_base(J_true, J_pred):
    return (J_true + 1 - min([J_true + 1, np.abs(J_true - J_pred)])) / \
        (J_true + 1)


def score1C_base(phi_true_values, phi_pred_values, weights_true=None,
                 weights_pred=None):
    return 1 - sp.stats.wasserstein_distance(phi_true_values, phi_pred_values,
                                             weights_true, weights_pred)


def score2A_base(true_cluster_assign, pred_cluster_assign):
    # co-clustering matrix - hard assignement
    # step1 build the co-clustering matrix
    N = len(true_cluster_assign)
    sim_mat = np.zeros((N + 1, N + 1))

    for u in np.unique(true_cluster_assign):
        sim_mat[tuple(zip(*list(
            itertools.product(np.where(true_cluster_assign == u)[0],
                              np.where(true_cluster_assign == u)[0]))))] = 1
    sim_mat[N, N] = 1
    pred_mat = np.zeros((N + 1, N + 1))
    for u in np.unique(pred_cluster_assign):
        pred_mat[tuple(zip(*list(
            itertools.product(np.where(pred_cluster_assign == u)[0],
                              np.where(pred_cluster_assign == u)[0]))))] = 1
    pred_mat[N, N] = 1
    p = sp.stats.pearsonr(sim_mat.flatten(), pred_mat.flatten())[0]
    m = matthews_corrcoef(sim_mat.flatten(), pred_mat.flatten())
    v = v_measure_score(true_cluster_assign, pred_cluster_assign)

    good_scen_matrix = sim_mat.copy()
    bad1 = np.identity(N + 1)
    bad2 = np.ones((N + 1, N + 1))
    bad2[:-1, -1] = 0
    bad2[-1, :-1] = 0

    p_good = sp.stats.pearsonr(sim_mat.flatten(),
                               good_scen_matrix.flatten())[0]
    p_bad1 = sp.stats.pearsonr(sim_mat.flatten(), bad1.flatten())[0]
    p_bad2 = sp.stats.pearsonr(sim_mat.flatten(), bad2.flatten())[0]
    p_bad = min(p_bad1, p_bad2)
    pn = max(0, - p / (p_bad - p_good) + p_bad / (p_bad - p_good))

    m_good = matthews_corrcoef(sim_mat.flatten(), good_scen_matrix.flatten())
    m_bad1 = matthews_corrcoef(sim_mat.flatten(), bad1.flatten())
    m_bad2 = matthews_corrcoef(sim_mat.flatten(), bad2.flatten())
    m_bad = min(m_bad1, m_bad2)
    mn = max(0, - m / (m_bad - m_good) + m_bad / (m_bad - m_good))

    v_good = v_measure_score(sim_mat.flatten(), good_scen_matrix.flatten())
    v_bad1 = v_measure_score(sim_mat.flatten(), bad1.flatten())
    v_bad2 = v_measure_score(sim_mat.flatten(), bad2.flatten())
    v_bad = min(v_bad1, v_bad2)
    vn = max(0, - v / (v_bad - v_good) + v_bad / (v_bad - v_good))
    return np.mean([pn, mn, vn])

# get all arguments (number of clones, phi_values, weights, cluster_assignment) for each method
# merge to the input.tsv to be sure of the order of mutations
result_dict = dict()
for folder in ('protected_hg38_vcf', 'public_hg38_vcf', 'protected_hg38_vcf_absolute'):
    input_df = pd.read_csv('results/{}/pyclone/{}/input.tsv'.format(patient_id, folder), sep='\t')
    input_df = input_df.assign(sub_id=input_df.mutation_id.str[:-3])
    input_df = input_df.drop_duplicates(subset=['sub_id'], keep='first')
    # get pyclone values
    pyclone_results = pd.read_csv('results/{}/pyclone/{}/tables/loci.tsv'.format(patient_id, folder), sep='\t')
    pyclone_results = pyclone_results.assign(sub_id=pyclone_results.mutation_id.str[:-3])
    pyclone_results = pyclone_results.drop_duplicates(subset=['sub_id'], keep='first')
    inp_pyclone = pd.merge(input_df, pyclone_results, on='sub_id')

    nb_clones = inp_pyclone.cluster_id.nunique()
    phi_values = pyclone_results.groupby('cluster_id').cellular_prevalence.mean().values
    weights = pyclone_results.groupby('cluster_id').cellular_prevalence.count().values
    cluster_assignment = inp_pyclone[['sub_id', 'cluster_id']]
    result_dict['pyclone_{}'.format(folder)] = [nb_clones, phi_values, weights, cluster_assignment]

    # get sciclone values
    sciclone_results = pd.read_csv('results/{}/sciclone/{}/clusters1'.format(patient_id, folder), sep='\t')
    sciclone_results = sciclone_results.assign(full_chrom='chr' + sciclone_results.chr.astype(str))
    sciclone_results = sciclone_results.assign(sub_id=sciclone_results.full_chrom + '_' + sciclone_results.st.astype(str).str.replace('\\.0', ''))
    sciclone_results = sciclone_results.drop_duplicates(subset=['sub_id'], keep='first')
    sciclone_results = sciclone_results.assign(cluster_id=sciclone_results.cluster)
    inp_sciclone = pd.merge(input_df, sciclone_results, on='sub_id', how='left')
    nb_clones = inp_sciclone.cluster_id.nunique()
    phi_values = inp_sciclone.groupby('cluster_id')['tumor.vaf'].mean().values / 100
    weights = inp_sciclone.groupby('cluster_id')['tumor.vaf'].count().values
    cluster_assignment = inp_sciclone[['sub_id', 'cluster_id']]
    result_dict['sciclone_{}'.format(folder)] = [nb_clones, phi_values, weights, cluster_assignment]

    # get phylowgs results
    phylowgs_input = pd.read_csv('results/{}/PhyloWGS/{}/parsed_ssm_data.txt'.format(patient_id, folder), sep='\t')
    phylowgs_input = phylowgs_input.assign(sub_id='chr' + phylowgs_input.gene)
    phylowgs_results = pd.read_csv('results/{}/PhyloWGS/{}/outputs/2A.txt'.format(patient_id, folder), header=None, names=['cluster_id'])
    phylowgs_input = phylowgs_input.assign(cluster_id=phylowgs_results.cluster_id)
    inp_phylowgs = pd.merge(input_df, phylowgs_input, on='sub_id', how='left')
    nb_clones = inp_phylowgs.cluster_id.nunique()
    phylowgs_clust = pd.read_csv('results/{}/PhyloWGS/{}/outputs/1C.txt'.format(patient_id, folder), sep='\t', header=None, names=['cluster_idx', 'cluster_size', 'cluster_vaf'])
    phi_values = phylowgs_clust.cluster_vaf.values
    weights = phylowgs_clust.cluster_size.values
    cluster_assignment = inp_phylowgs[['sub_id', 'cluster_id']]
    result_dict['phylowgs_{}'.format(folder)] = [nb_clones, phi_values, weights, cluster_assignment]

    # get expands results
    expands_results = pd.read_csv('results/{}/expands/{}/{}__{}.sps'.format(patient_id, folder, patient_id, folder), sep='\t', skiprows=1)
    expands_results = expands_results.assign(full_chrom='chr' + expands_results.chr.astype(str))
    expands_results = expands_results.assign(sub_id=expands_results.full_chrom + '_' + expands_results.startpos.astype(str).str.replace('\\.0', ''))
    expands_results = expands_results.assign(cluster_id=expands_results.SP.apply(lambda x: expands_results.SP.unique().tolist().index(x)))
    expands_results = expands_results.drop_duplicates(subset=['sub_id'], keep='first')
    inp_expands = pd.merge(input_df, expands_results, on='sub_id', how='left')
    nb_clones = inp_expands.cluster_id.nunique()
    phi_values = inp_expands.groupby('cluster_id').SP.mean().values
    weights = inp_expands.groupby('cluster_id').SP.count().values
    cluster_assignment = inp_expands[['sub_id', 'cluster_id']]
    result_dict['expands_{}'.format(folder)] = [nb_clones, phi_values, weights, cluster_assignment]

    # get CSR results
    csr_results = pd.read_csv('results/{}/CSR/{}/output/mutation_assignments.txt'.format(patient_id, folder),
                              sep='\t', header=None, names=['cluster_id'])
    mut_list = pd.read_csv('results/{}/CSR/{}/output/mutations_list.txt'.format(patient_id, folder),
                              sep='\t', header=None, names=['tmp_id'])
    csr_results = csr_results.assign(sub_id='chr' + mut_list.tmp_id)
    inp_csr = pd.merge(input_df, csr_results, on='sub_id', how='left')
    nb_clones = inp_csr.cluster_id.nunique()
    phi_values = None
    weights = None
    cluster_assignment = inp_csr[['sub_id', 'cluster_id']]
    result_dict['csr_{}'.format(folder)] = [nb_clones, phi_values, weights, cluster_assignment]


df_list = list()
for k1 in result_dict.keys():
    for k2 in result_dict.keys():
        print(k1, k2)
        score1B = score1B_base(result_dict[k1][0], result_dict[k2][0])
        try:
            score1C = score1C_base(result_dict[k1][1], result_dict[k2][1],
                                   result_dict[k1][2], result_dict[k2][2])
        except TypeError:
            score1C = np.nan
        mm = pd.merge(result_dict[k1][3], result_dict[k2][3], on='sub_id', how='outer')
        mm = mm.dropna(subset=['cluster_id_x', 'cluster_id_y'], how='any')
        try:
            score2A = score2A_base(mm.cluster_id_x.values, mm.cluster_id_y.values)
        except ZeroDivisionError:
            score2A = np.nan
        df_list.append([k1, k2, score1B, score1C, score2A])


final_res = pd.DataFrame(df_list, columns=['run1', 'run2', 'score1B', 'score1C', 'score2A']) 
final_res = final_res.assign(patient_id=patient_id)
final_res.to_csv('results/{}/metrics_ith.csv'.format(patient_id), sep='\t', index=False)



