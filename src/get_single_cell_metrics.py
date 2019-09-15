#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to analyse the results of the different ith methods
on the single cell dataset
"""
import pandas as pd
from util_functions import safe_mkdir
import numpy as np
import scipy as sp
import itertools
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from sklearn.metrics.cluster import v_measure_score
from scipy.spatial.distance import cosine, euclidean




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


ground_truth_co5 = pd.read_csv('external_data/CO5-PA-MA-exome.filtered.BSCITE.csv', sep='\t')
ground_truth_co5 = ground_truth_co5.assign(sub_id=ground_truth_co5.Chromosome + '_' + ground_truth_co5.Position.astype(str))
ground_truth_co5 = ground_truth_co5.assign(cluster_bscite=ground_truth_co5.cluster)

ground_truth_co8 = pd.read_csv('external_data/CO8-PA-MA-exome.BSCITE.filtered.csv', sep='\t')
ground_truth_co8 = ground_truth_co8.assign(sub_id=ground_truth_co8.Chromosome + '_' + ground_truth_co8.Position.astype(str))
ground_truth_co8 = ground_truth_co8.assign(cluster_bscite=ground_truth_co8.cluster)

ground_truth_tnbc = pd.read_csv('external_data/full-tnbc.csv', sep='\t')
ground_truth_tnbc = ground_truth_tnbc.assign(full_chrom='chr' + ground_truth_tnbc.Chromosome.astype(str))
ground_truth_tnbc = ground_truth_tnbc.assign(sub_id=ground_truth_tnbc.full_chrom + '_' + ground_truth_tnbc.liftover_pos.astype(str).str.replace('.0', ''))

ground_truth_p1 = pd.read_csv('external_data/gawad_p1_malilik.csv', sep='\t')
ground_truth_p1 = ground_truth_p1.assign(sub_id=ground_truth_p1.Chromosome + '_' + ground_truth_p1.Position.astype(str))
ground_truth_p1 = ground_truth_p1.assign(ID=ground_truth_p1.Gene)

ground_truth_p2 = pd.read_csv('external_data/gawad_p2_malilik.csv', sep='\t')
ground_truth_p2 = ground_truth_p2.assign(sub_id=ground_truth_p2.Chromosome + '_' + ground_truth_p2.Position.astype(str))
ground_truth_p2 = ground_truth_p2.assign(ID=ground_truth_p2.Gene)


patient_ground_truth = {'CO5_colon': ground_truth_co5,
                        'CO5_liver': ground_truth_co5,
                        'CO8_colon': ground_truth_co8,
                        'CO8_liver': ground_truth_co8,
                        'BRCA_wang_TNBC': ground_truth_tnbc,
                        'ALL_gawad_P1': ground_truth_p1,
                        'ALL_gawad_P2': ground_truth_p2}



final_list = list()
for patient in patient_ground_truth.keys():
    res_list = list()
    ground_truth = patient_ground_truth[patient]
    pyclone_results = pd.read_csv('results/{}/pyclone/protected_hg38_vcf/tables/loci.tsv'.format(patient), sep='\t')
    pyclone_results = pyclone_results.assign(sub_id=pyclone_results.mutation_id.str[:-4])
    pyclone_results = pyclone_results.assign(cluster_pyclone=pyclone_results.cluster_id)
    pyclone_results_merge = pd.merge(pyclone_results, ground_truth[['ID', 'sub_id', 'cluster_bscite']], on='sub_id', how='left')
    try:
        sciclone_results = pd.read_csv('results/{}/sciclone/protected_hg38_vcf/clusters1'.format(patient), sep='\t')
        sciclone_results = sciclone_results.assign(full_chrom='chr' + sciclone_results.chr.astype(str))
        sciclone_results = sciclone_results.assign(sub_id=sciclone_results.full_chrom + '_' + sciclone_results.st.astype(str))
        sciclone_results = sciclone_results.assign(cluster_sciclone=sciclone_results.cluster)
        pyclone_results_merge_sc = pd.merge(pyclone_results_merge, sciclone_results[['sub_id', 'cluster_sciclone']], on='sub_id', how='left')
    except:
        pyclone_results_merge_sc = pyclone_results_merge.copy()

    try:
        phylowgs_input = pd.read_csv('results/{}/PhyloWGS/protected_hg38_vcf/parsed_ssm_data.txt'.format(patient), sep='\t')
        phylowgs_input = phylowgs_input.assign(sub_id='chr' + phylowgs_input.gene)
        phylowgs_results = pd.read_csv('results/{}/PhyloWGS/protected_hg38_vcf/outputs/2A.txt'.format(patient), header=None, names=['cluster_phylowgs'])
        phylowgs_input = phylowgs_input.assign(cluster_phylowgs=phylowgs_results.cluster_phylowgs)
        pyclone_results_merge_sc_ph = pd.merge(pyclone_results_merge_sc, phylowgs_input[['sub_id', 'cluster_phylowgs']], on='sub_id', how='left')
    except:
        pyclone_results_merge_sc_ph = pyclone_results_merge_sc.copy()

    try:
        expands_results = pd.read_csv('results/{}/expands/protected_hg38_vcf/{}__protected_hg38_vcf.sps'.format(patient, patient), sep='\t', skiprows=1)
        expands_results = expands_results.assign(full_chrom='chr' + expands_results.chr.astype(str))
        expands_results = expands_results.assign(sub_id=expands_results.full_chrom + '_' + expands_results.startpos.astype(str))
        expands_results = expands_results.assign(cluster_expands=expands_results.SP.apply(lambda x: expands_results.SP.unique().tolist().index(x)))
        pyclone_results_merge_sc_ph_ex = pd.merge(pyclone_results_merge_sc_ph, expands_results[['sub_id', 'cluster_expands']], on='sub_id', how='left')
    except:
        pyclone_results_merge_sc_ph_ex = pyclone_results_merge_sc_ph.copy()
    try:
        csr_results = pd.read_csv('results/{}/CSR/protected_hg38_vcf/output/mutation_assignments.txt'.format(patient),
                                  sep='\t', header=None, names=['cluster_csr'])
        mut_list = pd.read_csv('results/{}/CSR/protected_hg38_vcf/output/mutations_list.txt'.format(patient),
                                  sep='\t', header=None, names=['tmp_id'])
        csr_results = csr_results.assign(sub_id='chr' + mut_list.tmp_id)
        pyclone_results_merge_sc_ph_ex_csr = pd.merge(pyclone_results_merge_sc_ph_ex, csr_results[['sub_id', 'cluster_csr']], on='sub_id', how='left')
    except:
        pyclone_results_merge_sc_ph_ex_csr = pyclone_results_merge_sc_ph_ex.copy()

    res_list.append(patient)
    res_list.append(ground_truth.cluster_bscite.nunique())
    res_list.append(pyclone_results.cluster_pyclone.nunique())
    if 'cluster_sciclone' in pyclone_results_merge_sc_ph_ex.columns:
        res_list.append(pyclone_results_merge_sc_ph_ex_csr.cluster_sciclone.nunique())
    else:
        res_list.append(np.nan)
    if 'cluster_phylowgs' in pyclone_results_merge_sc_ph_ex_csr.columns:
        res_list.append(pyclone_results_merge_sc_ph_ex_csr.cluster_phylowgs.nunique())
    else:
        res_list.append(np.nan)
    if 'cluster_expands' in pyclone_results_merge_sc_ph_ex_csr.columns:
        res_list.append(pyclone_results_merge_sc_ph_ex_csr.cluster_expands.nunique())
    else:
        res_list.append(np.nan)
    if 'cluster_csr' in pyclone_results_merge_sc_ph_ex_csr.columns:
        res_list.append(pyclone_results_merge_sc_ph_ex_csr.cluster_csr.nunique())
    else:
        res_list.append(np.nan)    
    score2A_df = pyclone_results_merge_sc_ph_ex_csr.dropna(subset=['cluster_bscite'])
    res_list.append(score2A_base(score2A_df.cluster_bscite, score2A_df.cluster_pyclone))
    if 'cluster_sciclone' in pyclone_results_merge_sc_ph_ex_csr.columns:
        res_list.append(score2A_base(score2A_df.cluster_bscite, score2A_df.cluster_sciclone))
    else:
        res_list.append(np.nan)
    if 'cluster_phylowgs' in pyclone_results_merge_sc_ph_ex_csr.columns:
        res_list.append(score2A_base(score2A_df.cluster_bscite, score2A_df.cluster_phylowgs))
    else:
        res_list.append(np.nan)
    if 'cluster_expands' in pyclone_results_merge_sc_ph_ex_csr.columns:
        res_list.append(score2A_base(score2A_df.cluster_bscite, score2A_df.cluster_expands))
    else:
        res_list.append(np.nan)
    if 'cluster_csr' in pyclone_results_merge_sc_ph_ex_csr.columns:
        res_list.append(score2A_base(score2A_df.cluster_bscite, score2A_df.cluster_csr))
    else:
        res_list.append(np.nan)
    res_list.append(100 * pyclone_results_merge_sc_ph_ex_csr.variant_allele_frequency.mad() / pyclone_results_merge_sc_ph_ex_csr.variant_allele_frequency.median())
    res_list.append(len(pyclone_results_merge_sc_ph_ex_csr))
    res_list.append(len(score2A_df))
    final_list.append(res_list)


res_df = pd.DataFrame(final_list, columns=['sample', 'nb_clones_bscite',
                                           'nb_clones_pyclone', 'nb_clones_sciclone',
                                           'nb_clones_phylowgs', 'nb_clones_expands',
                                           'nb_clones_csr',
                                           'score2A_pyclone', 'score2A_sciclone',
                                           'score2A_phylowgs', 'score2A_expands',
                                           'score2A_csr', 'MATH_score',
                                           'nb_mutations_ith', 'nb_mutations_metric'])

print(res_df.to_latex(float_format=lambda x: '%.3f' % x)) 


