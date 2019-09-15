#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to gather all ITH results in a final exploitable table.
An initial table is built by using clinical data, and other precomputed results, such
as purity, copy number.
ITH features are added by going over all result files.
Run times are extracted for each patient from torque logs, using file system timestamps.
Immune features are then extracted from RNAseq results.
"""
import os
import sys
import pandas as pd
import numpy as np


pd.set_option('display.max_columns', 200)

CANCER_LOCS = ['BRCA', 'HNSC', 'BLCA']
TORQUE_LOG_PATH_old = sys.argv[1]
TORQUE_LOG_PATH_new = sys.argv[2]

if __name__ == '__main__':
    ###################
    ## initial table ##
    ###################
    clinical_data_dict = dict()
    for cancer_loc in CANCER_LOCS:
        # load clinical data
        clinical_data = pd.read_csv(
            'data/{}/clinical/data_bcr_clinical_data_patient.txt'.format(cancer_loc),
            sep='\t', skiprows=4)
        clinical_data_sample = pd.read_csv(
            'data/{}/clinical/data_bcr_clinical_data_sample.txt'.format(cancer_loc),
            sep='\t', skiprows=4)
        clinical_data = clinical_data.assign(
            binary_vital_status=clinical_data.apply(
                lambda x: 1 if x.OS_STATUS == 'DECEASED' else
                (0 if x.OS_STATUS == 'LIVING' else np.nan), axis=1))
        clinical_data = clinical_data[clinical_data.OS_MONTHS != '[Not Available]']
        clinical_data = clinical_data[clinical_data.DAYS_TO_BIRTH != '[Not Available]']
        clinical_data = clinical_data.assign(survival_days=clinical_data.apply(lambda x: 30.5*float(x['OS_MONTHS']), axis=1))
        clinical_data = clinical_data.assign(age_at_diagnosis=clinical_data.apply(lambda x: int(np.round(-float(x.DAYS_TO_BIRTH)/365)), axis=1))
        clinical_data = clinical_data.assign(cancer_loc=cancer_loc)
        clinical_data_dict[cancer_loc] = clinical_data

    useful_final_qc_merge_purity = pd.read_csv('tmp/useful_final_qc_merge_cnv_purity_absolute.csv', sep='\t')
    useful_final_public_merge_purity = pd.read_csv('tmp/useful_final_public_merge_cnv_purity.csv', sep='\t')
    common_cols = ['PATIENT_ID', 'binary_vital_status', 'survival_days', 'age_at_diagnosis', 'cancer_loc']
    clinical_data_merge = pd.concat([clinical_data_dict[loc][common_cols] for loc in CANCER_LOCS], axis=0)

    # load absolute purity
    abs_purity = pd.read_csv(
        'data/pancancer/TCGA_mastercalls.abs_tables_JSedit.fixed.txt',
        sep='\t')
    abs_purity = abs_purity.assign(patient_id=abs_purity.array.str[:12])
    abs_purity = abs_purity.assign(sample_id=abs_purity['sample'].str[:16])

    # load ascat purity
    all_ascat_purity = pd.read_csv('data/pancancer/liftover_ASCAT_TCGA/filtered.combined.acf.ploidy.txt', sep='\t')
    all_ascat_purity_relevant = all_ascat_purity[all_ascat_purity.cancer_type.isin(CANCER_LOCS)]
    all_ascat_purity_samples = [s[:16] for s in all_ascat_purity_relevant.barcodeTumour.unique()]
    all_ascat_purity_patients = [s[:12] for s in all_ascat_purity_samples]

    # load ascat CNV data
    all_cnv_hg38 = pd.read_csv('data/pancancer/liftover_ASCAT_TCGA/cnasHg38.tsv', sep='\t')
    all_cnv_hg38 = all_cnv_hg38[all_cnv_hg38.participant.isin(all_ascat_purity['name'])]
    all_cnv_hg38 = all_cnv_hg38.assign(patient_id=all_cnv_hg38.participant.str[:12])
    all_cnv_hg38 = all_cnv_hg38.assign(weight=all_cnv_hg38.stop-all_cnv_hg38.start)

    # load absolute CNV data
    abs_data = pd.read_csv('data/pancancer/ABSOLUTE_cnasHg38.tsv', sep='\t')
    abs_data = abs_data.assign(patient_id=abs_data.Sample.str[:12])
    abs_data = abs_data.assign(key=abs_data.patient_id+'_'+abs_data.Chromosome.astype(str))
    abs_data = abs_data.assign(weight=abs_data.End-abs_data.Start)

    useful_final_merge_cnv_purity_abs = pd.merge(useful_final_qc_merge_purity, abs_purity[['sample_id', 'purity']], left_on='Sample ID', right_on='sample_id', how='left', suffixes=['', '_absolute'])
    useful_final_merge_cnv_purity_abs.purity_absolute = useful_final_merge_cnv_purity_abs.purity_absolute.fillna(1)
    useful_final_merge_cnv_purity_abs = useful_final_merge_cnv_purity_abs.assign(min_nonzero_cn_abs=useful_final_merge_cnv_purity_abs.apply(lambda x: x['Modal_HSCN_1'] if (x['Modal_HSCN_1']>0) and (x['Modal_HSCN_2']<=2) else x['Modal_HSCN_2'] if x['Modal_HSCN_2']>0 else x['Modal_Total_CN'], axis=1))
    useful_final_merge_cnv_purity_abs = useful_final_merge_cnv_purity_abs.assign(total_cn_abs=useful_final_merge_cnv_purity_abs.Modal_Total_CN)
    useful_final_merge_cnv_purity_abs = useful_final_merge_cnv_purity_abs[useful_final_merge_cnv_purity_abs.total_cn_abs > 0]
    useful_final_merge_cnv_purity_abs = useful_final_merge_cnv_purity_abs.assign(vaf_cn_abs=useful_final_merge_cnv_purity_abs.vaf * useful_final_merge_cnv_purity_abs['total_cn_abs'] / useful_final_merge_cnv_purity_abs['min_nonzero_cn_abs'])
    useful_final_merge_cnv_purity_abs = useful_final_merge_cnv_purity_abs.assign(vaf_purity_abs=useful_final_merge_cnv_purity_abs.apply(lambda x: x['vaf']/x.purity * ((1 - x.purity_absolute) * 2 + x.purity_absolute * x['total_cn_abs']) / x['min_nonzero_cn_abs'], axis=1))

    def key_barcode(b):
        return int(b[13:15])

    non_diploid_fraction = all_cnv_hg38[(all_cnv_hg38.major!=1) | (all_cnv_hg38.minor!=1)].groupby('patient_id').weight.sum()/ all_cnv_hg38.groupby('patient_id').weight.sum()
    non_diploid_fraction_absolute = abs_data[(abs_data.Modal_HSCN_2!=1) | (abs_data.Modal_HSCN_1!=1)].groupby('patient_id').weight.sum()/ abs_data.groupby('patient_id').weight.sum()
    clinical_data_merge_ascat = pd.merge(clinical_data_merge, all_ascat_purity_relevant, left_on='PATIENT_ID', right_on='tissue', how='left')
    clinical_data_merge_ascat = clinical_data_merge_ascat.assign(sample_id=clinical_data_merge_ascat.barcodeTumour.str[:16])

    clinical_data_merge_ascat_abs = pd.merge(clinical_data_merge_ascat, abs_purity[['sample_id', 'purity', 'ploidy', 'Subclonal genome fraction']], left_on='sample_id', right_on='sample_id', how='left', suffixes=['', '_absolute'])
    clinical_data_merge_ascat_abs_cnv = pd.merge(clinical_data_merge_ascat_abs, non_diploid_fraction.to_frame(), left_on='PATIENT_ID', right_index=True, how='left')
    clinical_data_merge_ascat_abs_cnv_abs = pd.merge(clinical_data_merge_ascat_abs_cnv, non_diploid_fraction_absolute.to_frame(), left_on='PATIENT_ID', right_index=True, how='left', suffixes=['', '_absolute'])

    clinical_data_merge_ascat_abs_cnv_protect = pd.merge(clinical_data_merge_ascat_abs_cnv_abs, useful_final_qc_merge_purity.groupby('patient_id_x')['mutation_id.1'].count().to_frame(), left_on='PATIENT_ID', right_index=True, how='inner')
    clinical_data_merge_ascat_abs_cnv_protect_pub = pd.merge(clinical_data_merge_ascat_abs_cnv_protect, useful_final_public_merge_purity.groupby('patient_id')['mutation_id.1'].count().to_frame(), left_on='PATIENT_ID', right_index=True, how='left')


    clinical_data_merge_ascat_abs_cnv_protect_pub.rename(columns={'mutation_id.1_x': 'mutation_count_protected',
                                                                  'mutation_id.1_y': 'mutation_count_public', 'weight': 'perc_non_diploid_ascat',
                                                                  'weight_absolute': 'perc_non_diploid_absolute'}, inplace=True)

    ######################
    ## add ITH features ##
    ######################
    # compute MATH score
    hh = useful_final_merge_cnv_purity_abs.groupby('patient_id_x')
    math_protected = (100 * hh.vaf.mad() / hh.vaf.median()).to_frame()
    math_protected_cn = (100 * hh.vaf_cn.mad() / hh.vaf_cn.median()).to_frame()
    math_protected_cn_abs = (100 * hh.vaf_cn_abs.mad() / hh.vaf_cn_abs.median()).to_frame()

    gg = useful_final_public_merge_purity.groupby('patient_id')
    math_public = (100 * gg.vaf.mad() / gg.vaf.median()).to_frame()
    math_public_cn = (100 * gg.vaf_cn.mad() / gg.vaf_cn.median()).to_frame()

    for df in (math_protected, math_protected_cn, math_protected_cn_abs, math_public, math_public_cn):
        clinical_data_merge_ascat_abs_cnv_protect_pub = pd.merge(clinical_data_merge_ascat_abs_cnv_protect_pub, df, left_on='PATIENT_ID', right_index=True, how='left')

    clinical_data_merge_ascat_abs_cnv_protect_pub.rename(columns={'vaf_x': 'math_protected',
                                                                  'vaf_y': 'math_public',
                                                                  'vaf_cn_x': 'math_cn_protected',
                                                                  'vaf_cn_y': 'math_cn_public',
                                                                  'vaf_cn_abs': 'math_cn_protected_abs'}, inplace=True)
    useful_final_qc_merge_purity = useful_final_qc_merge_purity.assign(patient_id=useful_final_qc_merge_purity.patient_id_y)
    # get result for other methods
    def get_ith_method_measures(ith_method, folder, patient):
        sub_protected = useful_final_qc_merge_purity[useful_final_qc_merge_purity.patient_id == patient]
        sub_protected = sub_protected.assign(mutation_id=sub_protected['mutation_id.1'].str[29:-2] + sub_protected['mutation_id.1'].str[-1])
        sub_protected.index = range(1, len(sub_protected)+1)
        if ith_method == 'sciclone':
            try:
                sci = pd.read_csv('results/{}/sciclone/{}/clusters1'.format(patient, folder), sep='\t')
                h = pd.DataFrame(sci.cluster.value_counts())
                sci['cluster_size'] = pd.merge(sci, h, left_on='cluster', right_index=True)[['cluster_y']]
                sci = sci.assign(cluster_name=lambda x: x['cluster'])
                sci = sci[(sci.cluster_size > 1) & ( ~pd.isnull(sci.cluster))]
                nb_clones = len(sci.cluster.unique())
                zz = sci.groupby('cluster')[['tumor.vaf', 'cluster_size']].mean()
                clonal_prop = float(zz[zz['tumor.vaf'] == zz['tumor.vaf'].max()].
                                    cluster_size/zz.cluster_size.sum()*100)
                smallest_vaf = zz['tumor.vaf'].min()/100
                zz = zz.assign(prop=zz.cluster_size/zz.cluster_size.sum())
                shannon_index = -np.sum(zz.prop*np.log(zz.prop))
                most_populated_clone_vaf = zz[zz['cluster_size'] == zz['cluster_size'].max()]['tumor.vaf'].values[0]/100
            except IOError:
                nb_clones, clonal_prop, smallest_vaf, shannon_index, most_populated_clone_vaf = [np.nan]*5

        elif ith_method == 'pyclone':
            try:
                pyclone = pd.read_csv('results/{}/pyclone/{}/tables/cluster.tsv'.format(patient, folder), sep='\t')
                if pyclone['size'].max() > 5:
                    pyclone = pyclone[pyclone['size'] > 5]
                else:
                    pyclone = pyclone[pyclone['size'] >= pyclone['size'].max()]
                nb_clones = len(pyclone)
                clonal_prop = float(pyclone[pyclone['mean'] == pyclone['mean'].max()]['size'].astype(float)/pyclone['size'].sum()*100)
                smallest_vaf = pyclone['mean'].min()
                pyclone = pyclone.assign(prop=pyclone['size']/pyclone['size'].sum())
                shannon_index = -np.sum(pyclone.prop*np.log(pyclone.prop))
                most_populated_clone_vaf = pyclone[pyclone['size']==pyclone['size'].max()]['mean'].values[0]
            except IOError:
                nb_clones, clonal_prop, smallest_vaf, shannon_index, most_populated_clone_vaf = [np.nan]*5

        elif ith_method == 'PhyloWGS':
            try:
                m = pd.read_csv('results/{}/{}/{}/outputs/1C.txt'.format(patient, ith_method, folder), sep='\t', names=['cluster_id', 'size', 'mean'], header=None)
                if m['size'].max() > 5:
                    m = m[m['size'] > 5]
                else:
                    m = m[m['size'] >= m['size'].max()]
                nb_clones = len(m)
                clonal_prop = float(m[m['mean'] == m['mean'].max()]['size'].astype(float)/m['size'].sum()*100)
                smallest_vaf = m['mean'].min()
                m = m.assign(prop=m['size']/m['size'].sum())
                shannon_index = -np.sum(m.prop*np.log(m.prop))
                most_populated_clone_vaf = m[m['size'] == m['size'].max()]['mean'].values[0]
            except IOError:
                nb_clones, clonal_prop, smallest_vaf, shannon_index, most_populated_clone_vaf = [np.nan]*5

        elif ith_method == 'baseline':
            try:
                m = pd.read_csv('results/{}/{}/{}/cluster_assignment_bis.csv'.format(patient, ith_method, folder), sep='\t')
                if m['cluster_size'].max() > 5:
                    m = m[m['cluster_size'] > 5]
                else:
                    m = m[m['cluster_size'] >= m['cluster_size'].max()]
                nb_clones = len(m)
                clonal_prop = float(m[m['mean'] == m['mean'].max()]['cluster_size'].astype(float)/m['cluster_size'].sum()*100)
                smallest_vaf = m['mean'].min()
                m = m.assign(prop=m['cluster_size']/m['cluster_size'].sum())
                shannon_index = -np.sum(m.prop*np.log(m.prop))
                most_populated_clone_vaf = m[m['cluster_size'] == m['cluster_size'].max()]['mean'].values[0]
            except IOError:
                nb_clones, clonal_prop, smallest_vaf, shannon_index, most_populated_clone_vaf = [np.nan]*5

        elif ith_method == 'expands':
            try:
                all_mut = pd.read_csv('results/{}/{}/{}/{}__{}.sps'.format(patient, ith_method, folder, patient, folder), sep='\t', skiprows=1)
                m = all_mut.SP.value_counts().to_frame()
                m = m.assign(cluster_size=m.SP)
                m = m.assign(mean=m.index)
                if m['cluster_size'].max() > 5:
                    m = m[m['cluster_size'] > 5]
                else:
                    m = m[m['cluster_size'] >= m['cluster_size'].max()]
                nb_clones = len(m)
                clonal_prop = float(m[m['mean'] == m['mean'].max()]['cluster_size'].astype(float)/m['cluster_size'].sum()*100)
                smallest_vaf = m['mean'].min()
                m = m.assign(prop=m['cluster_size']/m['cluster_size'].sum())
                shannon_index = - np.sum(m.prop*np.log(m.prop))
                most_populated_clone_vaf = m[m['cluster_size'] == m['cluster_size'].max()]['mean'].values[0]
            except IOError:
                nb_clones, clonal_prop, smallest_vaf, shannon_index, most_populated_clone_vaf = [np.nan]*5

        elif ith_method == 'CSR':
            try:
                summary = pd.read_csv('results/{}/{}/{}/output/summary_table.txt'.format(patient, ith_method, folder),
                                      sep='\t', header=None, names=['cluster_name', 'nb_mut', 'cluster_cp'])
                if summary['nb_mut'].max() > 5:
                    summary = summary[summary['nb_mut'] > 5]
                else:
                    summary = summary[summary['nb_mut'] >= summary['nb_mut'].max()]
                nb_clones = len(summary)
                clonal_prop = float(summary[summary['cluster_cp'] == summary['cluster_cp'].max()]['nb_mut'].astype(float)/summary['nb_mut'].sum()*100)
                smallest_vaf = summary['cluster_cp'].min()
                summary = summary.assign(prop=summary['nb_mut']/summary['nb_mut'].sum())
                shannon_index = - np.sum(summary.prop*np.log(summary.prop))
                most_populated_clone_vaf = summary[summary['nb_mut'] == summary['nb_mut'].max()]['cluster_cp'].values[0]
            except IOError:
                nb_clones, clonal_prop, smallest_vaf, shannon_index, most_populated_clone_vaf = [np.nan]*5
        return patient, nb_clones, clonal_prop, smallest_vaf, shannon_index, most_populated_clone_vaf

    tmp_df = clinical_data_merge_ascat_abs_cnv_protect_pub.copy()
    base_cols = ['nb_clones', 'clonal_prop', 'smallest_vaf', 'shannon_index', 'most_populated_clone_vaf']
    for ith_method in ['pyclone', 'PhyloWGS', 'sciclone', 'expands', 'CSR']:
        for folder in ['protected_hg38_vcf', 'public_hg38_vcf', 'protected_hg38_vcf_absolute']:
            col_list = ['PATIENT_ID'] + ['{}_{}_{}'.format(ith_method, folder, c) for c in base_cols]
            pilot_results = pd.DataFrame(columns=col_list)
            for patient in clinical_data_merge_ascat_abs_cnv_protect_pub.PATIENT_ID.tolist():
                out = get_ith_method_measures(ith_method, folder, patient)
                pilot_results = pilot_results.append(dict(zip(col_list, out)), ignore_index=True)
                # smg to build pilot_results
            tmp_df = pd.merge(tmp_df, pilot_results, left_on='PATIENT_ID', right_on='PATIENT_ID', how='left')
    tmp_df.to_csv('20190901_tmp_df_get_features.csv', sep='\t', index=False)
    tmp_df = pd.read_csv('20190901_tmp_df_get_features.csv', sep='\t')


    ###################
    ## add run times ##
    ###################

    def get_timestamp(ith_method, folder, patient):
        if ith_method=='sciclone':
            try:
                t = os.path.getmtime('results/{}/sciclone/{}/clusters1'.format(patient, folder))
            except IOError:
                t = np.nan

        elif ith_method=='pyclone':
            try:
                t = os.path.getmtime('results/{}/pyclone/{}/tables/cluster.tsv'.format(patient, folder))

            except IOError:
                t = np.nan

        elif ith_method == 'PhyloWGS':
            try:
                t = os.path.getmtime('results/{}/{}/{}/outputs/1C.txt'.format(patient, ith_method, folder))
            except IOError:
                t = np.nan

        elif ith_method == 'baseline':
            try:
                t = os.path.getmtime('results/{}/{}/{}/cluster_assignment_for_each_mutation.csv'.format(patient, ith_method, folder))
            except IOError:
                t = np.nan

        elif ith_method == 'expands':
            try:
                t = os.path.getmtime('results/{}/{}/{}/{}__{}.sps'.format(patient, ith_method, folder, patient, folder))
            except IOError:
                t = np.nan
        return t

    for ith_method in ['pyclone', 'PhyloWGS', 'sciclone', 'expands']:
        for folder in ['protected_hg38_vcf', 'public_hg38_vcf', 'protected_hg38_vcf_absolute']:
            out_list = list()
            for patient in tmp_df.PATIENT_ID.tolist():
                out = get_timestamp(ith_method, folder, patient)
                out_list.append(out)
            print(ith_method, folder)
            tmp_df = tmp_df.assign(**{'{}_{}_timestamp'.format(ith_method, folder): out_list})
    match_numbers_1 = pd.read_csv('official_patient_list.csv', sep=',', header=None, index_col=0, names=['patient_id'])
    match_numbers_1 = match_numbers_1.assign(line_nb=match_numbers_1.index)
    tmp_df_m = pd.merge(tmp_df, match_numbers_1, left_on='PATIENT_ID', right_on='patient_id', how='left')
    tmp_df_m.to_csv('20190902_tmp_df_get_features_timestamp.csv', sep='\t', index=False)
    tmp_df_m = pd.read_csv('20190902_tmp_df_get_features_timestamp.csv', sep='\t')

    # get big table of torque logs
    all_files = os.listdir(TORQUE_LOG_PATH_old)
    all_files_new = os.listdir(TORQUE_LOG_PATH_new)
    all_files_e_old = [f for f in all_files if '.E' in f]
    all_files_e_new = [f for f in all_files_new if '.E' in f]
    tables_e_list = list()
    for n in all_files_e_old:
        if n in all_files_e_old:
            e = pd.read_csv('{}/{}'.format(TORQUE_LOG_PATH_old, n), sep=' ', header=None)
        else:
            e = pd.read_csv('{}/{}'.format(TORQUE_LOG_PATH_new, n), sep=' ', header=None)

        e = e.assign(start=e[8].str.replace('start=', '').str.replace('start_count=', '').astype(int))
        if e.shape[1]-1==24:
            e = e.assign(end=e[18].str.replace('end=', '').astype(int))
            e = e.assign(exit_status=e[19].str.replace('Exit_status=', '').astype(int))
        elif e.shape[1]-1==23:
            e = e.assign(end=e[17].str.replace('end=', '').astype(int))
            e = e.assign(exit_status=e[18].str.replace('Exit_status=', '').astype(int))
        elif e.shape[1]-1==27:
            e = e.assign(end=e[20].str.replace('end=', '').astype(int))
            e = e.assign(exit_status=e[21].str.replace('Exit_status=', '').astype(int))
        else:
            print('pb', n)
            break
        e = e.assign(job_id=e[1].str.split(';').str[2].str.replace('.torque.curie.fr', '').str.replace('.torque6.curie.fr', ''))
        e = e.assign(jobname=e[3].str.replace('jobname=', ''))
        e = e.assign(job_nb=e.job_id.str.split('[').str[0])
        e = e[~pd.isnull(e.job_id.str.split('[').str[1])]
        e = e.assign(batch_nb=e.job_id.str.split('[').str[1].str.replace(']', '').astype(int))
        tables_e_list.append(e[['job_id', 'jobname', 'start', 'end', 'exit_status', 'job_nb', 'batch_nb']])
    big_e = pd.concat(tables_e_list, axis=0)
    big_e = big_e.assign(runtime=big_e.end - big_e.start)

    # merge with the timestamps table
    tolerance = 2  # to allow some latency
    tmp_merge = pd.merge(tmp_df_m, big_e, left_on='line_nb', right_on='batch_nb', how='left')
    for ith_method in ['pyclone', 'sciclone', 'expands', 'PhyloWGS']:
        # phylowgs is a special case
        for folder in ['protected_hg38_vcf', 'public_hg38_vcf']:
            tmp_ = tmp_merge[(tmp_merge['{}_{}_timestamp'.format(ith_method, folder)] >= tmp_merge.start)&(tmp_merge['{}_{}_timestamp'.format(ith_method, folder)]<=tmp_merge.end + tolerance)]
            maj_job_nb = tmp_.job_nb.value_counts().index[0]
            print(ith_method, folder, maj_job_nb, tmp_.job_nb.value_counts())
            tmp_ = tmp_merge[(tmp_merge['{}_{}_timestamp'.format(ith_method, folder)] >= tmp_merge.start) &
                             (tmp_merge['{}_{}_timestamp'.format(ith_method, folder)] <= tmp_merge.end + tolerance) &
                             (tmp_merge['job_nb'] == maj_job_nb)]
            tmp_ = tmp_.assign(**{'{}_{}_runtime'.format(ith_method, folder): tmp_.runtime})
            tmp_df_m = pd.merge(tmp_df_m, tmp_[['PATIENT_ID', '{}_{}_runtime'.format(ith_method, folder)]], left_on='PATIENT_ID', right_on='PATIENT_ID', how='left')



    # manually remove the cases that were allowed to run for too long by the cluster
    for ith_method in ['pyclone', 'sciclone', 'expands', 'PhyloWGS']:
        # phylowgs is a special case
        for folder in ['protected_hg38_vcf', 'public_hg38_vcf']:
            tmp_df_m.loc[tmp_df_m['{}_{}_runtime'.format(ith_method, folder)] > 15*3600, '{}_{}_nb_clones'.format(ith_method, folder)] = np.nan


    ################################
    ## add RNAseq immune features ##
    ################################

    bindea_entrez = pd.read_csv('external_data/bindea13_signatures_entrez.csv', sep=',')
    # adapted from the supplementary of the article
    ensembl_entrz = pd.read_csv('external_data/ensembl_entrez.txt', sep='\t')
    # downloaded from ensembl biomart version 91 on 03/26/2018
    new_bind = pd.merge(bindea_entrez, ensembl_entrz, left_on='EntrezGene', right_on='NCBI gene ID', how='left')
    new_bind.to_csv('data/pancancer/bindea_signatures_ensembl.csv', sep='\t', index=False)

    signatures_definition = pd.read_csv('data/pancancer/bindea_signatures_ensembl.csv', sep='\t')
    sig_list = list()
    for loc in CANCER_LOCS:
        rnaseq_data = pd.read_csv('data/{}/RNAseq/tcga_vst_matrix.csv'.format(loc), sep='\t')
        rna_seq_id = rnaseq_data.index.tolist()
        intersetc = set(rna_seq_id).intersection(set(signatures_definition['Gene stable ID']))
        print(signatures_definition[signatures_definition['Gene stable ID'].isin(intersetc)].shape)
        signatures_tcga = pd.merge(signatures_definition[['CellType', 'Gene stable ID']], rnaseq_data, left_on='Gene stable ID', right_index=True, how='left')
        signatures_tcga.drop(['Gene stable ID'], axis=1, inplace=True)
        signatures_tcga_agg = signatures_tcga.groupby('CellType').mean().T
        signatures_tcga_agg.index=[i[:16] for i in signatures_tcga_agg.index]
        sig_list.append(signatures_tcga_agg)
    all_sig = pd.concat(sig_list, axis=0)
    tmp_df_sig = pd.merge(tmp_df_m, all_sig, left_on='sample_id', right_index=True, how='left').drop_duplicates(subset=['barcodeTumour'])
    tmp_df_sig.to_csv('tmp/20190903_ith_method_metrics_final_runtime_immunity.csv', sep='\t', index=False)


