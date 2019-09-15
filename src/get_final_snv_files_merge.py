#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to conduct and exploratory analysis of the new tcga data.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from util_functions import safe_mkdir
import matplotlib.pyplot as plt

import seaborn as sns
import requests


pd.set_option('display.max_columns', 200)
part_number = 2

plt.style.use('seaborn-white')

cancer_locs = ['BRCA', 'HNSC', 'BLCA']

output_path = 'results/mutations_exploratory_analysis'
safe_mkdir(output_path)

# load absolute purity
all_absolute_purity = pd.read_excel('external_data/ncomms9971-s2.xlsx',
                                    skiprows=[0, 1, 2])
all_absolute_purity = all_absolute_purity.assign(
    patient_id=all_absolute_purity['Sample ID'].str[:12])

# load ascat purity
all_ascat_purity = pd.read_csv('data/pancancer/liftover_ASCAT_TCGA/filtered.combined.acf.ploidy.txt', sep='\t')
all_ascat_purity_relevant = all_ascat_purity[all_ascat_purity.cancer_type.isin(cancer_locs)]
all_ascat_purity_samples = [s[:16] for s in all_ascat_purity_relevant.barcodeTumour.unique()]
all_ascat_purity_patients = [s[:12] for s in all_ascat_purity_samples]


# load CNV data
all_cnv_hg38 = pd.read_csv('data/pancancer/liftover_ASCAT_TCGA/cnasHg38.tsv',
                           sep='\t')
all_cnv_hg38 = all_cnv_hg38.assign(patient_id=all_cnv_hg38.participant.str[:12])
all_cnv_hg38 = all_cnv_hg38.assign(weight=all_cnv_hg38.stop-all_cnv_hg38.start)
all_cnv_hg38 = all_cnv_hg38[all_cnv_hg38.type.isin(cancer_locs)]
all_cnv_hg38 = all_cnv_hg38[all_cnv_hg38.participant.isin(all_ascat_purity['name'])]
all_cnv_hg38 = all_cnv_hg38.assign(key=all_cnv_hg38.patient_id+'_'+all_cnv_hg38.chromosome)

patient_list_list = list()
for file_type in ['qc', 'public']:
    useful_final = dict()
    for cancer_loc in cancer_locs:
        # load snv data
        filename = 'tmp/useful_final_{}_{}.csv'.format(file_type, cancer_loc)
        useful_final[cancer_loc] = pd.read_csv(filename, sep='\t', index_col=0)
        useful_final[cancer_loc] = useful_final[cancer_loc]\
            .assign(patient_id=useful_final[cancer_loc].sample_id.str[:12])
        useful_final[cancer_loc] = useful_final[cancer_loc]\
            .assign(cancer_loc=cancer_loc)

        annotations = pd.read_csv('data/{}/annotations.txt'.format(cancer_loc),
                                  sep='\t')
        barcode_list = list()
        for ann_uuid in annotations['id'].tolist():
            r = requests.get('https://api.gdc.cancer.gov/annotations/{uuid}'.format(uuid=ann_uuid))
            barcode_list.append(r.json()['data']['entity_submitter_id'])
        annotations = annotations.assign(**{'Entity Barcode': barcode_list})
        exclude = [b[:12] for b in annotations['Entity Barcode'].unique().tolist()]
        useful_final[cancer_loc] = useful_final[cancer_loc][~useful_final[cancer_loc].patient_id.isin(exclude)]

    useful_final_merge = pd.concat(useful_final.values(), axis=0)

    useful_final_merge = useful_final_merge\
        .assign(sample_id_short=useful_final_merge.sample_id.str[:-12])

    useful_final_merge = useful_final_merge.assign(key=useful_final_merge.patient_id+'_'+useful_final_merge['mutation_id.1'].str.split('_').str[1].str[3:])
    useful_final_merge = useful_final_merge.assign(position=useful_final_merge['mutation_id.1'].str.split('_').str[2].astype(int))

    tmp_useful_final_merge_cnv = pd.merge(useful_final_merge, all_cnv_hg38, left_on='key', right_on='key', how='left')
    tmp_useful_final_merge_cnv = tmp_useful_final_merge_cnv[((tmp_useful_final_merge_cnv.position >= tmp_useful_final_merge_cnv.start) & (tmp_useful_final_merge_cnv.position <= tmp_useful_final_merge_cnv.stop))]
    useful_final_merge_cnv = pd.merge(useful_final_merge,
                                      tmp_useful_final_merge_cnv[['mutation_id.1', 'major', 'minor']],
                                      left_on='mutation_id.1', right_on='mutation_id.1', how='left')
    useful_final_merge_cnv = useful_final_merge_cnv[(useful_final_merge_cnv.sample_id_short.isin(all_ascat_purity_samples)) |
                                                    (~useful_final_merge_cnv.patient_id.isin(all_ascat_purity_patients))]
    useful_final_merge_cnv_purity_tmp = pd.merge(useful_final_merge_cnv,
                                                 all_ascat_purity[['tissue', 'purity', 'ploidy']],
                                                 left_on='patient_id', right_on='tissue', how='left')
    useful_final_merge_cnv_purity = pd.merge(useful_final_merge_cnv_purity_tmp,
                                             all_absolute_purity[['Sample ID', 'CPE']],
                                             left_on='sample_id_short', right_on='Sample ID', how='left')
    useful_final_merge_cnv_purity = useful_final_merge_cnv_purity.assign(vaf=useful_final_merge_cnv_purity.t_alt_count.astype(float)/useful_final_merge_cnv_purity.t_depth.astype(float))
    useful_final_merge_cnv_purity = useful_final_merge_cnv_purity.assign(min_nonzero_cn=useful_final_merge_cnv_purity.apply(lambda x: x['minor'] if (x['minor']>0) and (x['major']<=2) else x['major'], axis=1))
    useful_final_merge_cnv_purity = useful_final_merge_cnv_purity.assign(total_cn=lambda x: x['minor'] + x['major'])
    useful_final_merge_cnv_purity = useful_final_merge_cnv_purity[useful_final_merge_cnv_purity.total_cn > 0]
    useful_final_merge_cnv_purity = useful_final_merge_cnv_purity.assign(vaf_cn=useful_final_merge_cnv_purity.vaf * useful_final_merge_cnv_purity['total_cn'] / useful_final_merge_cnv_purity['min_nonzero_cn'])
    useful_final_merge_cnv_purity = useful_final_merge_cnv_purity.assign(vaf_purity=useful_final_merge_cnv_purity.apply(lambda x: x['vaf']/x.purity * ((1 - x.purity) * 2 + x.purity * x['total_cn']) / x['min_nonzero_cn'], axis=1))
    useful_final_merge_cnv_purity = useful_final_merge_cnv_purity.assign(vaf_CFE=useful_final_merge_cnv_purity.apply(lambda x: x['vaf']/x.CPE * ((1 - x.CPE) * 2 + x.CPE * x['total_cn']) / x['min_nonzero_cn'], axis=1))

    patient_list_list.append(useful_final_merge_cnv_purity.patient_id.unique().tolist())
    print(useful_final_merge_cnv_purity.shape, useful_final_merge_cnv_purity['mutation_id.1'].nunique())
    useful_final_merge_cnv_purity.drop_duplicates(subset=['mutation_id.1'], inplace=True)
    useful_final_merge_cnv_purity.to_csv('tmp/useful_final_{}_merge_cnv_purity.csv'.format(file_type), sep='\t', index=False)

union_patients = list(set(patient_list_list[0]).union(set(patient_list_list[1])))
udf = pd.DataFrame(union_patients, columns=['patient_id'], index=range(1, len(union_patients)+1))
udf.to_csv('official_patient_list.csv', header=False)

# adding absolute copy number
useful_final_merge_cnv_purity = pd.read_csv('tmp/useful_final_qc_merge_cnv_purity.csv', sep='\t')
abs_data = pd.read_csv('data/pancancer/ABSOLUTE_cnasHg38.tsv', sep='\t')
abs_data = abs_data.assign(patient_id=abs_data.Sample.str[:12])
abs_data = abs_data.assign(key=abs_data.patient_id+'_'+abs_data.Chromosome.astype(str))
tmp_useful_final_merge_cnv = pd.merge(useful_final_merge_cnv_purity, abs_data, left_on='key', right_on='key', how='left')
tmp_useful_final_merge_cnv = tmp_useful_final_merge_cnv[((tmp_useful_final_merge_cnv.position >= tmp_useful_final_merge_cnv.Start) & (tmp_useful_final_merge_cnv.position <= tmp_useful_final_merge_cnv.End))]
tmp_useful_final_merge_cnv.to_csv('tmp/useful_final_qc_merge_cnv_purity_absolute.csv', sep='\t', index=False)


