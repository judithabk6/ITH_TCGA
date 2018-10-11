#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to get pyclone input for each
TCGA patient
"""

import sys
import pandas as pd
import numpy as np
from util_functions import safe_mkdir


p = sys.argv[1]
folder = sys.argv[2]
input_filename = sys.argv[3]

useful_final_merge_cnv_purity = pd.read_csv(input_filename, sep='\t')
pancancer_path = 'data/pancancer/liftover_ASCAT_TCGA'
all_ascat_purity = pd.read_csv('{}/filtered.combined.acf.ploidy.txt'
                               .format(pancancer_path), sep='\t')

all_cnv_hg38 = pd.read_csv('{}/cnasHg38.tsv'.format(pancancer_path), sep='\t')
# deal with the fact that there may be several samples per patient, choose
# the one where purity is available to ensure consistency
all_cnv_hg38 = all_cnv_hg38[all_cnv_hg38.participant.isin(all_ascat_purity['name'])]
all_cnv_hg38 = all_cnv_hg38.assign(patient_id=all_cnv_hg38.participant.str[:12])
all_cnv_hg38 = all_cnv_hg38.assign(weight=all_cnv_hg38.stop-all_cnv_hg38.start)
all_cnv_hg38 = all_cnv_hg38[all_cnv_hg38.weight > 0]


pyclone_fields = ['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn']
cnv_cols = ['chromosome', 'start', 'stop', 'major', 'minor', 'patient_id']

sub_protected = useful_final_merge_cnv_purity[useful_final_merge_cnv_purity.patient_id == p]
sub_protected = sub_protected.assign(ref_counts=sub_protected['t_ref_count'].astype(int))
sub_protected = sub_protected.assign(var_counts=sub_protected['t_alt_count'].astype(int))
sub_protected = sub_protected.assign(minor_cn=sub_protected['minor'].astype(int))
sub_protected = sub_protected.assign(major_cn=sub_protected['major'].astype(int))
sub_protected = sub_protected.assign(normal_cn=2)
sub_protected = sub_protected.assign(mutation_id=sub_protected['mutation_id.1']
                                     .str[29:-2] +
                                     sub_protected['mutation_id.1'].str[-1])

safe_mkdir('results/{}'.format(p))
safe_mkdir('results/{}/pyclone'.format(p))
safe_mkdir('results/{}/pyclone/{}'.format(p, folder))

sub_protected[pyclone_fields].to_csv('results/{}/pyclone/{}/input.tsv'
                                     .format(p, folder), sep='\t', index=False)
sub_cnv = all_cnv_hg38[all_cnv_hg38.patient_id == p][cnv_cols]
sub_cnv.to_csv('results/{}/pyclone/ascat_hg38_cnv_table.csv'.format(p),
               sep='\t', index=False)


patient_purity = sub_protected.purity.unique()[0]
if np.isnan(patient_purity):
    patient_purity = 1.0

with open('results/{}/pyclone/ascat_purity.txt'.format(p), 'w') as f:
    f.write('{}'.format(patient_purity))
