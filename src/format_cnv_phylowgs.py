#!/usr/bin/env python
# -*- coding:utf-8 -*-

# reformat copy number tcga for phyloWGS
import pandas as pd
import sys
from util_functions import safe_mkdir

output_fields = ['chromosome', 'start', 'end', 'copy_number', 'minor_cn', 'major_cn', 'cellular_prevalence']
# no idea of anything so defining
# minor cn as min between 1 and max(copy number - 1, 0)
# cellular prevalence is approximated as purity
patient = sys.argv[1]
# patient = "TCGA-BH-A18U"  # for testing purpose
for method in ('ascat', 'absolute'):
    cnv_table = pd.read_csv('results/{}/pyclone/{}_hg38_cnv_table.csv'.format(patient, method), sep='\t')
    with open('results/{}/pyclone/{}_purity.txt'.format(patient, method), 'r') as f:
        purity = float(f.read())
    cnv_table = cnv_table.assign(start=lambda x: x.start)
    cnv_table = cnv_table.assign(end=lambda x: x.stop)
    cnv_table = cnv_table.assign(copy_number=lambda x: x.major+x.minor)
    cnv_table = cnv_table.assign(cellular_prevalence=lambda x: purity)
    cnv_table = cnv_table.assign(major_cn= cnv_table.major)
    cnv_table = cnv_table.assign(minor_cn= cnv_table.minor)

    safe_mkdir('results/{}/PhyloWGS'.format(patient))
    cnv_table[output_fields].to_csv('results/{}/PhyloWGS/{}_cnv_table.csv'.format(patient, method), sep='\t')

