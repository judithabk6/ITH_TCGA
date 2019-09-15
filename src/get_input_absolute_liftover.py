#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to get input in the right format for segment
liftover.
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')


pd.set_option('display.max_columns', 200)
part_number = 2


absolute_input = pd.read_csv(
    'external_data/TCGA_mastercalls.abs_segtabs.fixed.txt', sep='\t')
absolute_input = absolute_input.dropna(subset=['Chromosome'])
absolute_input = absolute_input.assign(
    Chromosome=absolute_input.Chromosome.astype(int))

absolute_input[['Sample', 'Chromosome', 'Start', 'End', 'Modal_HSCN_1',
                'Modal_HSCN_2', 'Modal_Total_CN']]\
    .to_csv(
        'external_data/TCGA_mastercalls.abs_segtabs.fixed.liftover_input.txt',
        sep='\t', index=False)
