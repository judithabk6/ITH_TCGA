#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to gather the results of alternative metrics 
for the TCGA samples
"""

import pandas as pd
import os
import time



inp = pd.read_csv('intersection_patient_list.csv', header=None, names=['idx', 'patient'])
tables = list()
for i, p in enumerate(inp.patient):
    try:
        filename = 'results/{}/metrics_ith.csv'.format(p)
        #print(f, filename)
        timestamp = os.path.getmtime(filename)
        df = pd.read_csv(filename, sep='\t')
        df = df.assign(timestamp=timestamp)
        tables.append(df)
    except FileNotFoundError:
        print(i+1, p)
big_table2 = pd.concat(tables, axis=0)
big_table2.to_csv('tmp/ith_metrics_tcga.csv', sep='\t', index=False)

