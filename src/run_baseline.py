#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to perform an baseline very simple
ITH detection.
"""

import matplotlib
matplotlib.use('Agg')
import sys
import pandas as pd
import numpy as np
import peakutils
import matplotlib.pyplot as plt
import seaborn as sns
from util_functions import safe_mkdir

pd.set_option('display.max_columns', 200)

plt.style.use('seaborn-white')

output_path = 'results/20170118_exploratory_analysis'
safe_mkdir(output_path)

patient = sys.argv[1]
folder = sys.argv[2]
input_filename = sys.argv[3]

snv_merge = pd.read_csv(input_filename, sep='\t')


def detect_elbows(y, indexes, derivative_thresh=0.02, height_thresh=0.15, min_dist=5):
    abs_height_thresh = height_thresh * max(abs(y))
    dy = np.hstack((y[1:] - y[:-1], [max(y)]))
    potential_idx = np.hstack((peakutils.indexes(dy, thres=0), peakutils.indexes(-dy, thres=0)))
    peaks = list()
    for p in potential_idx:
        if (np.abs(dy[p]) < derivative_thresh) & (np.abs(y[p]) > abs_height_thresh):
            peaks.append(p)
    peaks = np.array(peaks).astype(int)
    peaks = np.hstack((peaks, indexes)).astype(int)
    if peaks.size > 1 and min_dist > 1:
        highest = peaks[np.argsort(dy[peaks])][::-1]
        rem = np.ones(dy.size, dtype=bool)
        rem[peaks] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peak] = False

        peaks = np.arange(dy.size)[~rem]
    return peaks


output_path_patients = '{}/20180801_patient_detail_analysis_{}'.format(output_path, folder)
safe_mkdir(output_path_patients)
# np.random.seed(12)
# list_of_patients = uu.PATIENT_ID.unique().tolist()
with open('{}/{}_baseline_naive_ith_method.csv'.format(output_path_patients, patient), 'w') as outfile:
    try:
        rr = list()
        sub_snv_merge = snv_merge[snv_merge.patient_id == patient]
        fig, axes = plt.subplots(1, 5, figsize=(20, 5))  # change figsize, get more examples, add peak detection
        for i, c in enumerate((sub_snv_merge.vaf, sub_snv_merge.vaf_cn, sub_snv_merge.vaf_purity, sub_snv_merge.vaf_CFE)):
            f = sns.distplot(c, rug=True, ax=axes[i])
            x, y = f.lines[0].get_data()
            indexes = peakutils.indexes(y, thres=0.2, min_dist=0.1)
            axes[i].plot(x[indexes], y[indexes], 'bo')
            rr.append(str(len(indexes)))
            rr.append(str(x[indexes[-1]]))
            idx = detect_elbows(y, indexes)
            axes[i].plot(x[idx], y[idx], 'c*', mec='c', ms=6)
            rr.append(str(len(idx)))
            rr.append(str(x[idx[-1]]))
            if i == 2:
                output_df = pd.DataFrame(columns=['cluster_id', 'cluster_size', 'mean'])
                cluster_assignment = np.argmin(np.sqrt((np.repeat(x[idx], len(sub_snv_merge)).reshape(len(idx), len(sub_snv_merge)) - np.repeat(np.array([sub_snv_merge.vaf_purity]), len(idx), axis=0))**2), axis=0)
                for j in range(len(idx)):
                    output_df = output_df.append({'cluster_id': j, 'cluster_size': sum(cluster_assignment == j), 'mean': x[idx[j]]}, ignore_index=True)
                output_path_results = 'results/{}/baseline/{}'.format(patient, folder)
                safe_mkdir('results/{}/baseline'.format(patient))
                safe_mkdir(output_path_results)
                output_df = output_df.assign(cluster_id=output_df.cluster_id.astype(int))
                output_df = output_df.assign(cluster_size=output_df['cluster_size'].astype(int))
                output_df.to_csv('{}/cluster_assignment_bis.csv'.format(output_path_results), sep='\t', index=False)
                sub_snv_merge = sub_snv_merge.assign(cluster=cluster_assignment)
                sub_snv_merge = sub_snv_merge.assign(mutation_id=sub_snv_merge['mutation_id.1'].str[29:-2] + sub_snv_merge['mutation_id.1'].str[-1])
                sub_snv_merge[['mutation_id', 'vaf_cn', 'cluster']].to_csv('{}/cluster_assignment_for_each_mutation.csv'.format(output_path_results), sep='\t', index=False)

        dy = np.hstack((y[1:] - y[:-1], [0]))
        ddy = np.hstack(([0], dy[1:] - dy[:-1]))
        axes[4].plot(x, y, 'r')
        axes[4].plot(x, dy, 'b')
        axes[4].plot(x, ddy, 'g')
        axes[4].plot([x[0], x[-1]], [0, 0], 'k')
        axes[3].plot([], [], 'bo', label='peakutils on f')
        axes[3].plot([], [], 'c*', mec='c', ms=6, label='peakutils on f, df -df')
        axes[3].legend(loc='best')
        plt.savefig('{}/{}_vaf_hist.pdf'.format(output_path_patients, patient), bbox_inches='tight')
        plt.clf()

    except ValueError:
        print('nan for patient {}'.format(patient))
        rr = [str(np.nan)]*16
    outfile.write('{}\t{}\n'.format(patient, '\t'.join(rr)))
