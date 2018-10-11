#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to convert all methods outputs to the input of CSR
https://github.com/kaixiany/CSR
that is tab separated file with 4 cols: chr pos cluster cp
"""

import sys
import pandas as pd


patient = sys.argv[1]
#folder = sys.argv[2]

"""
# for testing purpose
patient = "TCGA-3C-AALK"
folder = "protected_hg38_vcf"


former baseline overwritten
In [48]: m.head() (old result)
Out[48]: 
   cluster_id  cluster_size      mean
0           0           205  0.437047

In [49]: m = pd.read_csv('results/{}/{}/{}/cluster_assignment.csv'.format(patient, ith_method, folder), sep='\t')

In [50]: m (new result)
Out[50]: 
   cluster_id  cluster_size      mean
0           0            98  0.429730
1           1           107  0.642925
"""
out_cols = ['chr', 'pos', 'cluster', 'cp']


def write_csr_input(ith_method, folder, patient):
    if ith_method=='sciclone':
        try:
            sci = pd.read_csv('results/{}/sciclone/{}/clusters1'.format(patient, folder), sep='\t')
            sci = sci.assign(pos=sci.st)
            sci = sci.assign(cp=sci['tumor.vaf']/100)
            sci[out_cols].to_csv('results/{}/sciclone/{}/input_csr.tsv'.format(patient, folder), sep='\t', index=False)
        except IOError:
            print('no output available for sciclone patient {} and folder {}'.format(patient, folder))

    elif ith_method=='pyclone':
        try:
            pyloci = pd.read_csv('results/{}/pyclone/{}/tables/loci.tsv'.format(patient, folder), sep='\t')
            pyloci = pyloci.assign(chr=pyloci.mutation_id.str.split('_').str[0].str.replace('chr', ''))
            pyloci = pyloci.assign(pos=pyloci.mutation_id.str.split('_').str[1])
            pyloci = pyloci.assign(cluster=pyloci.cluster_id)
            pyloci = pyloci.assign(cp=pyloci.cellular_prevalence)
            pyloci[out_cols].to_csv('results/{}/pyclone/{}/input_csr.tsv'.format(patient, folder), sep='\t', index=False)
        except IOError:
            print('no output available for pyclone patient {} and folder {}'.format(patient, folder))

    elif ith_method=='PhyloWGS':
        try:
            m = pd.read_csv('results/{}/{}/{}/outputs/1C.txt'.format(patient, ith_method, folder), sep='\t', names=['cluster_id', 'size', 'mean'], header=None)
            input_mut_phylo = pd.read_csv('results/{}/{}/{}/parsed_ssm_data.txt'.format(patient, ith_method, folder), sep='\t')
            c_assignment = pd.read_csv('results/{}/{}/{}/outputs/2A.txt'.format(patient, ith_method, folder), sep='\t', names=['cluster'], header=None)
            input_mut_phylo = input_mut_phylo.assign(cluster=c_assignment.cluster)
            input_mut_phylo = pd.merge(input_mut_phylo, m, left_on='cluster', right_on='cluster_id', how='left')
            input_mut_phylo = input_mut_phylo.assign(chr=input_mut_phylo.gene.str.split('_').str[0])
            input_mut_phylo = input_mut_phylo.assign(pos=input_mut_phylo.gene.str.split('_').str[1])
            input_mut_phylo = input_mut_phylo.assign(cp=input_mut_phylo['mean'])
            input_mut_phylo[out_cols].to_csv('results/{}/{}/{}/input_csr.tsv'.format(patient, ith_method, folder), sep='\t', index=False)

        except IOError:
            print('no output available for PhyloWGS patient {} and folder {}'.format(patient, folder))

    elif ith_method == 'baseline':
        try:
            m = pd.read_csv('results/{}/{}/{}/cluster_assignment_for_each_mutation.csv'.format(patient, ith_method, folder), sep='\t')
            m = m.assign(chr=m.mutation_id.str.split('_').str[0].str.replace('chr', ''))
            m = m.assign(pos=m.mutation_id.str.split('_').str[1])
            m = m.assign(cp=m.vaf_cn)
            m[out_cols].to_csv('results/{}/{}/{}/input_csr.tsv'.format(patient, ith_method, folder), sep='\t', index=False)
        except IOError:
            print('no output available for baseline patient {} and folder {}'.format(patient, folder))

    elif ith_method == 'expands':
        try:
            all_mut = pd.read_csv('results/{}/{}/{}/{}__{}.sps'.format(patient, ith_method, folder, patient, folder), sep='\t', skiprows=1)
            all_mut = all_mut.assign(pos=all_mut.startpos)
            all_mut = all_mut.assign(cp=all_mut.AF_Tumor)
            m = all_mut.SP.value_counts().to_frame()
            m = m.assign(cluster=range(len(m)))
            all_mut = pd.merge(all_mut, m[['cluster']], left_on='SP', right_index=True)
            all_mut[out_cols].to_csv('results/{}/{}/{}/input_csr.tsv'.format(patient, ith_method, folder), sep='\t', index=False)
        except IOError:
            print('no output available for expands patient {} and folder {}'.format(patient, folder))
for ith_method in ['pyclone', 'PhyloWGS', 'sciclone', 'baseline', 'expands']:
    for folder in ['protected_hg38_vcf', 'public_hg38_vcf']:
        write_csr_input(ith_method, folder, patient)
