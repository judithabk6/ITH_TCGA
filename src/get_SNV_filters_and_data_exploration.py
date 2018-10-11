#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to conduct and exploratory analysis of the new tcga data.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import sys
#from util_fonctions import safe_mkdir
import errno
from scipy.stats import pearsonr

import seaborn as sns


pd.set_option('display.max_columns', 200)

plt.style.use('seaborn-white')

def safe_mkdir(path):
    """
    this function takes a path as input, and creates a directory without raising
    if the directory exists (and displaying that it exists)
    """
    try:
        os.mkdir(path)
    # except FileExistsError:  # for python 3.3 and more
    except OSError as e:
        if e.errno == errno.EEXIST:
            print('Directory {} not created. already exists'.format(path))
        else:
            raise

output_path = 'results/mutations_exploratory_analysis'
safe_mkdir('results')
safe_mkdir(output_path)
safe_mkdir('tmp')
folder_list = ['protected_maf_hg38', 'public_somatic_maf_hg38']
#folder_list = ['public_somatic_maf_legacy_hg19', 'public_somatic_maf_hg38']
cancer_loc = sys.argv[1]

filename_dict = dict()
maf_dict = dict()
mut_count_dict = dict()
for folder in folder_list:
    filename_dict[folder] = os.listdir('data/{}/{}'.format(cancer_loc, folder))
    for filename in filename_dict[folder]:
        short_name = filename.split('.')[2]
        maf_dict['{}_{}'.format(folder, short_name)] = pd.read_csv('data/{}/{}/{}'.format(cancer_loc, folder, filename), sep='\t', comment='#', encoding='latin1')
        if 'Start_position' in maf_dict['{}_{}'.format(folder, short_name)].columns:
            maf_dict['{}_{}'.format(folder, short_name)] = maf_dict['{}_{}'.format(folder, short_name)].assign(Start_Position=maf_dict['{}_{}'.format(folder, short_name)].Start_position)
        print(filename, maf_dict['{}_{}'.format(folder, short_name)].shape)
        maf_dict['{}_{}'.format(folder, short_name)] = maf_dict['{}_{}'.format(folder, short_name)].drop_duplicates(subset=['Tumor_Sample_Barcode', 'Start_Position', 'Chromosome', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2'])
        print(filename, maf_dict['{}_{}'.format(folder, short_name)].shape)
        maf_dict['{}_{}'.format(folder, short_name)] = maf_dict['{}_{}'.format(folder, short_name)][maf_dict['{}_{}'.format(folder, short_name)].Variant_Type=='SNP']
        print(filename, maf_dict['{}_{}'.format(folder, short_name)].shape)
        aa = maf_dict['{}_{}'.format(folder, short_name)].groupby('Tumor_Sample_Barcode')
        mut_count_dict['{}_{}'.format(folder, short_name)] = aa.Tumor_Sample_UUID.count()
        maf_dict['{}_{}'.format(folder, short_name)] = maf_dict['{}_{}'.format(folder, short_name)].assign(origin='{}_{}'.format(folder, short_name))
        maf_dict['{}_{}'.format(folder, short_name)] = maf_dict['{}_{}'.format(folder, short_name)].assign(mutation_id=maf_dict['{}_{}'.format(folder, short_name)].apply(lambda x: x['Tumor_Sample_Barcode']+'_'+str(x['Chromosome'])+'_'+str(x['Start_Position'])+'_'+x['Tumor_Seq_Allele1']+'_'+x['Tumor_Seq_Allele2'], axis=1))

# for mutect, possible values of filter are : t_lod_fstar, PASS,
# germline_risk, panel_of_normals, bSeq, alt_allele_in_normal,
# clustered_events, bPcr, homologous_mapping_event,
# multi_event_alt_allele_in_normal, oxog, triallelic_site
# so apart PASS and panel_of_normals (kept in SOMATIC_MAF generation)
# the rest seems to be artifacts        
#maf_dict['protected_maf_hg38_mutect'] = maf_dict['protected_maf_hg38_mutect'][(maf_dict['protected_maf_hg38_mutect'].FILTER=='panel_of_normals') | (maf_dict['protected_maf_hg38_mutect'].FILTER=='PASS')

# for somatic sniper and varscan, already done in the MAF generation from
# some things have been flagged, but it is kind of stringent. But it seems
# reasonable and better than designing custom filters

MAF_cols = [u for u in maf_dict['protected_maf_hg38_mutect'].columns if 'MAF' in u]
protected_maf_hg38_maf_names = [i for i in maf_dict.keys() if "protected_maf_hg38" in i]
for n in protected_maf_hg38_maf_names:
    maf_dict[n] = maf_dict[n][(maf_dict[n].FILTER == 'panel_of_normals') |
                              (maf_dict[n].FILTER == 'PASS') |
                              (maf_dict[n].FILTER == 'Tier1') |
                              (maf_dict[n].FILTER == 'Tier2') |
                              (maf_dict[n].FILTER == 'Tier3') |
                              (maf_dict[n].FILTER == 'Tier4') |
                              (maf_dict[n].FILTER == 'Tier5')]
    maf_dict[n] = maf_dict[n].assign(
        global_maf=maf_dict[n][MAF_cols].max(axis=1))
    mut_count_dict[n] = maf_dict[n].groupby('Tumor_Sample_Barcode').Tumor_Sample_UUID.count()

very_big_one = pd.concat([maf_dict[i] for i in protected_maf_hg38_maf_names], axis=0)

aa = very_big_one.groupby('mutation_id')#.origin.count()
useful_final = pd.DataFrame(dict({'nb_caller': aa.origin.count(),
                                  'db_snp_status': aa.Existing_variation.first(),
                                  'exac_freq': aa.ExAC_AF.mean(),
                                  '1000G_freq': aa.GMAF.mean(),
                                  'COSMIC': aa.COSMIC.first(),
                                  'PHENO': aa.PHENO.first(),
                                  'IMPACT': aa.IMPACT.first(),
                                  'SIFT': aa.SIFT.first(),
                                  'PolyPhen': aa.PolyPhen.first(),
                                  'Hugo_Symbol': aa.Hugo_Symbol.first(),
                                  'Gene': aa.Gene.first(),
                                  'SYMBOL': aa.SYMBOL.first(),
                                  'SOMATIC': aa.SOMATIC.first(),
                                  'GDC_FILTER': aa.GDC_FILTER.first(),
                                  'CONTEXT': aa.CONTEXT.first()},
                            **{i: aa[i].mean() for i in ['t_depth', 't_ref_count', 't_alt_count', 'n_depth', 'n_ref_count', 'n_alt_count']}))
useful_final = useful_final.assign(t_vaf=useful_final.t_alt_count/useful_final.t_depth)
useful_final = useful_final.assign(n_vaf=useful_final.n_alt_count/useful_final.n_depth)
useful_final = useful_final.assign(mutation_id=useful_final.index)
useful_final = useful_final.assign(sample_id=useful_final.mutation_id.str.split('_').str[0])
useful_final_qc = useful_final[(((~(useful_final['1000G_freq'] > 0.01)) &
                                 (~(useful_final['exac_freq'] > 0.01))) | (~pd.isnull(useful_final.COSMIC)) | (useful_final.SOMATIC==1) | (useful_final.PHENO==1))&
                                 (useful_final.n_depth>=6) &
                                 ((useful_final.n_alt_count<=1) | (useful_final.n_vaf<=0.01)) & 
                                 (useful_final.t_depth>=8) &
                                 ((useful_final.t_alt_count>=3) | (useful_final.t_vaf>0.2))]

useful_final_2 = useful_final[useful_final.nb_caller>=2]
useful_final_2_qc = useful_final[(useful_final.nb_caller>=2) & (((~(useful_final['1000G_freq'] > 0.01)) &
                                 (~(useful_final['exac_freq'] > 0.01))) | (~pd.isnull(useful_final.COSMIC)) | (useful_final.SOMATIC==1) | (useful_final.PHENO==1)) &
                                 (useful_final.n_depth>=6) &
                                 ((useful_final.n_alt_count<=1) | (useful_final.n_vaf<=0.01)) & 
                                 (useful_final.t_depth>=8) &
                                 ((useful_final.t_alt_count>=3) | (useful_final.t_vaf>0.2))]
useful_final_4 = useful_final[useful_final.nb_caller>=4]


mut_count_dict['protected_union'] = useful_final.groupby('sample_id').mutation_id.count()
mut_count_dict['protected_union_qc'] = useful_final_qc.groupby('sample_id').mutation_id.count()
mut_count_dict['protected_union_2'] = useful_final_2.groupby('sample_id').mutation_id.count()
mut_count_dict['protected_union_2_qc'] = useful_final_2_qc.groupby('sample_id').mutation_id.count()
mut_count_dict['protected_union_4'] = useful_final_4.groupby('sample_id').mutation_id.count()

public_somatic_maf_hg38_maf_names = [i for i in maf_dict.keys()
                                     if "public_somatic_maf_hg38" in i]
very_big_one_public = pd.concat([maf_dict[i] for i in public_somatic_maf_hg38_maf_names], axis=0)
aa_public = very_big_one_public.groupby('mutation_id')
cols_to_consider = ['t_depth', 't_ref_count', 't_alt_count', 'n_depth',
                    'n_ref_count', 'n_alt_count']
useful_final_public = pd.DataFrame(dict({'nb_caller': aa_public.origin.count(),
                                         'db_snp_status': aa_public.Existing_variation.first(),
                                         'exac_freq': aa_public.ExAC_AF.mean(),
                                         '1000G_freq': aa_public.GMAF.mean(),
                                         'COSMIC': aa_public.COSMIC.first(),
                                         'PHENO': aa_public.PHENO.first(),
                                         'IMPACT': aa_public.IMPACT.first(),
                                         'SIFT': aa_public.SIFT.first(),
                                         'PolyPhen': aa_public.PolyPhen.first(),
                                         'Hugo_Symbol': aa_public.Hugo_Symbol.first(),
                                         'Gene': aa_public.Gene.first(),
                                         'SYMBOL': aa_public.SYMBOL.first(),
                                         'SOMATIC': aa_public.SOMATIC.first(),
                                         'GDC_FILTER': aa_public.GDC_FILTER.first(),
                                         'CONTEXT': aa_public.CONTEXT.first()},
                                     **{i: aa_public[i].mean() for i in cols_to_consider}))
useful_final_public = useful_final_public.assign(mutation_id=useful_final_public.index)
useful_final_public = useful_final_public.assign(sample_id=useful_final_public.mutation_id.str.split('_').str[0])
mut_count_dict['public_union'] = useful_final_public.groupby('sample_id').mutation_id.count()

zz = pd.DataFrame(data=mut_count_dict)
zz_reasonable = zz[zz.apply(lambda x: np.abs(x - x.mean()) / x.std() < 3).all(axis=1)]
zz_nona = zz_reasonable.dropna()

def corrfunc(x, y, **kws):
    (r, p) = pearsonr(x, y)
    ax = plt.gca()
    ax.annotate("r = {:.2f}".format(r),
                xy=(.1, .9), xycoords=ax.transAxes)
    ax.annotate("p = {:.3f}".format(p),
                xy=(.4, .9), xycoords=ax.transAxes)
sns.set_context("notebook", font_scale=2.5, rc={"lines.linewidth": 2.5})
g = sns.PairGrid(zz_nona, diag_sharey=False, size=8)
g.map_lower(sns.regplot)
g.map_lower(corrfunc)
g.map_upper(sns.kdeplot, cmap="Blues_d")
g.map_diag(sns.kdeplot, lw=3)
plt.xticks(rotation=70)

plt.savefig('{}/{}_nb_mutations_all_maf_correlation.pdf'.format(output_path, cancer_loc), bbox_inches='tight')
plt.close()

zz_nona.sum(axis=0).plot.bar(fontsize=20)
plt.savefig('{}/{}_nb_mutations_total.pdf'.format(output_path, cancer_loc), bbox_inches='tight')
plt.close()

a = zz_nona.stack()
a.index = a.index.droplevel(0)
a = pd.DataFrame(a, columns=['nb_mut'])
a = a.assign(caller=a.index)
a = a.assign(log_nb_mut=np.log10(a.nb_mut))

sns.set(context='notebook', font_scale=5, rc={"axes.facecolor": (0, 0, 0, 0)})
pal = sns.cubehelix_palette(15, rot=-.25, light=.7)
g = sns.FacetGrid(a, row="caller", hue="caller", aspect=15, size=3, palette=pal)

# Draw the densities in a few steps
g.map(sns.kdeplot, "log_nb_mut", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
g.map(sns.kdeplot, "log_nb_mut", clip_on=False, color="w", lw=2, bw=.2)

g.map(plt.axhline, y=0, lw=2, clip_on=False)

def label(x, color, label):
    ax = plt.gca()
    ax.text(1, .2, label, fontweight="bold", color=color,
            ha="left", va="center", transform=ax.transAxes)

g.map(label, "log_nb_mut")

# Set the subplots to overlap
g.fig.subplots_adjust(hspace=-.2)

# Remove axes details that don't play will with overlap
g.set_titles("")
g.set(yticks=[])
ax = plt.gca()
labels = [item.get_text() for item in ax.get_xticklabels()]
new_labels=[str(10**int(i.replace('−', '-'))) for i in labels]

ax.set_xticklabels(new_labels)
g.despine(bottom=True, left=True)
plt.savefig('{}/{}_nb_mutations_densities.pdf'.format(output_path, cancer_loc), bbox_inches='tight')
plt.close()

useful_final_public.to_csv('tmp/useful_final_public_{}.csv'.format(cancer_loc), sep='\t')
useful_final_qc.to_csv('tmp/useful_final_qc_{}.csv'.format(cancer_loc), sep='\t')

