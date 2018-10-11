#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to make plots from results for the main article.
"""
from scipy import stats
import seaborn as sns
import matplotlib as mpl
from matplotlib.text import Text
from matplotlib.patches import Patch
from scipy.stats import combine_pvalues
import matplotlib.colors as colors
from matplotlib.legend import Legend
import venn
import matplotlib.pyplot as plt
import pandas as pd
from statsmodels.stats.multitest import multipletests
import numpy as np
# get variables without duplicating code.
from perform_survival_regression import restricted_df_dict, tmp_df_c, nb_clones_col, cols_restricted, tmp_df

plt.style.use('seaborn-white')

cancer_locs = ['BRCA', 'HNSC', 'BLCA']

colors_protected = {'Baseline': '#A1E4C5', 'SciClone': '#AE88FF',
                    'MATH': '#C9AB3C', 'PhyloWGS': '#FF7DA4', 'CSR': '#9B9B9B',
                    'PyClone': '#84B8FF', 'Expands': '#FFC89B',
                    'Combination': '#FFFFFF'}
output_path = 'results/20180801_results_rerun_newsksurv'

class Handler(object):
    def __init__(self, color1, color2):
        self.color1 = color1
        self.color2 = color2

    def legend_artist(self, legend, orig_handle, fontsize, handlebox):
        x0, y0 = handlebox.xdescent, handlebox.ydescent
        width, height = handlebox.width, handlebox.height
        patch = plt.Rectangle([x0, y0], width, height, facecolor=self.color1,
                              edgecolor='k', transform=handlebox.get_transform())
        patch2 = plt.Rectangle([x0+width/2., y0], width/2., height, facecolor=self.color2,
                               edgecolor='k', transform=handlebox.get_transform())
        handlebox.add_artist(patch)
        handlebox.add_artist(patch2)
        return patch

mpl.rcParams['hatch.linewidth'] = 2  # previous pdf hatch linewidth

# new plots for the paper
for k in restricted_df_dict.keys():

    # build regression_df
    loc_df_list = list()
    for loc in ['BRCA', 'BLCA', 'HNSC']:
        for fold in range(1, 6):
            filename = 'tmp/20180801_{}_test_results/{}_fold{}_test_res_cindex_r.csv'.format(k, loc, fold)
            sub_reg_df = pd.read_csv(filename, sep='\t')
            sub_reg_df = sub_reg_df.assign(cancer_loc=loc)
            loc_df_list.append(sub_reg_df)
    regression_df = pd.concat(loc_df_list, axis=0)
    regression_df = regression_df.assign(ith_method=pd.Categorical(regression_df.ith_method.str.replace('pyclone', 'PyClone')
                                                                                           .str.replace('sciclone', 'SciClone')
                                                                                           .str.replace('math_score', 'MATH')
                                                                                           .str.replace('combination', 'Combination')
                                                                                           .str.replace('baseline', 'Baseline')
                                                                                           .str.replace('expands', 'Expands'),
                                                                   categories=['Baseline', 'Expands', 'MATH', 'CSR', 'PhyloWGS', 'PyClone', 'SciClone', 'Combination'], ordered=True))
    regression_df = regression_df.assign(folder=regression_df.folder.str.split('_').str[0])
    regression_df = regression_df.assign(new_method_folder=regression_df.ith_method.astype(str)+'+'+regression_df.folder)
    regression_df = regression_df.assign(conditions=regression_df.cancer_loc + ' & ' + regression_df.new_method_folder+' & '+regression_df.additional_cols)
    gr = regression_df.groupby('conditions').agg({'test_score_cindex_r':[np.mean, np.std, np.median], 'pval_test_cindex_r':lambda x: combine_pvalues(x)[1]})
    gr = gr.assign(corrected_pvalue=multipletests(gr['pval_test_cindex_r'].values.ravel(), method='fdr_bh')[1])
    gr = gr.assign(cancer_loc=gr.index.str.split(' & ').str[0])
    gr = gr.assign(new_method_folder=gr.index.str.split(' & ').str[1])
    gr = gr.assign(additional_cols=gr.index.str.split(' & ').str[2])

pval_cols = ['pvalue_fold{}'.format(i) for i in range(1, 6)]
loc_df_list = list()
for k in restricted_df_dict.keys():
    for loc in ['BRCA', 'BLCA', 'HNSC']:
            filename = 'tmp/20180801_{}_test_results/{}_comparison_cindex_r.csv'.format(k, loc)
            sub_reg_df = pd.read_csv(filename, sep='\t')
            sub_reg_df = sub_reg_df.assign(cancer_loc=loc)
            sub_reg_df = sub_reg_df.assign(setting=k)
            sub_reg_df = sub_reg_df.assign(combined_pval=sub_reg_df.apply(lambda x: combine_pvalues(np.array(x[pval_cols].values).astype(float))[1], axis=1))
            loc_df_list.append(sub_reg_df)

comp_df = pd.concat(loc_df_list, axis=0)
comp_df_clinical = comp_df[(comp_df.setting1 == 'clonality+clinical') & (comp_df.setting2 == 'clinical_alone')]
comp_df_clinical = comp_df_clinical.assign(corrected_pvalue=multipletests(comp_df_clinical['combined_pval'].values.ravel(), method='fdr_bh')[1])
# nothing is significant for clinical after multiple testing correction
comp_df = comp_df.assign(corrected_pvalue=multipletests(comp_df['combined_pval'].values.ravel(), method='fdr_bh')[1])
comp_df_metrics = comp_df[(comp_df.setting1.isin(['clonality', 'nb_clones', 'clonality_squared'])) & (comp_df.setting2.isin(['nb_clones', 'clonality']))]
comp_df_metrics = comp_df_metrics.assign(corrected_pvalue=multipletests(comp_df_metrics['combined_pval'].values.ravel(), method='fdr_bh')[1])


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


cmap = plt.get_cmap('PiYG')
new_cmap = truncate_colormap(cmap, 0.2, 0.8)

#########################
# merged plot for paper #
#########################

# venn diagrams for the number of patients
protected_nb_clones_cols = relevant_nb_clones = [c for c in nb_clones_col if 'protected' in c]
for c in protected_nb_clones_cols:
    public_col = c.replace('protected', 'public')
    merged_col = c.replace('protected', 'merged')
    tmp_df = tmp_df.assign(**{merged_col: tmp_df.apply(lambda x: np.nan if pd.isnull(x[c]) or pd.isnull(x[public_col]) else 1, axis=1)})

fig, axes = plt.subplots(3, 4, figsize=(35, 20))
for ii, folder in enumerate(['public', 'protected', 'merged']):
    relevant_nb_clones = [c for c in tmp_df.columns if (folder in c) and ('nb_clones' in c) and ('CSR' not in c) and ('baseline' not in c) and ('norm' not in c)]
    for i, loc in enumerate(cancer_locs):
        relevant_tmp_df = tmp_df[tmp_df.cancer_loc==loc]
        labels = venn.get_labels([relevant_tmp_df[c].dropna().index for c in relevant_nb_clones], fill=['number', 'logic'])
        labels_short = {k: v.split(': ')[1] for k,v in labels.items()}
        label_names = [c.split('_')[0].replace('pyclone', 'PyClone').replace('sciclone', 'SciClone').replace('baseline', 'Baseline').replace('expands', 'Expands') for c in relevant_nb_clones]
        label_names_fig = ['' for l in label_names]
        ax = venn.venn4(labels_short, names=label_names_fig, colors=[colors_protected[i] for i in label_names], fontsize_text=22, fontsize_number=25, ax=axes[ii, i])
        nb_samples_total = relevant_tmp_df[relevant_nb_clones].dropna(axis='index', how='all').shape[0]
        ax.set_title('{loc} (n={n})'.format(loc=loc, n=nb_samples_total), fontsize=40)
        if i == 0:
            ax.text(-0.1, 0.4, '{}'.format(folder).replace('merged', 'public and protected\nintersection'),
                    fontsize=40, rotation=90, va='center', ha='center')


    labels = venn.get_labels([tmp_df[c].dropna().index for c in relevant_nb_clones], fill=['number', 'logic'])
    labels_short = {k: v.split(': ')[1] for k, v in labels.items()}
    label_names = [c.split('_')[0].replace('pyclone', 'PyClone').replace('sciclone', 'SciClone').replace('baseline', 'Baseline').replace('expands', 'Expands') for c in relevant_nb_clones]
    label_names_fig = ['' for l in label_names]
    ax = venn.venn4(labels_short, names=label_names_fig, colors=[colors_protected[i] for i in label_names], fontsize_text=22, fontsize_number=25, ax=axes[ii, 3])
    if ii==0:
        ax.legend(label_names, loc='center left', bbox_to_anchor=(1.1, 0.5), fontsize=30)
    nb_samples_total = tmp_df[relevant_nb_clones].dropna(axis='index', how='all').shape[0]
    ax.set_title('All 3 types (n={n})'.format(n=nb_samples_total), fontsize=40)
plt.savefig('{}/total_venn_total.pdf'.format(output_path), bbox_inches='tight')
plt.close()

# results of clonality alone
plt.style.use('seaborn-white')
fig, axes = plt.subplots(1, 3, figsize=(23, 5))
for ii, loc in enumerate(['BRCA', 'BLCA', 'HNSC']):
    relevant_regression_df = regression_df[(regression_df.cancer_loc == loc)]

    color_dict = {'{}+public'.format(k): v for k, v in colors_protected.items()}
    color_dict.update({'{}+protected'.format(k): v for k, v in colors_protected.items()})
    # ax=sns.boxplot(x='new_method_folder', y='test_score', data=relevant_regression_df[relevant_regression_df.additional_cols=='clonality'], palette=color_dict, order=sorted(relevant_regression_df.new_method_folder.unique().tolist())[:-1], whis='range')
    sns.boxplot(x='new_method_folder', y='test_score_cindex_r',
                data=relevant_regression_df[relevant_regression_df.additional_cols == 'clonality'],
                palette=color_dict, order=sorted(relevant_regression_df.new_method_folder.unique().tolist())[:-1],
                whis='range', ax=axes[ii])
    n = len(axes[ii].artists)
    for i in range(1, n, 2):
        g = axes[ii].artists[i]
        g.set_hatch('//')
    new = list()
    axes[ii].set_xticks(new)
    axes[ii].set_xlabel('')
    axes[ii].set_ylabel('test C index')
    axes[ii].yaxis.grid(which="major", color='k', linestyle='--', linewidth=1)
    for i, method_folder in enumerate(sorted(relevant_regression_df.new_method_folder.unique().tolist())[:-1]):
        try:
            p = gr[(gr.new_method_folder==method_folder)&(gr.additional_cols=='clonality')&(gr.cancer_loc==loc)]['corrected_pvalue'].values[0]
            if p >= 0.05:
                displaystring = r'n.s.'
            elif p < 0.001:
                displaystring = r'***'
            elif p < 0.01:
                displaystring = r'**'
            elif p < 0.05:
                displaystring = r'*'
            if p < 0.05:
                print(p, method_folder, loc, k)
                axes[ii].text(i, 0.7, displaystring, ha='center', va='bottom',
                              color='r', fontsize=20)
        except:
            pass
    axes[ii].plot(axes[ii].get_xlim(), [0.5, 0.5], 'k', lw=2)
    axes[ii].set_ylim([0.2, 0.8])
    axes[ii].set_title(loc, fontsize=30)

legend_items = list()
for key in sorted(colors_protected.keys()):
    if relevant_regression_df.ith_method.value_counts()[key]:
        legend_items.append(Patch(facecolor=colors_protected[key], edgecolor='black', label=key))

# build custom legend
other_handles = [Patch(facecolor='white', edgecolor='black', label='protected'),
                 Patch(facecolor='white', edgecolor='black', label='public', hatch='//')]
other_ = axes[-1].add_artist(plt.legend(handles=other_handles,
                             bbox_to_anchor=(1.05, 0.2), fontsize=15,
                             title='input mutations', loc=2))

first_legend = axes[-1].legend(handles=legend_items, bbox_to_anchor=(1.05, 1),
                               loc=2, fontsize=15, title='ITH methods')


plt.savefig('{}/merged_clonality.pdf'.format(output_path), bbox_inches='tight',
            bbox_extra_artists=(first_legend,other_))
plt.close()


# results of clonality + clinical variables
fig, axes = plt.subplots(1, 3, figsize=(23, 5))
for ii, loc in enumerate(['BRCA', 'BLCA', 'HNSC']):
    relevant_regression_df = regression_df[(regression_df.cancer_loc == loc)]

    color_dict = {'{}+public'.format(k):v for k, v in colors_protected.items()}
    color_dict.update({'{}+protected'.format(k):v for k, v in colors_protected.items()})
    color_dict['no_ITH+public'] = '#FFFFFF'
    color_dict['no_ITH+protected'] = '#FFFFFF'
    color_dict['no_ITH'] = '#FFFFFF'
    relevant_regression_df.loc[relevant_regression_df.new_method_folder == 'nan+none', 'new_method_folder'] = 'no_ITH'
    sns.boxplot(x='new_method_folder', y='test_score_cindex_r',
                data=relevant_regression_df[relevant_regression_df.additional_cols=='clonality+clinical'],
                palette=color_dict, order=['no_ITH']+sorted(relevant_regression_df.new_method_folder.unique().tolist())[:-1],
                whis='range', ax=axes[ii])
    for i in range(1):
        g = axes[ii].artists[i]
        g.set_edgecolor('red')
        for j in range(i*6, i*6+6):
            line = axes[ii].lines[j]
            line.set_color('red')
            line.set_mfc('red')
            line.set_mec('red')
    n = len(axes[ii].artists)
    for i in range(2, n, 2):
        g = axes[ii].artists[i]
        g.set_hatch('//')
    new = list()
    axes[ii].set_xticks(new)
    axes[ii].set_xlabel('')
    axes[ii].yaxis.grid(which="major", color='k', linestyle='--', linewidth=1)
    axes[ii].set_ylabel('test C index')

    axes[ii].plot(axes[ii].get_xlim(), [0.5, 0.5], 'k', lw=2)
    axes[ii].set_ylim([0.1, 0.9])
    axes[ii].set_title(loc, fontsize=30)

legend_items = list()
for key in sorted(colors_protected.keys()):
    if relevant_regression_df.ith_method.value_counts()[key]:
        legend_items.append(Patch(facecolor=colors_protected[key], edgecolor='black', label=key))

legend_items.append(Patch(facecolor='white', edgecolor='red', lw=2, label='no_ITH'))
other_handles = [Patch(facecolor='white', edgecolor='black', label='protected'),
                 Patch(facecolor='white', edgecolor='black', label='public', hatch='//')]
other_ = axes[-1].add_artist(plt.legend(handles=other_handles,
                                        bbox_to_anchor=(1.05, 0.2),
                                        fontsize=15, title='input mutations',
                                        loc=2))

first_legend = axes[-1].legend(handles=legend_items, bbox_to_anchor=(1.05, 1),
                               loc=2, fontsize=15, title='ITH methods')


plt.savefig('{}/merged_clonality_clinical.pdf'.format(output_path),
            bbox_inches='tight', bbox_extra_artists=(first_legend, other_))
plt.close()


# other metrics relevance
fig, axes = plt.subplots(3, 1, figsize=(9, 18))
for ii, loc in enumerate(cancer_locs):
    relevant_regression_df = regression_df[(regression_df.cancer_loc == loc)]
    extract = relevant_regression_df[relevant_regression_df.additional_cols.isin(['clonality', 'nb_clones', 'clonality_squared', 'nb_clones_squared'])]

    sns.boxplot(x='new_method_folder', y='test_score_cindex_r',
                  data=extract,
                  palette='Greys',
                  order=sorted(relevant_regression_df.new_method_folder.unique().tolist())[:-1],
                  hue='additional_cols', zorder=3, ax=axes[ii], whis='range')
    n = len(axes[ii].get_xticks())
    t = axes[ii].get_xticks()

    for i in range(1, n, 2):
        axes[ii].axvspan(i-1.5, i-0.5, color=colors_protected[axes[ii].get_xticklabels()[i].get_text().split('+')[0]], zorder=-5, alpha=0.7)
        axes[ii].axvspan(i-0.5, i+0.5, color=colors_protected[axes[ii].get_xticklabels()[i].get_text().split('+')[0]], zorder=-4, alpha=0.7)
        axes[ii].axvspan(i-0.5, i+0.5, color='none', zorder=-3, alpha=1, hatch='//', edgecolor='gray')
    axes[ii].set_xticks([])
    axes[ii].set_xlabel('')
    axes[ii].plot(axes[ii].get_xlim(), [0.5, 0.5], 'k', lw=2)
    axes[ii].set_ylim([0.1, 0.9])
    axes[ii].set_ylabel('test C index')
    plt.setp(axes[ii].get_xticklabels(), rotation=40, ha='right', fontsize=10)
    axes[ii].yaxis.grid(which="major", color='k', linestyle='--', linewidth=1)
    axes[ii].legend_.remove()
    axes[ii].set_title(loc, fontsize=30)
    if ii == 2:
        ith_legend = axes[ii].add_artist(plt.legend(handles=legend_items, bbox_to_anchor=(1.05, 1), loc=2, fontsize=12, title='ITH methods'))
        other_ = axes[ii].add_artist(plt.legend(handles=other_handles, bbox_to_anchor=(1.05, 0.4), fontsize=12, title='input mutations', loc=2))
        clinical_legend = axes[ii].legend(bbox_to_anchor=(1.05, 0.2), loc=2, borderaxespad=0., fontsize=12, title='model features')

plt.savefig('{}/merged_clonality_other_metrics.pdf'.format(output_path, loc, k),
            bbox_inches='tight',
            bbox_extra_artists=(ith_legend,other_, clinical_legend), aspect=0.5)
plt.close()


# ith_composition
fig, axes = plt.subplots(2, 2, figsize=(15, 12))
bkp_nb_clones_col = ['{}_bkp'.format(c) for c in nb_clones_col]
for ii, loc in enumerate(cancer_locs):
    ax = axes[ii//2, ii % 2]
    df = tmp_df_c.copy()
    sub_df = df[df.cancer_loc == loc].dropna(subset=nb_clones_col)
    for ith_method in ['pyclone', 'CSR', 'PhyloWGS', 'sciclone', 'baseline', 'expands']:
        for folder in ['protected_hg38_vcf', 'public_hg38_vcf']:
            sub_df['{}_{}_nb_clones_bkp'.format(ith_method, folder)] = sub_df['{}_{}_nb_clones'.format(ith_method, folder)].copy()
            sub_df.loc[sub_df['{}_{}_nb_clones'.format(ith_method, folder)]>=6, '{}_{}_nb_clones'.format(ith_method, folder)] = '6+'

    aa = pd.melt(sub_df, value_vars=nb_clones_col).\
        pivot_table(index='value', columns=['variable'], aggfunc=len).fillna(0).T
    aa.columns = [int(c) if type(c) != str else c for c in aa.columns]
    aa.plot.bar(stacked=True, colormap=new_cmap, edgecolor='gray', ax=ax)
    ee = pd.melt(sub_df, value_vars=bkp_nb_clones_col).\
        pivot_table(index='value', columns=['variable'],
                    aggfunc=len).fillna(0).apply(lambda x:  x/sum(x), axis=0).T
    x_values = sorted(list(set([ax.patches[i].get_x() for i in range(len(ax.patches))])))
    for j in range(1, len(x_values), 2):
        relevant_patches = [ax.patches[i] for i in range(len(ax.patches)) if
                            ax.patches[i].get_x() == x_values[j]]
        for p in relevant_patches:
            p.set_x(p.get_x()-0.4)
            p.set_hatch('//')
    ax.set_xticklabels([Text(0, 0, ax.get_xticklabels()[i].get_text().split('_')[0].replace('pyclone', 'PyClone').replace('sciclone', 'SciClone').replace('baseline', 'Baseline').replace('expands', 'Expands')) for i in range(0, len(ax.get_xticklabels()), 2)], fontsize=17)
    ax.set_xticks(ax.get_xticks()[range(0, len(ax.get_xticks()), 2)]+0.3)
    plt.setp(ax.get_xticklabels(), rotation=40, ha='right')
    handles, labels = ax.get_legend_handles_labels()
    ax.set_ylabel('number of samples', fontsize=17)
    ax.set_xlabel('')
    ax.legend_.remove()
    ax.set_title(loc, fontsize=30)
    if ii == 1:
        other_handles = [Patch(facecolor='white', edgecolor='gray', label='protected'),
                         Patch(facecolor='white', edgecolor='gray', label='public', hatch='//')]
        first_legend = ax.add_artist(ax.legend(handles=other_handles, bbox_to_anchor=(1., 1), title='input mutations', borderaxespad=0., fontsize=17))
        other_ = ax.legend(handles=handles, labels=labels, bbox_to_anchor=(1., 0.7), title='Number\nof clones', borderaxespad=0., fontsize=17)


ax = axes[1, 1]
df = tmp_df_c.copy()
sub_df = df.dropna(subset=nb_clones_col)
for ith_method in ['pyclone', 'CSR', 'PhyloWGS', 'sciclone', 'baseline', 'expands']:
    for folder in ['protected_hg38_vcf', 'public_hg38_vcf']:
        sub_df['{}_{}_nb_clones_bkp'.format(ith_method, folder)] = sub_df['{}_{}_nb_clones'.format(ith_method, folder)].copy()
        sub_df.loc[sub_df['{}_{}_nb_clones'.format(ith_method, folder)] >= 6, '{}_{}_nb_clones'.format(ith_method, folder)] = '6+'
aa = pd.melt(sub_df, value_vars=nb_clones_col).\
    pivot_table(index='value', columns=['variable'], aggfunc=len).fillna(0).T
aa.columns = [int(c) if type(c)!=str else c for c in aa.columns]
aa.plot.bar(stacked=True, colormap=new_cmap, edgecolor='gray', ax=ax)
ee = pd.melt(sub_df, value_vars=bkp_nb_clones_col).\
    pivot_table(index='value', columns=['variable'], aggfunc=len).fillna(0).apply(lambda x:  x/sum(x), axis=0).T
x_values = sorted(list(set([ax.patches[i].get_x() for i in range(len(ax.patches))])))
for j in range(1, len(x_values), 2):
    relevant_patches = [ax.patches[i] for i in range(len(ax.patches)) if ax.patches[i].get_x()==x_values[j]]
    for p in relevant_patches:
        p.set_x(p.get_x()-0.4)
        p.set_hatch('//')
ax.set_xticklabels([Text(0, 0, ax.get_xticklabels()[i].get_text().split('_')[0].replace('pyclone', 'PyClone').replace('sciclone', 'SciClone').replace('baseline', 'Baseline').replace('expands', 'Expands')) for i in range(0, len(ax.get_xticklabels()), 2)], fontsize=17)
ax.set_xticks(ax.get_xticks()[range(0, len(ax.get_xticks()), 2)]+0.3)
plt.setp(ax.get_xticklabels(), rotation=40, ha='right')

handles, labels = ax.get_legend_handles_labels()

ax.set_ylabel('number of samples', fontsize=17)
ax.set_xlabel('')
ax.legend_.remove()
ax.set_title('All 3 types', fontsize=30)
fig.subplots_adjust(hspace=0.4)
plt.savefig('{}/total_ith_composition.pdf'.format(output_path),
            bbox_extra_artists=(other_,), bbox_inches='tight')
plt.close()

# runtime of methods
runtime_cols = [c for c in tmp_df_c.columns if c[-4:] == 'time']

melted = pd.melt(tmp_df_c,
                 id_vars=['mutation_count_protected', 'mutation_count_public'],
                 value_vars=runtime_cols)
protected_mask = melted.variable.str.contains('protected')
melted.loc[~protected_mask, 'mutation_count_protected'] = np.nan
melted.loc[protected_mask, 'mutation_count_public'] = np.nan
melted = melted.assign(nb_mutations=melted.mutation_count_protected)
melted['nb_mutations'].fillna(melted.mutation_count_public, inplace=True)
melted = melted.assign(ith_method=pd.Categorical(melted.variable.str.split('_').str[0].str.replace('pyclone', 'PyClone').str.replace('sciclone', 'SciClone').replace('baseline', 'Baseline').replace('expands', 'Expands'),
                                                 categories=['Baseline', 'Expands', 'MATH', 'PhyloWGS', 'PyClone', 'SciClone'], ordered=True))
melted = melted.assign(value_hour=melted.value/3600)


sns.set_style('whitegrid')
g = sns.lmplot(x="nb_mutations", y="value_hour", hue="ith_method",
               truncate=True, size=5, data=melted[(melted.value_hour < 15)],
               scatter_kws={"s": 1.5}, order=2, legend=False,
               palette=colors_protected)
plt.legend(bbox_to_anchor=(1.05, 0.7), markerscale=8)

g.set_axis_labels("number of mutations", "runtime (hours)")
g.ax.set_ylim([0, 17])
plt.savefig('{}/ith_method_runtime.pdf'.format(output_path),
            bbox_inches='tight')
plt.close()


# correlation between variables
cols_plot = ['purity', 'perc_non_diploid', 'inv_T_cells', 'mutation_count_public', 'mutation_count_protected', 'math_public', 'math_protected'] + \
            [c for c in restricted_df_dict['no_purity1'] if (c[-9:]=='nb_clones') and ('norm' not in c) and ('PhyloSub' not in c)]
alt_cols_plot = [c.replace('_hg38_vcf_nb_clones', '').replace('pyclone', 'PyClone').replace('sciclone', 'SciClone').replace('baseline', 'Baseline').replace('expands', 'Expands').replace('math', 'MATH') for c in cols_plot]
color_list = [color_dict[c.replace('_', '+')] if c.replace('_', '+') in color_dict.keys() else 'white' for c in alt_cols_plot]
network_colors = pd.Series(color_list, index=alt_cols_plot)
restricted_df_dict['no_purity1'] = restricted_df_dict['no_purity1'].assign(inv_T_cells=100-restricted_df_dict['no_purity1']['T cells'])
methods = [c.split('_')[0] for c in alt_cols_plot]


#follow_labels = ['purity'] + ['{}_{}'.format(m, f) for m in ('PyClone', 'PhyloWGS', 'CSR', 'SciClone') for f in ['public', 'protected']]
follow_labels = ['{}_{}'.format(m, f) for m in ('Expands', 'mutation_count') for f in ['public', 'protected']]
mpl.rcParams['hatch.linewidth'] = 4.0  # previous svg hatch linewidth
mpl.rcParams['lines.linewidth'] = 10
for ii, loc in enumerate(cancer_locs):
    sub_df = restricted_df_dict['no_purity1'][(restricted_df_dict['no_purity1'].cancer_loc==loc)].dropna(subset=cols_restricted)
    sub_df.rename(columns=dict(zip(cols_plot, alt_cols_plot)), inplace=True)
    sns.set(font_scale=3)
    cg = sns.clustermap(sub_df[alt_cols_plot].corr(), figsize=(20, 20), cmap="vlag", vmin=-1, vmax=1, row_colors=network_colors, col_colors=network_colors, cbar_kws={"ticks":[-1, 0, 1]})
    order_labels = [t.get_text() for t in cg.ax_heatmap.get_xticklabels()]
    aa = cg.ax_row_colors
    bb = cg.ax_col_colors

    for i, n in enumerate(order_labels):
        if 'public' in n:
            aa.add_patch(mpl.patches.Rectangle([0, i], 1, 1, fill=False, hatch='/', ec='k', lw=3))
            bb.add_patch(mpl.patches.Rectangle([i, 0], 1, 1, fill=False, hatch='/', ec='k', lw=3))
        else:
            aa.add_patch(mpl.patches.Rectangle([0, i], 1, 1, fill=False, hatch='', ec='k', lw=3))
            bb.add_patch(mpl.patches.Rectangle([i, 0], 1, 1, fill=False, hatch='', ec='k', lw=3))
    cg.fig.suptitle(loc, fontsize=80)
    idx = [order_labels.index(l) for l in follow_labels]
    cg.ax_heatmap.add_patch(mpl.patches.Rectangle((min(idx), min(idx)), max(idx)-min(idx)+1, max(idx)-min(idx)+1, fill=False, edgecolor='black', lw=4))
    cg.ax_heatmap.add_patch(mpl.patches.Rectangle((0, min(idx)), min(idx), max(idx)-min(idx)+1, fill=False, edgecolor='black', lw=4, alpha=0.9, ls='--'))
    cg.ax_heatmap.add_patch(mpl.patches.Rectangle((min(idx), 0), max(idx)-min(idx)+1, min(idx), fill=False, edgecolor='black', lw=4, alpha=0.9, ls='--'))
    if ii == 0:
        legend_items = list()
        for key, v in colors_protected.items():
            if key in methods:
                legend_items.append(Patch(facecolor=v, edgecolor='black', label=key))
        legend_items.append(Patch(facecolor='#FFFFFF', edgecolor='black', label='NGS property'))

        other_handles = [Patch(facecolor='white', edgecolor='black', label='protected'),
                         Patch(facecolor='white', edgecolor='black', label='public', hatch='/', lw=3)]
        other_ = cg.ax_heatmap.add_artist(Legend(cg.ax_heatmap, handles=other_handles, bbox_to_anchor=(1.5, 0.4), fontsize=30, title='input mutations', loc=2, labels=['protected', 'public']))
        first_legend = cg.ax_heatmap.legend(handles=legend_items, bbox_to_anchor=(1.5, 1), loc=2, fontsize=30, title='ITH methods')
    plt.savefig('{}/{}_clustermap_corr_all_variables.pdf'.format(output_path, loc), bbox_inches='tight')

    plt.close()


