#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the objective of this script is to take as input a copy number file and a vcf file and to write input for
pyclone, sciclone, phylowgs and expands
"""
import pandas as pd
from util_functions import safe_mkdir
import numpy as np


pyclone_fields = ['mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn']
cnv_cols = ['chromosome', 'start', 'stop', 'major', 'minor', 'patient_id']


def write_ith_input(patient_name, folder, vcf_filename, cnv_filename,
                    doc_normal_filename=None, doc_tumor_filename=None):
    if (patient_name == 'CO8_liver'):
        vcf = pd.read_csv(vcf_filename, sep='\t', skiprows=77, header=None,
                          names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                                 'FILTER', 'INFO', 'FORMAT', 'TUMOR', 'NORMAL'])
    else:
        vcf = pd.read_csv(vcf_filename, sep='\t', skiprows=77, header=None,
                          names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL',
                                 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR'])

    vcf = vcf.assign(ref_counts=vcf.TUMOR.str.split(':').str.get(1).str.split(',').str.get(0).astype(int))
    vcf = vcf.assign(var_counts=vcf.TUMOR.str.split(':').str.get(1).str.split(',').str.get(1).astype(int))
    vcf = vcf.assign(ref_counts_normal=vcf.NORMAL.str.split(':').str.get(1).str.split(',').str.get(0).astype(int))
    vcf = vcf.assign(var_counts_normal=vcf.NORMAL.str.split(':').str.get(1).str.split(',').str.get(1).astype(int))
    vcf['ref_l'] = vcf.REF.str.len() 
    vcf['alt_l'] = vcf.ALT.str.len()
    vcf = vcf.assign(mutation_id=vcf.CHROM + '_' + vcf.POS.astype(str) + '_' + vcf.REF + '_' + vcf.ALT)
    vcf = vcf.assign(sub_id=vcf.CHROM + '_' + vcf.POS.astype(str))
    vcf = vcf.assign(n_vaf=vcf.var_counts_normal / (vcf.var_counts_normal + vcf.ref_counts_normal))
    vcf = vcf.assign(n_depth=vcf.var_counts_normal + vcf.ref_counts_normal)
    vcf = vcf.assign(t_vaf=vcf.var_counts / (vcf.var_counts + vcf.ref_counts))
    vcf = vcf.assign(t_depth=vcf.var_counts + vcf.ref_counts)

    #if ('ALL_gawad' not in patient_name):
    cnv_table = pd.read_csv(cnv_filename, sep='\t')
    cnv_table = cnv_table.dropna(subset=['Geno'])
    cnv_table = cnv_table.assign(chromosome=cnv_table.chrom)
    cnv_table = cnv_table.assign(start=cnv_table['loc.start'])
    cnv_table = cnv_table.assign(stop=cnv_table['loc.end'])
    cnv_table = cnv_table.assign(total_cn=cnv_table['CNt'].astype(int))
    cnv_table = cnv_table.assign(major=cnv_table['Geno'].str.count('A').astype(int))
    cnv_table = cnv_table.assign(minor=cnv_table['Geno'].str.count('B').astype(int))
    cnv_table = cnv_table.assign(patient_id=cnv_table['ID'])

    cnv_table = cnv_table.assign(full_chrom='chr' + cnv_table.chrom.astype(str))

    normal_doc_raw = pd.read_csv(doc_normal_filename, sep='\t')
    normal_doc = normal_doc_raw[normal_doc_raw.Total_Depth > 0]
    baselist = 'ACGTND'
    base_count_normal_colname = [c for c in normal_doc.columns if 'base_counts' in c][0] 
    for i, base in enumerate(baselist):
        normal_doc = normal_doc.assign(**{'{}_count'.format(base): normal_doc[base_count_normal_colname].str.split(' ').str.get(i).str.split(':').str.get(1).astype(int)})
    count_cols = ['{}_count'.format(base) for base in baselist]
    normal_doc = normal_doc.assign(ref_base=np.array(list(baselist))[normal_doc[count_cols].values.argmax(axis=1)])
    normal_doc = normal_doc.assign(ref_counts_normal=normal_doc[count_cols].values.max(axis=1))

    tumor_doc = pd.read_csv(doc_tumor_filename, sep='\t')
    tumor_doc = tumor_doc[normal_doc_raw.Total_Depth > 0]
    base_count_tumor_colname = [c for c in tumor_doc.columns if 'base_counts' in c][0]
    for i, base in enumerate(baselist):
        tumor_doc = tumor_doc.assign(**{'{}_count'.format(base): tumor_doc[base_count_tumor_colname].str.split(' ').str.get(i).str.split(':').str.get(1).astype(int)})
    tumor_doc = tumor_doc.assign(ref_base=normal_doc.ref_base)
    tumor_doc = tumor_doc.assign(ref_counts=tumor_doc.apply(lambda x: x['{}_count'.format(x['ref_base'])], axis=1))
    for idx, row in tumor_doc.iterrows():
        tumor_doc.loc[idx, '{}_count'.format(row['ref_base'])] = 0
    tumor_doc = tumor_doc.assign(var_base=np.array(list(baselist))[tumor_doc[count_cols].values.argmax(axis=1)])
    tumor_doc = tumor_doc.assign(var_counts=tumor_doc[count_cols].values.max(axis=1))

    normal_doc = normal_doc.assign(var_base=tumor_doc.var_base)
    normal_doc = normal_doc.assign(var_counts_normal=normal_doc.apply(lambda x: x['{}_count'.format(x['var_base'])], axis=1))

    tumor_doc = tumor_doc.assign(ref_counts_normal=normal_doc.ref_counts_normal)
    tumor_doc = tumor_doc.assign(var_counts_normal=normal_doc.var_counts_normal)
    tumor_doc = tumor_doc.assign(t_depth=tumor_doc.var_counts + tumor_doc.ref_counts)
    tumor_doc = tumor_doc.assign(n_depth=tumor_doc.var_counts_normal + tumor_doc.ref_counts_normal)
    tumor_doc = tumor_doc.assign(t_vaf=tumor_doc.var_counts / tumor_doc.t_depth)
    tumor_doc = tumor_doc.assign(n_vaf=tumor_doc.var_counts_normal / tumor_doc.n_depth)
    tumor_doc = tumor_doc.assign(ref_l=1)
    tumor_doc = tumor_doc.assign(alt_l=1)
    tumor_doc = tumor_doc.assign(sub_id=tumor_doc.Locus.str.replace(':', '_'))
    tumor_doc = tumor_doc.assign(CHROM=tumor_doc.Locus.str.split(':').str.get(0))
    tumor_doc = tumor_doc.assign(POS=tumor_doc.Locus.str.split(':').str.get(1).astype(int))
    tumor_doc = tumor_doc.assign(mutation_id=tumor_doc.Locus.str.replace(':', '_') + '_' + tumor_doc.ref_base + '_' + tumor_doc.var_base)

    vcf_c = pd.concat((vcf, tumor_doc), sort=True)
    vcf_c_dedup = vcf_c.drop_duplicates(subset=['sub_id'], keep='first')


    # merge
    vcf_cnv = pd.merge(vcf_c_dedup, cnv_table, left_on='CHROM', right_on='full_chrom')
    vcf_cnv_filter = vcf_cnv[(vcf_cnv.POS >= vcf_cnv.start)&(vcf_cnv.POS < vcf_cnv.stop)]
    vcf_cnv_filter = vcf_cnv_filter.assign(normal_cn=2)
    vcf_cnv_filter = vcf_cnv_filter.assign(minor_cn=vcf_cnv_filter.minor.astype(int))
    vcf_cnv_filter = vcf_cnv_filter.assign(major_cn=vcf_cnv_filter.major.astype(int))

    # filter to write mutations
    vcf_cnv_filter_cov = vcf_cnv_filter[(vcf_cnv_filter.n_depth>=6) &
                   ((vcf_cnv_filter.var_counts_normal<=1) | (vcf_cnv_filter.n_vaf<=0.01)) & 
                   (vcf_cnv_filter.t_depth>=8) &
                   ((vcf_cnv_filter.var_counts>=3) | (vcf_cnv_filter.t_vaf>0.2)) &
                   (vcf_cnv_filter.n_depth>=6) &
                   ((vcf_cnv_filter.var_counts_normal<=1) | (vcf_cnv_filter.n_vaf<=0.01)) & 
                   (vcf_cnv_filter.t_depth>=8) &
                   ((vcf_cnv_filter.var_counts>=3) | (vcf_cnv_filter.t_vaf>0.2)) &
                   (vcf_cnv_filter.ref_l == 1) & (vcf_cnv_filter.alt_l == 1) &
                   (vcf_cnv_filter.total_cn > 0)]

    safe_mkdir('results/{}'.format(patient_name))
    safe_mkdir('results/{}/pyclone'.format(patient_name))
    safe_mkdir('results/{}/pyclone/{}'.format(patient_name, folder))
    cnv_table[cnv_cols].to_csv('results/{}/pyclone/ascat_hg38_cnv_table.csv'.format(patient_name),
                               sep='\t', index=False)
    try:
        patient_purity = float(cnv_filename.split('allele_spe_cnv_')[1].split('cellularity')[0])
    except:
        patient_purity = 1

    with open('results/{}/pyclone/ascat_purity.txt'.format(patient_name), 'w') as f:
        f.write('{}'.format(patient_purity))
    print(patient_name, len(vcf_cnv_filter_cov))
    vcf_cnv_filter_cov[pyclone_fields].to_csv('results/{}/pyclone/{}/input.tsv'
                                            .format(patient_name, folder), sep='\t', index=False)

# TNBC
folder = 'protected_hg38_vcf'
vcf_filename = 'single_cell_dataset/BRCA_wang/ewok_results/BRCA_wang/Variant_Calling/Mutect/BRCA_wang.tumor_constit.mutect.pass.vcf'
cnv_filename = 'single_cell_dataset/BRCA_wang/ewok_results/BRCA_wang/Copy_Number/Facets/BRCA_wang_subclonal_allele_spe_cnv_0.7cellularity_2ploidy.transformed.txt'
patient_name = 'BRCA_wang_TNBC'
doc_normal_filename = 'single_cell_dataset/BRCA_wang/ewok_results/BRCA_wang/doc_testBRCA_wang_constit.doc'
doc_tumor_filename = 'single_cell_dataset/BRCA_wang/ewok_results/BRCA_wang/doc_testBRCA_wang_tumor.doc'
write_ith_input(patient_name, folder, vcf_filename, cnv_filename, doc_normal_filename, doc_tumor_filename)

# ALL_gawad p1
folder = 'protected_hg38_vcf'
vcf_filename = 'single_cell_dataset/ALL_gawad/ewok_results/ALL_gawad/patient1-Tumor/Variant_Calling/Mutect/patient1-Tumor.tumor_constit.mutect.pass.vcf'
cnv_filename = 'single_cell_dataset/ALL_gawad/ewok_results/patient1-Tumor/ALL_gawad_patient1/Copy_Number/Facets/ALL_gawad_patient1_subclonal_allele_spe_cnv_NAcellularity_2ploidy.transformed.txt'
doc_normal_filename = 'single_cell_dataset/ALL_gawad/doc_SRR1517761.doc'
doc_tumor_filename = 'single_cell_dataset/ALL_gawad/doc_SRR1517762.doc'
patient_name = 'ALL_gawad_P1'
write_ith_input(patient_name, folder, vcf_filename, cnv_filename, doc_normal_filename, doc_tumor_filename)

# # ALL_gawad p2
folder = 'protected_hg38_vcf'
vcf_filename = 'single_cell_dataset/ALL_gawad/ewok_results/patient2-Tumor/ALL_gawad_patient2/Variant_Calling/Mutect/ALL_gawad_patient2.tumor_constit.mutect.pass.vcf'
cnv_filename = 'single_cell_dataset/ALL_gawad/ewok_results/patient2-Tumor/ALL_gawad_patient2/Copy_Number/Facets/ALL_gawad_patient2_subclonal_allele_spe_cnv_0.49cellularity_2ploidy.transformed.txt'
patient_name = 'ALL_gawad_P2'
doc_normal_filename = 'single_cell_dataset/ALL_gawad/doc_SRR1517763.doc'
doc_tumor_filename = 'single_cell_dataset/ALL_gawad/doc_SRR1517764.doc'
write_ith_input(patient_name, folder, vcf_filename, cnv_filename, doc_normal_filename, doc_tumor_filename)

# CO5_colon
folder = 'protected_hg38_vcf'
vcf_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO5/CO5_colon-Tumor/Variant_Calling/Mutect/CO5_colon-Tumor.tumor_constit.mutect.pass.vcf'
cnv_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO5/CO5_colon-Tumor/Copy_Number/Facets/CO5_colon-Tumor_subclonal_allele_spe_cnv_0.55cellularity_3ploidy.transformed.txt'
patient_name = 'CO5_colon'
doc_normal_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO5/doc_testSRR3472566.doc'
doc_tumor_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO5/doc_testSRR3472569.doc'
write_ith_input(patient_name, folder, vcf_filename, cnv_filename, doc_normal_filename, doc_tumor_filename)

# CO5_liver
folder = 'protected_hg38_vcf'
vcf_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO5/CO5_liver-Tumor/Variant_Calling/Mutect/CO5_liver-Tumor.tumor_constit.mutect.pass.vcf'
cnv_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO5/CO5_liver-Tumor/Copy_Number/Facets/CO5_liver-Tumor_subclonal_allele_spe_cnv_0.8cellularity_3ploidy.transformed.txt'
patient_name = 'CO5_liver'
doc_normal_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO5/doc_testSRR3472567.doc'
doc_tumor_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO5/doc_testSRR3472571.doc'
write_ith_input(patient_name, folder, vcf_filename, cnv_filename, doc_normal_filename, doc_tumor_filename)

# CO8_colon
folder = 'protected_hg38_vcf'
vcf_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO8/CO8_colon-Tumor/Variant_Calling/Mutect/CO8_colon-Tumor.tumor_constit.mutect.pass.vcf'
cnv_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO8/CO8_colon-Tumor/Copy_Number/Facets/CO8_colon-Tumor_subclonal_allele_spe_cnv_0.89cellularity_4ploidy.transformed.txt'
patient_name = 'CO8_colon'
doc_normal_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO8/doc_testSRR3472798.doc'
doc_tumor_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO8/doc_testSRR3472800.doc'
write_ith_input(patient_name, folder, vcf_filename, cnv_filename, doc_normal_filename, doc_tumor_filename)

# CO8_liver
folder = 'protected_hg38_vcf'
vcf_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO8/CO8_liver-Tumor/Variant_Calling/Mutect/CO8_liver-Tumor.tumor_constit.mutect.pass.vcf'
cnv_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO8/CO8_liver-Tumor/Copy_Number/Facets/CO8_liver-Tumor_subclonal_allele_spe_cnv_0.93cellularity_3ploidy.transformed.txt'
patient_name = 'CO8_liver'
doc_normal_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO8/doc_testSRR3472799.doc'
doc_tumor_filename = 'single_cell_dataset/CRC_leung/ewok_results/CRC_leung_CO8/doc_testSRR3472796.doc'
write_ith_input(patient_name, folder, vcf_filename, cnv_filename, doc_normal_filename, doc_tumor_filename)

























