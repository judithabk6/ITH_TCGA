Sources for the files in this folder

**ncomms9971-s2.xlsx**
is a supplementary table from
[Aran, D. et al. Systematic pan-cancer analysis of tumour purity. Nat. Commun. 6:8971 doi: 10.1038/ncomms9971 (2015)](https://www.nature.com/articles/ncomms9971)
and reports sample purity for the whole TCGA using different approaches. It was used for quality control on ASCAT outputs, but not further used for analyses.

**bindea13_signatures_entrez.csv**
is a supplementary table from 
[Bindea, G., Mlecnik, B., Tosolini, M., Kirilovsky, A., Waldner, M., Obenauf, A.C., Angell, H., Fredriksen, T., Lafontaine, L., Berger, A. and Bruneval, P., 2013. Spatiotemporal dynamics of intratumoral immune cells reveal the immune landscape in human cancer. Immunity, 39(4), pp.782-795.](https://www.cell.com/immunity/fulltext/S1074-7613(13)00437-8?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS1074761313004378%3Fshowall%3Dtrue)
and reports immune signatures

**ensembl_entrez.txt**
was downloaded from Ensembl biomart on March 26th 2018 to match Ensembl Gene Ids (from TCGA RNAseq data) with gene Entrez Id (from Bindea's immune signatures table)

**TCGA_mastercalls.abs_tables_JSedit.fixed.txt and TCGA_mastercalls.abs_segtabs.fixed.txt**
are results from ABSOLUTE for purity and ploidy, and allele-specific copy number for the TCGA, downloaded from [the TCGA pancancer data portal](https://gdc.cancer.gov/about-data/publications/pancanatlas) on August 18th 2019.