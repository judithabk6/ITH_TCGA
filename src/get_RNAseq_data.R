#!/usr/bin/env Rscript



args = commandArgs(trailingOnly=TRUE)


require(TCGAbiolinks)
require(DESeq2)
require(dplyr)
require(ggplot2)

#missing location commands


#for (loc in c('BRCA', 'BLCA', 'HNSC')) {
loc=args[1]
setwd(paste('data', loc, 'RNAseq', sep='/'))
query <- GDCquery(project = paste0('TCGA-', loc),
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts")

GDCdownload(query,method='client')
data <- GDCprepare(query,directory='GDCdata')

require(SummarizedExperiment)
dim(assay(data))
write.table(assay(data),file=paste0('data_TCGA_', loc, '.csv'), row.names=TRUE, col.names=TRUE,sep='\t', quote=FALSE)
write.table(data.frame(rowData(data)),file='gene_Ensemble_Symbol.csv', row.names=TRUE, col.names=TRUE,sep='\t', quote=FALSE)
write.table(data.frame(colData(data)),file='sample_Ensemble_Symbol.csv', row.names=TRUE, col.names=TRUE,sep='\t', quote=FALSE)

data_tcga <- read.table(paste0('data_TCGA_', loc, '.csv'), sep='\t', header=TRUE, check.names=FALSE)
gene_tcga <- read.table('gene_Ensemble_Symbol.csv', sep='\t', header=TRUE, check.names=FALSE)
sampleinfo_tcga <- read.table('sample_Ensemble_Symbol.csv', sep='\t', header=TRUE, check.names=FALSE)

tcga_obj <- DESeqDataSetFromMatrix(countData=data_tcga, colData=sampleinfo_tcga,design=~sample)
# vst/rlog calculate size factors internally
tcga_vst_obj <- varianceStabilizingTransformation(tcga_obj)

colnames(data_tcga) <- gsub('\\.','-',colnames(data_tcga))
colnames(tcga_vst_obj) <- colnames(data_tcga)

save(tcga_vst_obj,file = 'tcga_vst_obj.Rdata')
write.table(assay(tcga_vst_obj),file = 'tcga_vst_matrix.csv', row.names=TRUE, col.names=TRUE,sep='\t', quote=FALSE)



