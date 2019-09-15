#!/usr/bin/env Rscript

library(sciClone)
library(optparse)

# get arguments
option_list <- list(
    make_option(c("-p", "--patient"), type="character", default=NULL, 
              help="patient name (format TCGA-XX-XXXX)", metavar="character"),
    make_option(c("-f", "--folder"), type="character", default=NULL, 
              help="folder indicating the filter", metavar="character")
); 
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);
if (is.null(opt$patient)){
  print_help(opt_parser)
  stop("The patient name is a mandatory argument.\n", call.=FALSE)
}
if (is.null(opt$folder)){
  print_help(opt_parser)
  stop("The folder name is a mandatory argument.\n", call.=FALSE)
}

patient <- opt$patient
folder <- opt$folder

if (grepl("absolute", folder)) {
cnv_file = read.table(paste('results', patient, 'pyclone', 'absolute_hg38_cnv_table.csv', sep='/'), sep='\t', header=TRUE)
} else {
cnv_file = read.table(paste('results', patient, 'pyclone', 'ascat_hg38_cnv_table.csv', sep='/'), sep='\t', header=TRUE)
}
#cnv_to_use = cnv_file[,c('chromosome', 'Start', 'End', 'total_cn')]
#cnv_file$abs_mean = 2^(cnv_file$Segment_Mean)*2
cnv_file$total_cn = cnv_file$major + cnv_file$minor
cnv_to_use = cnv_file[,c('chromosome', 'start', 'stop', 'total_cn')]
tt = read.table(paste('results', patient, 'pyclone', folder, 'input.tsv', sep='/'), sep='\t', header=TRUE)
hh = data.frame(do.call('rbind', strsplit(as.character(tt$mutation_id), '_')))
colnames(hh) <- c('chromosome', 'position', 'nt_change')
uu = cbind(tt,hh)
uu$vaf = uu$var_counts/(uu$var_counts + uu$ref_counts) * 100
gg  = uu[,c('chromosome', 'position', 'ref_counts', 'var_counts', 'vaf')]
gg$position = as.numeric(as.character(gg$position))
gg$chromosome = gsub('chr', '', gg$chromosome)

sc = sciClone(vafs=list(gg), copyNumberCalls=list(cnv_to_use), sampleNames=c('tumor'), minimumDepth=3, clusterMethod="bmm", clusterParams="no.apply.overlapping.std.dev.condition", cnCallsAreLog2=FALSE, useSexChrs=TRUE, doClustering=TRUE, verbose=TRUE,  copyNumberMargins=0.25, maximumClusters=10, annotation=NULL, doClusteringAlongMargins=TRUE, plotIntermediateResults=0)
out_path = paste('results', patient, 'sciclone', folder, sep='/')
dir.create(out_path, showWarnings = FALSE,recursive=TRUE)
writeClusterTable(sc, paste(out_path, 'clusters1', sep='/'))
sc.plot1d(sc,paste(out_path, "clusters1_1d.pdf", sep='/'))
