#!/usr/bin/env Rscript


library(expands)
library(optparse)

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


out_path = paste('results', patient, 'expands', folder, sep='/')
dir.create(out_path, showWarnings = FALSE,recursive=TRUE)

tt = read.table(paste('results', patient, 'pyclone', folder, 'input.tsv', sep='/'), sep='\t', header=TRUE)
hh = data.frame(do.call('rbind', strsplit(as.character(tt$mutation_id), '_')))
colnames(hh) <- c('chromosome', 'position', 'nt_change')
hh$position = as.character(hh$position)
uu = cbind(tt,hh)
uu$AF_Tumor = uu$var_counts/(uu$var_counts+uu$ref_counts)
uu$short_chrm = gsub('chr', '', uu$chromosome)
extract = uu[(uu$short_chrm!='X')&(uu$short_chrm!='Y'),c('short_chrm', 'position', 'AF_Tumor')]
extract$PN_B=0
extract$short_chrm = as.numeric(extract$short_chrm)
colnames(extract) <- c('chr', 'startpos', 'AF_Tumor', 'PN_B')
ssm_input = data.matrix(extract)

if (grepl("absolute", folder)) {
cnv_file = read.table(paste('results', patient, 'pyclone', 'absolute_hg38_cnv_table.csv', sep='/'), sep='\t', header=TRUE)
} else {
cnv_file = read.table(paste('results', patient, 'pyclone', 'ascat_hg38_cnv_table.csv', sep='/'), sep='\t', header=TRUE)
}
cnv_file$total_cn = cnv_file$major + cnv_file$minor
extract_cnv = cnv_file[(cnv_file$chromosome!='X')&(cnv_file$chromosome!='Y'),c('chromosome', 'start', 'stop', 'total_cn')]
extract_cnv$chromosome = as.numeric(extract_cnv$chromosome)
colnames(extract_cnv) = c('chr', 'startpos', 'endpos', 'CN_Estimate')
cbs_input = data.matrix(extract_cnv)
#set.seed(6); idx=sample(1:nrow(ssm_input), 80, replace=FALSE); snv=ssm_input[idx,]
setwd(out_path)
out = runExPANdS(ssm_input, cbs_input, snvF=paste(patient, folder, sep='__'))
save(out, file='expands_output.RData')



pdf(paste(paste(patient, folder, sep='__'), 'subpopulations_raw.pdf', sep='_'), width=12, height=12)
plotSPs(out$dm, cex=1, sampleID=paste(patient, folder, sep='__'), rawAF=T)
dev.off()
pdf(paste(paste(patient, folder, sep='__'), 'subpopulations_estimated.pdf', sep='_'), width=12, height=12)
plotSPs(out$dm, cex=1, sampleID=paste(patient, folder, sep='__'), rawAF=F)
dev.off()

