#!/usr/bin/env Rscript

#------------------------------------------------------------#
# This script takes care of the preprocessing for CSR
# The input files for each sample should be in one directory
# Different samples should be in different folder
# Initialized by Kaixian Yu
# Date: 04/04/2018
# Email: kaixiany@gmail.com
#------------------------------------------------------------#
# The script takes two commandline arguments: input_dir output_dir mutation_size filtering_cutoff
args = commandArgs(trailingOnly=TRUE)
# debug use
# args <- c('/workspace/kyu2/CSR/sample1/', '/workspace/kyu2/CSR/CSR_File/sample1', '3800', '0.35')
input.dir       <- args[1]
output.dir      <- args[2]
downsample.size <- as.numeric(args[3])
filtering.cut   <- as.numeric(args[4])
suppressMessages(library(dummies))
# if the input directory does not exist then quit
if(!dir.exists(input.dir)){
    stop(sprintf('The Input directory %s does not exist.',input.dir))
}
# if the output directory does not exist, then try to create the full path
if(!dir.exists(output.dir)){
    dir.create(output.dir, recursive = TRUE)
}

# list all input files
files        <- list.files(input.dir, full.names = TRUE)
inputData    <- lapply(files,function(x){
                    tmp <- read.table(x, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
                    return(tmp)
                    })

mut.names    <- lapply(inputData,function(x){
                    return(paste(x[,1], x[,2], sep = '_'))
                    }) 
# construct mutations to be included in the consensus
mut.freq         <- table(unlist(mut.names))/length(mut.names)
mut.selected     <- sort(mut.freq[mut.freq >= filtering.cut], decreasing = TRUE) 
if(length(mut.selected) > downsample.size){
    tmp.cut      <- mut.selected[downsample.size]
    tmp.included <- names(mut.selected)[mut.selected > tmp.cut]
    mut.included <- sort(union(tmp.included, sample(names(mut.selected)[mut.selected == tmp.cut], downsample.size - length(tmp.included))))
} else {
    mut.included <- sort(names(mut.selected))
}

# construct consensus number of cluster
ncls.methods <- unlist(lapply(1:length(inputData), function(x){
                return(sum(!is.na(unique(inputData[[x]][match(mut.included,mut.names[[x]])  , 3]))))
                }))
NCL          <- ceiling(median(ncls.methods, na.rm = TRUE))

# construct the assignment file
mut.assign <- lapply(1:length(inputData), function(x){
               tmp.mut.int <- intersect(mut.included, mut.names[[x]])
               tmp.assign  <- dummy(inputData[[x]][match(tmp.mut.int, mut.names[[x]]) ,3])
               tmp <- matrix(0, nrow = length(mut.included), ncol = ncol(tmp.assign))
               rownames(tmp) <- mut.included
               tmp[tmp.mut.int,] <- tmp.assign
               return(tmp)
    })

# construct the CP file
mut.CP       <- Reduce('cbind',lapply(1:length(inputData), function(x){
                    tmp <- inputData[[x]][match(mut.included, mut.names[[x]]) ,4]
                    return(tmp)
    }))
consensus.CP <- apply(mut.CP, 1, function(x){median(x, na.rm = TRUE)})

# taking the method names for assignment output file
out.name <- lapply(files,function(x){
                tmp <- unlist(strsplit(x,'/'))
                return(tmp[length(tmp)])
    })

# output all processed files
write.table(mut.included, file = sprintf('%s/mutations_list.tsv', output.dir), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(NCL, file = sprintf('%s/number_cluster.tsv', output.dir), quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(consensus.CP, file = sprintf('%s/CellularPrevalence.tsv', output.dir), quote = FALSE, row.names = FALSE, col.names = FALSE)
invisible(lapply(1:length(inputData),function(x){
    write.table(mut.assign[[x]], file = sprintf('%s/%s_mutation_assignment.tsv', output.dir, out.name[[x]]), quote = FALSE, row.names = FALSE, col.names = FALSE, sep ='\t')
    }))
