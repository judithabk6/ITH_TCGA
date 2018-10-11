#!/usr/bin/env Rscript

library(survcomp)
#for (eval_setting in c('no_purity1', 'no_purity1_noexpands')) {
for (eval_setting in c('no_purity1')) {
    for (loc in c('BRCA', 'BLCA', 'HNSC')) {
        for (fold in seq(1,5)) {
            filename = paste('tmp/20180801_', eval_setting, '_test_results/', loc, '_fold', fold, '_test_res.csv', sep='')
            results = read.table(filename, header=TRUE, sep='\t')
            ground_truth_time = results[1,8:dim(results)[2]]
            ground_truth_event = results[2,8:dim(results)[2]]

            my_concordance_index <- function(x) {
                c<-concordance.index(-as.numeric(as.vector(x)), surv.time=as.numeric(as.vector(ground_truth_time)), surv.event=as.numeric(as.vector(ground_truth_event)), method='noether', alternative='g')
                #return(c(c$c.index, c$p.value))
                return(list(c))
            }
            my_concordance_index_todf <- function(x) {
                c<-concordance.index(-as.numeric(as.vector(x)), surv.time=as.numeric(as.vector(ground_truth_time)), surv.event=as.numeric(as.vector(ground_truth_event)), method='noether', alternative='g')
                #return(c(c$c.index, c$p.value))
                return(c(c$c.index, c$p.value))
            }
            true_results = results[3:dim(results)[1],]
            # concordance.index(x=age, surv.time=stime, surv.event=sevent, strat=strat,
            #   weights=weight, method="noether", comppairs=comppairs)
            all_ci_tmp = apply(true_results[8:dim(results)[2]], 1, my_concordance_index)
            u = apply(true_results[8:dim(results)[2]], 1, my_concordance_index_todf)
            u = t(u)
            colnames(u) = c('test_score_cindex_r', 'pval_test_cindex_r')
            new_df = cbind(true_results, u)[,c(seq(1,7), seq(dim(results)[2]+1, dim(results)[2]+2))]
            save_filename = paste('tmp/20180801_', eval_setting, '_test_results/', loc, '_fold', fold, '_test_res_cindex_r.csv', sep='')
            write.table(new_df, save_filename, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
            if (fold==1) {
                all_ci = all_ci_tmp
            }
            else {
                all_ci = mapply(append, all_ci, all_ci_tmp, SIMPLIFY=FALSE)
            }
        }

        comparison_test_df <- data.frame(ith_method=character(),
                                              folder=character(),
                                              setting1=character(),
                                              setting2=character(),
                                              pvalue_fold1=double(),
                                              pvalue_fold2=double(),
                                              pvalue_fold3=double(),
                                              pvalue_fold4=double(),
                                              pvalue_fold5=double(),
                                              stringsAsFactors=FALSE)
        methods = c('baseline', 'pyclone', 'PhyloWGS', 'math_score', 'combination', 'sciclone', 'expands', 'CSR')
        folder = c('protected_hg38_vcf', 'public_hg38_vcf')
        clinical_cindex_index = rownames(results[(results$ith_method=='none')&(results$folder=='none')&(results$additional_cols=='clonality+clinical'),])
        clinical_ngs_cindex_index = rownames(results[(results$ith_method=='none')&(results$folder=='none')&(results$additional_cols=='clonality+clinical+NGS'),])
        for (m in methods) {
            for (f in folder) {
                # nb clones vs clonality
                clonality_cindex_index = rownames(results[(results$ith_method==m)&(results$folder==f)&(results$additional_cols=='clonality'),])
                nb_clones_cindex_index = rownames(results[(results$ith_method==m)&(results$folder==f)&(results$additional_cols=='nb_clones'),])
                clonality_squared_cindex_index = rownames(results[(results$ith_method==m)&(results$folder==f)&(results$additional_cols=='clonality_squared'),])
                clonality_clinical_cindex_index = rownames(results[(results$ith_method==m)&(results$folder==f)&(results$additional_cols=='clonality+clinical'),])
                clonality_clinical_ngs_cindex_index = rownames(results[(results$ith_method==m)&(results$folder==f)&(results$additional_cols=='clonality+clinical+NGS'),])
                if (length(clonality_cindex_index)&length(nb_clones_cindex_index)) {
                    ll <- list()
                    for (i in seq(1:length(all_ci[[clonality_cindex_index]]))) {
                        ll[[i]] = cindex.comp(all_ci[[clonality_cindex_index]][[i]], all_ci[[nb_clones_cindex_index]][[i]])$p.value
                    }
                    comp <- cindex.comp.meta(all_ci[[clonality_cindex_index]], all_ci[[nb_clones_cindex_index]])
                    comparison_test_df <- rbind(comparison_test_df, data.frame(ith_method=m, folder=f, setting1='clonality', setting2='nb_clones', pvalue_fold1=ll[[1]], pvalue_fold2=ll[[2]], pvalue_fold3=ll[[3]], pvalue_fold4=ll[[4]], pvalue_fold5=ll[[5]]))
                }
                if (length(clonality_cindex_index)&length(nb_clones_cindex_index)) {
                    ll <- list()
                    for (i in seq(1:length(all_ci[[clonality_cindex_index]]))) {
                        ll[[i]] = cindex.comp(all_ci[[nb_clones_cindex_index]][[i]], all_ci[[clonality_cindex_index]][[i]])$p.value
                    }
                    comparison_test_df <- rbind(comparison_test_df, data.frame(ith_method=m, folder=f, setting1='nb_clones', setting2='clonality', pvalue_fold1=ll[[1]], pvalue_fold2=ll[[2]], pvalue_fold3=ll[[3]], pvalue_fold4=ll[[4]], pvalue_fold5=ll[[5]]))
                }
                if (length(clonality_cindex_index)&length(clonality_squared_cindex_index)) {
                    ll <- list()
                    for (i in seq(1:length(all_ci[[clonality_cindex_index]]))) {
                        ll[[i]] = cindex.comp(all_ci[[clonality_squared_cindex_index]][[i]], all_ci[[clonality_cindex_index]][[i]])$p.value
                    }
                    comparison_test_df <- rbind(comparison_test_df, data.frame(ith_method=m, folder=f, setting1='clonality_squared', setting2='clonality', pvalue_fold1=ll[[1]], pvalue_fold2=ll[[2]], pvalue_fold3=ll[[3]], pvalue_fold4=ll[[4]], pvalue_fold5=ll[[5]]))
                }
                if (length(clinical_cindex_index)&length(clonality_clinical_cindex_index)) {
                    ll <- list()
                    for (i in seq(1:length(all_ci[[clonality_cindex_index]]))) {
                        ll[[i]] = cindex.comp(all_ci[[clonality_clinical_cindex_index]][[i]], all_ci[[clinical_cindex_index]][[i]])$p.value
                    }
                    comparison_test_df <- rbind(comparison_test_df, data.frame(ith_method=m, folder=f, setting1='clonality+clinical', setting2='clinical_alone', pvalue_fold1=ll[[1]], pvalue_fold2=ll[[2]], pvalue_fold3=ll[[3]], pvalue_fold4=ll[[4]], pvalue_fold5=ll[[5]]))
                }
                if (length(clinical_ngs_cindex_index)&length(clonality_clinical_ngs_cindex_index)) {
                    ll <- list()
                    for (i in seq(1:length(all_ci[[clonality_cindex_index]]))) {
                        ll[[i]] = cindex.comp(all_ci[[clonality_clinical_ngs_cindex_index]][[i]], all_ci[[clinical_ngs_cindex_index]][[i]])$p.value
                    }
                    comparison_test_df <- rbind(comparison_test_df, data.frame(ith_method=m, folder=f, setting1='clonality+clinical+NGS', setting2='clinical+ngs_alone', pvalue_fold1=ll[[1]], pvalue_fold2=ll[[2]], pvalue_fold3=ll[[3]], pvalue_fold4=ll[[4]], pvalue_fold5=ll[[5]]))
                }
            }
        }
        save_filename = paste('tmp/20180801_', eval_setting, '_test_results/', loc, '_comparison_cindex_r.csv', sep='')
        write.table(comparison_test_df, save_filename, sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
    }
}