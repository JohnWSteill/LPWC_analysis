library(ggplot2)
library(methods)
library(reshape2)
library(here)
source(here('LPWC_Analysis_Utilities.R'))


draw_plot <- function(tgt_gene){
    SUB_391_DATA = 'Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9'
    EC_DATA_FILE <- here('Raw_Data', SUB_391_DATA, 'genes.no_mt.ec.tab')
    TPM_DATA_FILE <- here('Raw_Data', SUB_391_DATA, 'genes.no_mt.tpm.rescale.tab')
    # Used to discover parameters of 175 & 40
    # dat = tuner(dataFile=EC_DATA_FILE, sd_cv_thresh=175)
    # dat = tuner(dataFile=TPM_DATA_FILE, sd_cv_thresh=40)
    high_var_ec_data = get_high_var_genes(dataFile=EC_DATA_FILE, sd_cv_thresh=175)
    high_var_tpm_data = get_high_var_genes(dataFile=TPM_DATA_FILE, sd_cv_thresh=40)
    load(here('Unproc_Rdata','wkspace_all.RData_20180825_164425'))
    dist = as.matrix(1 - tpm_hi$solution$corr)
    rownames(dist)=rownames(high_var_tpm_data)
    clust = hclust(as.dist(dist))
    k = 0
    refine = TRUE
    while (refine){
        k = k+ 1
        groups = cutree(clust,k)
        tgt_grp_id = groups[[tgt_gene]]
        tgt_grp = groups[groups == tgt_grp_id]
        refine = length(tgt_grp) > 10
        #refine =  refine & groups[['T']] == groups[['MIXL1']] &  groups[['T']] ==  groups[['EOMES']]
        #refine =  refine & groups[['SOX17']] == groups[['CXCR4']] &  groups[['SOX17']] ==  groups[['GATA6']]
    }
    k = k-1
    groups = cutree(clust,k)
    tgt_grp_id = groups[[tgt_gene]]
    tgt_grp = groups[groups == tgt_grp_id]
    genes = attributes(tgt_grp)$names
    df = high_var_tpm_data[genes,]
    df <- melt(df)
    df$genes <- genes

    time_from_sampname <- function(s){ return(strtoi(strsplit(s, '_')[[1]][2], base=10))}
    df['t'] = sapply(as.matrix(df)[,'variable'],time_from_sampname)
    # the first sample, with time 0 is spelled differently so timefrom sampname doesnt work
    df[is.na(df['t']),'t'] <- 0
    df['lg'] = log10(df['value']+1)
    for (i in seq_along(genes)){
        g_indx = df['genes']==genes[i]
        lag = tpm_hi$solution$lags[i]
        df[g_indx,'lag_t'] = df[g_indx,'t'] + lag
    }

    ggplot(df, aes(t, lg, group=factor(genes))) + geom_smooth(aes(color=factor(genes)))
    ggplot(df, aes(lag_t, lg, group=factor(genes))) + geom_smooth(aes(color=factor(genes)), se=FALSE) +
        geom_point(aes(color=factor(genes)), alpha = .2) +
        xlab('Time(h) with Lag adjustment') +
        ylab('Log10 TPM')
}



runrun_all <-function(){
    DATA_DIR <- 'Data/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9'
    EC_DATA_FILE <- file.path(DATA_DIR, 'genes.no_mt.ec.tab')
    TPM_DATA_FILE <- file.path(DATA_DIR, 'genes.no_mt.tpm.rescale.tab')


    # Used to discover parameters of 175 & 40
    # dat = tuner(dataFile=EC_DATA_FILE, sd_cv_thresh=175)
    # dat = tuner(dataFile=TPM_DATA_FILE, sd_cv_thresh=40)
    high_var_ec_data = get_high_var_genes(dataFile=EC_DATA_FILE, sd_cv_thresh=175)
    high_var_tpm_data = get_high_var_genes(dataFile=TPM_DATA_FILE, sd_cv_thresh=40)
    try(load('wkspace_all'))
    if (!("tpm_hi" %in% ls())){
        tpm_lo =run_lo_lpwc(dat=high_var_tpm_data, timepoints = 0:96)
        save.image(file='wkspace_all')
        save.image(file=get_time_stamp_wksp_file('wkspace_all'))
    }

    try(load('wkspace_all'))
    if (!("ec_hi" %in% ls())){
        ec_hi = run_hi_lpwc(dat=high_var_ec_data, timepoints = 0:96)
        save.image(file='wkspace_all')
        save.image(file=get_time_stamp_wksp_file('wkspace_all'))
    }

    try(load('wkspace_all'))
    if (!("tpm_hi" %in% ls())){
        tpm_hi = run_hi_lpwc(dat=high_var_tpm_data, timepoints = 0:96)
        save.image(file='wkspace_all')
        save.image(file=get_time_stamp_wksp_file('wkspace_all'))
    }

    try(load('wkspace_all'))
    if (!("ec_lo" %in% ls())){
        ec_lo =run_lo_lpwc(dat=high_var_ec_data, timepoints = 0:96)
        save.image(file='wkspace_all')
        save.image(file=get_time_stamp_wksp_file('wkspace_all'))
    }
}


