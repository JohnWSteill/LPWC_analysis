
library(ggplot2)
library(methods)
library(reshape2)
library(here)
#install.packages("LPWC", repos = "http://cran.us.r-project.org")
library(LPWC)

SUB_391_DATA = 'Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9'
EC_DATA_FILE <- here('Raw_Data', SUB_391_DATA, 'genes.no_mt.ec.tab')
TPM_DATA_FILE <- here('Raw_Data', SUB_391_DATA, 'genes.no_mt.tpm.rescale.tab')
# Parameters chosen to give an "interesting" number of high-variance genes
#     I wanted all 6 of ('T', 'MIXL1', 'EOMES', 'SOX17', 'CXCR4', 'GATA6')
#     I wanted at least 400.
TPM_SD_CV_THRESH <- 40
EC_SD_CV_THRESH <- 175

draw_plot <- function(tgt_gene){

    # Used to discover parameters of 175 & 40
    # dat = tuner(dataFile=EC_DATA_FILE, sd_cv_thresh=175)
    # dat = tuner(dataFile=TPM_DATA_FILE, sd_cv_thresh=40)
    high_var_ec_data = get_high_var_genes(
        dataFile=EC_DATA_FILE, sd_cv_thresh=EC_SD_CV_THRESH)
    high_var_tpm_data = get_high_var_genes(
        dataFile=TPM_DATA_FILE, sd_cv_thresh=TPM_SD_CV_THRESH)
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


dbg_fn <- function(){
    return (list(solution=1, dist=2, clust=3))
}

get_high_var_genes <- function(dataFile, sd_cv_thresh=500){
    # The true coefficient of variance is unitless,
    # I squared the sd to get larger counts.
    # YB: Try log transforming 1st, then take highest var
    # Look for packages in Bioconductor
    my_data <- read.table(dataFile, header=TRUE, row.names=1, sep = '\t',  quote = "", comment.char = "")
    my_data <-  my_data[,-ncol(my_data)] # Drop description
    my_data$sd_coef_var <-  apply(my_data, 1, function(x) sd(x) * sd(x)/ mean(x))
    my_data$sd_coef_var[is.na(my_data$sd_coef_var)]<-0 #YB why na's?
    return(subset(my_data[my_data$sd_coef_var>sd_cv_thresh,],
                  select=-c(sd_coef_var)))
}

get_appropriate_neighbors <- function(t, degree){
    mid = round(t)
    if (mid == floor(t)){
        bot = mid - floor((degree-1)/2)
        top = mid + ceil((degree-1)/2)
    } else {
        bot = mid - ceil((degree-1)/2)
        top = mid + floor((degree-1)/2)
    }
    return(seq(bot,top))
}

get_imputed_value <- function(data_row, time, degree){
    data_t = get_appropriate_neighbors(time, degree)
}

downsample_equally_spaced_data <- function(data, n_timepoints){
    degree = int(ncol(data)/n_timepoints) + 1
}

run_hi_lpwc <- function(dat, timepoints){
    solution = LPWC::corr.bestlag(dat, timepoints = timepoints, max.lag = 20,
                                  penalty = "high", iter = 10)
    dist <- 1 - solution$corr
    distm = as.matrix(dist)
    rownames(distm) = rownames(dat)
    dist = as.dist(distm)
    clust <- hclust(dist)
    return (list(solution=solution, dist=dist, clust=clust))
}

run_lo_lpwc <- function(dat, timepoints){
    solution = LPWC::corr.bestlag(dat, timepoints = timepoints, max.lag = 20,
                                  penalty = "low", iter = 10)
    dist <- 1 - solution$corr
    distm = as.matrix(dist)
    rownames(distm) = rownames(dat)
    dist = as.dist(distm)
    clust <- hclust(dist)
    return (list(solution=solution, dist=dist, clust=clust))
}

get_h_clust <- function(solution, dat){
    distm <- as.matrix(1 - solution$corr)
    rownames(distm) = rownames(dat)
    dist = as.dist(distm)
    return(hclust(dist))
}

tuner <- function(dataFile, sd_cv_thresh){
    # for TPM, 40 gets me all 6 with only 458 tot
    # for EC, 150 gets me all 6 with 1428 tot
    genes_of_interest = c('T', 'MIXL1', 'EOMES',
                          'SOX17', 'CXCR4', 'GATA6')
    dat = get_high_var_genes(dataFile=dataFile, sd_cv_thresh=sd_cv_thresh)
    n = 0
    for (i in seq_along(genes_of_interest)){
        if (genes_of_interest[i] %in% rownames(dat)){
            n = n+ 1
            print(genes_of_interest[i])
        }
    }
    print(c(n, dim(dat)))
}

get_time_stamp_wksp_file <- function(basename){
    return(paste0("wkspace_all.RData_", format(Sys.time(), "%Y%m%d_%H%M%S")))
}

# groups = cutree(clust, k=10)
# solution_lo = LPWC::corr.bestlag(top_var, timepoints = timepoints, max.lag = 20, penalty = "low", iter = 10)

