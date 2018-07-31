library(ggplot2)
source('LPWC_Analysis_Utilities.R')
if (file.exists('wkspace_lo.RData')){
    load(file='wkspace_lo.RData')
}
if (!("solution_lo" %in% ls())){
    dataFile = 'Data/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9/genes.no_mt.ec.tab' 
    high_var_data = get_high_var_genes(infile=dataFile, sd_cv_thresh=500)
    c(solution_lo,dist_lo,clust_lo)=run_lo_lpwc(dat=high_var_data, timepoints = 0:96)
}
if (!("sol_tpm_lo" %in% ls())){
    dataFile = 'Data/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9/genes.no_mt.tpm.rescale.tab' 
    high_var_data = get_high_var_genes(infile=dataFile, sd_cv_thresh=500)
    c(sol_tpm_lo,dist_tpm_lo,clust_tpm_lo)=run_lo_lpwc(dat=high_var_data, timepoints = 0:96)
}
save.image(file='wkspace_lo.RData')