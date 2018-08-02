library(ggplot2)
source('LPWC_Analysis_Utilities.R')

EC_DATA_FILE = 'Data/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9/genes.no_mt.ec.tab' 
TPM_DATA_FILE = 'Data/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9/genes.no_mt.tpm.rescale.tab'

if (file.exists('wkspace_lo.RData')){
    load(file='wkspace_lo.RData')
}
if (!("sol_ec_hi" %in% ls())){
    high_var_data = get_high_var_genes(infile=EC_DATA_FILE, sd_cv_thresh=500)
    c(sol_ec_hi,dist_ec_hi,clust_ec_hi)=run_hi_lpwc(dat=high_var_data, 
                                                timepoints = 0:96)
}
save.image(file='wkspace_lo.RData')
save.image(file=paste0("wkspace_lo.RDatadata_set.csv",
                       format(Sys.time(), "%Y%m%d_%H%M%S_")))
if (!("sol_tpm_hi" %in% ls())){
    high_var_data = get_high_var_genes(infile=TPM_DATA_FILE, sd_cv_thresh=500)
    c(sol_tpm_hi,dist_tpm_hi,clust_tpm_hi)=run_hi_lpwc(dat=high_var_data, 
                                                       timepoints = 0:96)
}
save.image(file='wkspace_lo.RData')
save.image(file= paste0( "wkspace_lo.RDatadata_set.csv",
                        format(Sys.time(), "%Y%m%d_%H%M%S_")))
if (!("sol_ec_lo" %in% ls())){
    c(sol_ec_lo,dist_ec_lo,clust_ec_lo)=run_lo_lpwc(dat=high_var_data, 
                                                timepoints = 0:96)
}
save.image(file='wkspace_lo.RData')
save.image(file= paste0( "wkspace_lo.RDatadata_set.csv",
                        format(Sys.time(), "%Y%m%d_%H%M%S_")))
if (!("sol_tpm_lo" %in% ls())){
    c(sol_tpm_lo,dist_tpm_lo,clust_tpm_lo)=run_lo_lpwc(dat=high_var_data, 
                                                timepoints = 0:96)
}
save.image(file= paste0( "wkspace_lo.RDatadata_set.csv",format(Sys.time(), "%Y%m%d_%H%M%S_")))
save.image(file='wkspace_lo.RData')
