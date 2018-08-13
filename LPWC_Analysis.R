library(ggplot2)
source('LPWC_Analysis_Utilities.R')

DATA_DIR <- 'Data/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9'
EC_DATA_FILE <- file.path(DATA_DIR, 'genes.no_mt.ec.tab')
TPM_DATA_FILE <- file.path(DATA_DIR, 'genes.no_mt.tpm.rescale.tab')

# Used to discover parameters of 175 & 40
# dat = tuner(dataFile=EC_DATA_FILE, sd_cv_thresh=175)
# dat = tuner(dataFile=TPM_DATA_FILE, sd_cv_thresh=40)

high_var_ec_data = get_high_var_genes(dataFile=EC_DATA_FILE, sd_cv_thresh=175)
high_var_tpm_data = get_high_var_genes(dataFile=TPM_DATA_FILE, sd_cv_thresh=40)

c(sol_ec_hi,dist_ec_hi,clust_ec_hi)=run_hi_lpwc(dat=high_var_ec_data, 
                                                timepoints = 0:96)
save.image(file='wkspace_all')
save.image(file=paste0("wkspace_all.RData_",
                       format(Sys.time(), "%Y%m%d_%H%M%S_")))

c(sol_tpm_hi,dist_tpm_hi,clust_tpm_hi)=run_hi_lpwc(dat=high_var_tpm_data, 
                                                       timepoints = 0:96)
save.image(file='wkspace_all')
save.image(file= paste0( "wkspace_all.RData_",
                        format(Sys.time(), "%Y%m%d_%H%M%S_")))

#if (!("sol_ec_lo" %in% ls())){
c(sol_ec_lo,dist_ec_lo,clust_ec_lo)=run_lo_lpwc(dat=high_var_ec_data, 
                                                timepoints = 0:96)
#}
save.image(file='wkspace_all.RData')
save.image(file= paste0( "wkspace_all.RData_",
                        format(Sys.time(), "%Y%m%d_%H%M%S_")))
#if (!("sol_tpm_lo" %in% ls())){
    c(sol_tpm_lo,dist_tpm_lo,clust_tpm_lo)=run_lo_lpwc(dat=high_var_tpm_data, 
                                                timepoints = 0:96)
#}
save.image(file= paste0( "wkspace_all.RData_",format(Sys.time(), "%Y%m%d_%H%M%S_")))
save.image(file='wkspace_all.RData')
