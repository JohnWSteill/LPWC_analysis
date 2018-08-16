library(ggplot2)
library(methods)
source('LPWC_Analysis_Utilities.R')

DATA_DIR <- 'Data/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9'
EC_DATA_FILE <- file.path(DATA_DIR, 'genes.no_mt.ec.tab')
TPM_DATA_FILE <- file.path(DATA_DIR, 'genes.no_mt.tpm.rescale.tab')

try(load('wkspace_all'))

# Used to discover parameters of 175 & 40
# dat = tuner(dataFile=EC_DATA_FILE, sd_cv_thresh=175)
# dat = tuner(dataFile=TPM_DATA_FILE, sd_cv_thresh=40)
high_var_ec_data = get_high_var_genes(dataFile=EC_DATA_FILE, sd_cv_thresh=175)
high_var_tpm_data = get_high_var_genes(dataFile=TPM_DATA_FILE, sd_cv_thresh=40)

if (!("ec_hi" %in% ls())){
    ec_hi = run_hi_lpwc(dat=high_var_ec_data, timepoints = 0:96)
    save.image(file='wkspace_all')
    save.image(file=get_time_stamp_wksp_file('wkspace_all'))
}

if (!("tpm_hi" %in% ls())){
    tpm_hi = run_hi_lpwc(dat=high_var_tpm_data, timepoints = 0:96) 
    save.image(file='wkspace_all')
    save.image(file=get_time_stamp_wksp_file('wkspace_all'))
}

if (!("ec_lo" %in% ls())){
    ec_lo =run_lo_lpwc(dat=high_var_ec_data, timepoints = 0:96)
    save.image(file='wkspace_all')
    save.image(file=get_time_stamp_wksp_file('wkspace_all'))
}

if (!("tpm_hi" %in% ls())){ 
    tpm_lo =run_lo_lpwc(dat=high_var_tpm_data, timepoints = 0:96)
    save.image(file='wkspace_all')
    save.image(file=get_time_stamp_wksp_file('wkspace_all'))
}

