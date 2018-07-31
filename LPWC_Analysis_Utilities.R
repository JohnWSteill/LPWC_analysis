#install.packages("LPWC", repos = "http://cran.us.r-project.org")
library(LPWC)

get_high_var_genes <- function(infile, sd_cv_thresh=500){
    # The true coefficient of variance is unitless, I squared the sd to get 
    # larger counts.
    # dataFile = 'Data/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm_\
    #                                      _32a860ee85db9ac9/genes.no_mt.ec.tab'
    data = read.table(dataFile, header=TRUE, row.names=1, sep = '\t')  
    my_data = my_data[,-ncol(my_data)] # Drop description
    my_data$coefVar = apply(my_data, 1, function(x) sd(x) * sd(x)/ mean(x))
    my_data$coefVar[is.na(my_data$coefVar)]<-0
    return( subset(my_data[my_data$coefVar>500,],select=-c(coefVar)))
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
    rownames(distm) = rownames(top_var) 
    dist = as.dist(distm) 
    clust <- hclust(dist)
    return (list(solution=solution, dist=dist, clust=clust))
}
run_lo_lpwc <- function(dat, timepoints){ 
    solution = LPWC::corr.bestlag(dat, timepoints = timepoints, max.lag = 20, 
                                  penalty = "low", iter = 10) 
    dist <- 1 - solution$corr 
    distm = as.matrix(dist) 
    rownames(distm) = rownames(top_var) 
    dist = as.dist(distm) 
    clust <- hclust(dist)
    return (list(solution=solution, dist=dist, clust=clust))
}
get_high_var_genes <- function(dataFile, sd_cv_thresh=500){
    # The true coefficient of variance is unitless, 
    # I squared the sd to get larger counts. 
    my_data = read.table(dataFile, header=TRUE, row.names=1, sep = '\t')
    my_data = my_data[,-ncol(my_data)] # Drop description
    my_data$sd_coef_var = apply(my_data, 1, function(x) sd(x) * sd(x)/ mean(x))
    my_data$sd_coef_var[is.na(my_data$sd_coef_var)]<-0
    return(subset(my_data[my_data$sd_coef_var>sd_cv_thresh,], 
                  select=-c(sd_coef_var)))
}
if (!interactive()){
    if (file.exists('wkspace_lo.RData')){
            load(file='wkspace_lo.RData')
    }
    if (!("solution_lo" %in% ls())){
            dataFile = 'Data/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9/genes.no_mt.ec.tab' 
        high_var_data = get_high_var_genes(dataFile=dataFile, sd_cv_thresh=500)
            c(solution_lo,dist_lo,clust_lo)=run_lo_lpwc(dat=high_var_data, timepoints = 0:96)
    }
    if (!("sol_tpm_lo" %in% ls())){
        dataFile = 'Data/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9/genes.no_mt.tpm.rescale.tab' 
        high_var_data = get_high_var_genes(dataFile=dataFile, sd_cv_thresh=500)
        c(sol_tpm_lo,dist_tpm_lo,clust_tpm_lo)=run_lo_lpwc(dat=high_var_data, timepoints = 0:96)
    }
    save.image(file='wkspace_lo.RData')
    print("finished.")
}

#c(lat, lng) %<-% list(38.061944, -122.643889)
# '''
# solution = LPWC::corr.bestlag(top_var, timepoints = timepoints, max.lag = 20, penalty = "high", iter = 10)
# 
# dist <- 1 - solution$corr
# distm = as.matrix(dist)
# rownames(distm) = rownames(top_var)
# dist = as.dist(distm)
# clust <- hclust(dist)
# groups = cutree(clust, k=10)
# save.image(file='wkspace2.RData')
# #plot(clust, cex = .5)
# # ---------------  END OF SCRIPT --------------------------------------------

# 
# 
# load(file='wkspace.RData')
# #library(ggplot2)
# # dataFile = '/isiseqruns/.GUP_HOME/RUNS/C7DCVACXX/Collate/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9/genes.no_mt.ec.tab'

# my_data = read.table(dataFile, header=TRUE, row.names=1, sep = '\t')
# 
#  
# 
# timepoints = 0:96
# solution_lo = LPWC::corr.bestlag(top_var, timepoints = timepoints, max.lag = 20, penalty = "low", iter = 10)
# dist_lo <- 1 - solution_lo$corr
# distm_lo = as.matrix(dist_lo)
# rownames(distm_lo) = rownames(top_var)
# dist_lo = as.dist(distm_lo)

