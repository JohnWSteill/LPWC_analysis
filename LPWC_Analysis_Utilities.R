#install.packages("LPWC", repos = "http://cran.us.r-project.org")
library(LPWC)

dbg_fn <- function(){
    return (list(solution=1, dist=2, clust=3))
}

get_high_var_genes <- function(dataFile, sd_cv_thresh=500){
    # The true coefficient of variance is unitless, 
    # I squared the sd to get larger counts. 
    my_data = read.table(dataFile, header=TRUE, row.names=1, sep = '\t',  quote = "", comment.char = "")
    my_data = my_data[,-ncol(my_data)] # Drop description
    my_data$sd_coef_var = apply(my_data, 1, function(x) sd(x) * sd(x)/ mean(x))
    my_data$sd_coef_var[is.na(my_data$sd_coef_var)]<-0
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

