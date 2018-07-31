library(LPWC)
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
# '''

# #install.packages("LPWC", repos = "http://cran.us.r-project.org")
# 
# load(file='wkspace.RData')
# #library(ggplot2)
# # dataFile = '/isiseqruns/.GUP_HOME/RUNS/C7DCVACXX/Collate/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9/genes.no_mt.ec.tab'
# dataFile = 'Data/Sub_391_2_pool_22adb_hg19_SOX2_plus__foregut_endoderm__32a860ee85db9ac9/genes.no_mt.ec.tab'
# my_data = read.table(dataFile, header=TRUE, row.names=1, sep = '\t')
# my_data = my_data[,-ncol(my_data)] # Drop description
# # The true coefficient of variance is unitless, I squared the sd to get larger counts. 
# my_data$coefVar = apply(my_data, 1, function(x) sd(x) * sd(x)/ mean(x))
# my_data$coefVar[is.na(my_data$coefVar)]<-0
# top_var = subset(my_data[my_data$coefVar>500,],select=-c(coefVar))
# timepoints = 0:96
# solution_lo = LPWC::corr.bestlag(top_var, timepoints = timepoints, max.lag = 20, penalty = "low", iter = 10)
# dist_lo <- 1 - solution_lo$corr
# distm_lo = as.matrix(dist_lo)
# rownames(distm_lo) = rownames(top_var)
# dist_lo = as.dist(distm_lo)
# clust_lo <- hclust(dist_lo)
# groups_lo = cutree(clust_lo, k=10)
# save.image(file='wkspace_lo.RData')

