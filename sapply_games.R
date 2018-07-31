v = c(1:10)
v
dim(v) = c(2,5)
v

mat = matrix(v, ncol=2, nrow=5, byrow = T)
mat

sapply(mat, function(x) x + 2)
lapply(mat, function(x) x + 2)
apply(mat, 2, function(x) x + 2)
mat
