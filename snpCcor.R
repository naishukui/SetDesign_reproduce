# For continuous outcome and correlated genotypes, simulate the power for true and misspecified  models using SKAT,
# and also calculate the analytical power using  Lee's method

library(SKAT)
library(bindata)
library(expm)
output_path<-"/path_to_the_output_file/"
output_file_name<-"snpCcor_"

index <- as.numeric(commandArgs(TRUE)[1])

beta1list <- list(
  c(-0.5, -0.4, -0.3, -0.2, -0.1, -0.5, -0.4, -0.3, -0.2, -0.1,
    -0.5, -0.4, -0.3, -0.2, -0.1, -0.5, -0.4, -0.3, -0.2, -0.1,
    -0.5, -0.4, -0.3, -0.2, -0.1),
  c(-0.9, -0.75, -0.5, -0.25, -0.1, -0.9, -0.75, -0.5, -0.25, -0.1,
    -0.9, -0.75, -0.5, -0.25, -0.1, -0.9, -0.75, -0.5, -0.25, -0.1,
    -0.9, -0.75, -0.5, -0.25, -0.1),
  c(-1.5, -1.25, -1, -0.75, -0.5, -1.5, -1.25, -1, -0.75, -0.5,
    -1.5, -1.25, -1, -0.75, -0.5, -1.5, -1.25, -1, -0.75, -0.5,
    -1.5, -1.25, -1, -0.75, -0.5)
)

beta2list <- list(
  c(rep(0.5, 5), rep(0.4, 5),  rep(0.3, 5), rep(0.2, 5),  rep(0.1, 5)),
  c(rep(0.9, 5), rep(0.75, 5), rep(0.5, 5), rep(0.25, 5), rep(0.1, 5)),
  c(rep(1.5, 5), rep(1.25, 5), rep(1, 5), rep(0.75, 5), rep(0.5, 5)),
  c(rep(3, 5), rep(2.5, 5), rep(2, 5), rep(1.5, 5), rep(1, 5))
)

result<-powerC_meanCor(n = 2000, k = 50, maf_alt = 0.05,rho=0.15, runs = 1000, list1 = beta1list[[index]], list2 = beta2list[[index]])
save(result,file=paste0(output_path,output_file_name,index,".Rdata"))
