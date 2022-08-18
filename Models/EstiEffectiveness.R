#Proportion of reduction  
#To generate Table S5
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

fit_1para <- readRDS(paste(folder,"/fit_1para.rds", sep = ""))
y00<-rstan::extract(fit_1para)$y00
Re<-seq(10000)
for (i in 1:10000) {
  Re[i]<-(rnorm(1,mean = 0.29,sd = 0.02)-sample((1-y00),1))/rnorm(1,mean = 0.29,sd = 0.02)
}
quantile(Re,c(0.05,0.5,0.95))

#Proportion of reduction from Model 2
fit_2para <- readRDS(paste(folder,"/fit_2para.rds", sep = ""))
y00<-rstan::extract(fit_2para)$y00
Re<-seq(10000)
for (i in 1:10000) {
  Re[i]<-(rnorm(1,mean = 0.29,sd = 0.02)-sample((1-y00),1))/rnorm(1,mean = 0.29,sd = 0.02)
}
quantile(Re,c(0.05,0.5,0.95))

#Proportion of reduction from Model 3
fit_3para <- readRDS(paste(folder,"/fit_3para.rds", sep = ""))
y00<-rstan::extract(fit_3para)$y00
Re<-seq(10000)
for (i in 1:10000) {
  Re[i]<-(rnorm(1,mean = 0.29,sd = 0.02)-sample((1-y00),1))/rnorm(1,mean = 0.29,sd = 0.02)
}
quantile(Re,c(0.05,0.5,0.95))

#Proportion of reduction from Model 4
fit_4para <- readRDS(paste(folder,"/fit_4para.rds", sep = ""))
y00<-rstan::extract(fit_4para)$y00
Re<-seq(10000)
for (i in 1:10000) {
  Re[i]<-(rnorm(1,mean = 0.29,sd = 0.02)-sample((1-y00),1))/rnorm(1,mean = 0.29,sd = 0.02)
}
quantile(Re,c(0.05,0.5,0.95))