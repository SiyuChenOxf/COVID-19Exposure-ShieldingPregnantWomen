#This code is to run 2-lambda model
rm(list = ls())
set.seed(100)

library(rstan)
library(ggpubr)
library(ggplot2)
library(bayesplot)

options (mc.cores = parallel::detectCores ())

## common values for plotting
ymax_kft <- 0.012
ymax_exposure <- 0.30
library(RColorBrewer)
colors_Dark<-brewer.pal(7,"Dark2")
colors_Spectral<-brewer.pal(7,"Spectral")
font_size = 14
font_size_title = 16
lwd = 0.5
pt_size = 0.4
right_margin=1

#preparing input data for Stan
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

y <- readRDS(paste(folder,"/data.rds", sep = ""))

start_week<-1   #2020-4-20
end_week<-36    #2020-12-28
t <- seq(start_week, end_week ,by=1) # No. of week that we want to make predictions
n_week<-length(t)                    
n_week2<-length(y[,1])
t2<-c(1:5,14:21,24:36)               #No. of week that raw data exist
t0 = 0                               #initial time point

data_pw <- list(t0 = t0, t = t,t2=t2,n_week=n_week,n_week2=n_week2,y=y)

niter <- 2000 #number of MCMC steps

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Models"
folder = paste(folder_strings, sep = "", collapse = "/")

fit_model <- stan(paste(folder,"/MfPW_2para.stan", sep = ""),
                      data = data_pw,
                      iter = niter,
                      chains = 4,
                      seed = 0,
                      control = list(adapt_delta = 0.99))

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

saveRDS(fit_model, file=paste(folder,"/fit_2para.rds", sep = ""))

##Rstan posterior estimations of parameters
tau <- rstan::extract(fit_model)$tau
sigma <- rstan::extract(fit_model)$sigma
beta <- rstan::extract(fit_model)$beta
y00<-rstan::extract(fit_model)$y00
k1<-rstan::extract(fit_model)$k1 
k2<-rstan::extract(fit_model)$k2
k01<-rstan::extract(fit_model)$k01
k11<-rstan::extract(fit_model)$k11
k10<-rstan::extract(fit_model)$k10

quantile(7/tau,c(0.05,0.25,0.5,0.75,0.95))
quantile(7/sigma,c(0.05,0.25,0.5,0.75,0.95))
quantile(7/beta,c(0.05,0.25,0.5,0.75,0.95))
quantile(y00,c(0.05,0.25,0.5,0.75,0.95))

mcmc_intervals(fit_model, pars = c("k10", "k11", "k01"))
mcmc_intervals(fit_model, pars = c("k1", "k2"))

# pairs(fit_model, pars = c("sigma", "tau", "beta", "y00","k1","k2"))
# traceplot(fit_model, pars = c("sigma", "tau", "beta", "y00","k1","k2"))
# traceplot(fit_model, pars = c("sigma", "tau", "beta", "y00"), inc_warmup = TRUE, nrow = 2)
# 
# #parameter posterior density
# stan_dens(fit_model, pars = c("sigma", "tau", "beta", "y00"), separate_chains = TRUE)
# stan_dens(fit_model, pars = c("sigma", "tau", "beta", "y00","k1","k2","k11","k10","k01"), separate_chains = TRUE)
# 
# mcmc_areas(fit_model,pars = c("k1", "k2"),prob = 0.95)     #95% credible interval 
# mcmc_areas(fit_model,pars = c("k01", "k11", "k10"),prob = 0.95)     #95% credible interval 

Result_median<-readRDS(file=paste(folder,"/Result_median_4.rds", sep = ""))
Result_upp<-readRDS(file=paste(folder,"/Result_upp_4.rds", sep = ""))
Result_low<-readRDS(file=paste(folder,"/Result_low_4.rds", sep = ""))

##Rstan posterior estimations of proportions of pregnant women in 4 compartments (merging x00 with z00 into one compartment)
res00<-res10<-res11<-res01<-matrix(0,nrow = 36,ncol = 5)
for (i in 1:36) {
  res00[i,]<-quantile(unlist(extract(fit_model, pars = paste0("out[", i,",","1", "]"), inc_warmup = FALSE)),c(0.05,0.25,0.5,0.75,0.95))
  res10[i,]<-quantile(unlist(extract(fit_model, pars = paste0("out[", i,",","2", "]"), inc_warmup = FALSE)),c(0.05,0.25,0.5,0.75,0.95))
  res11[i,]<-quantile(unlist(extract(fit_model, pars = paste0("out[", i,",","3", "]"), inc_warmup = FALSE)),c(0.025,0.25,0.5,0.75,0.975))
  res01[i,]<-quantile(unlist(extract(fit_model, pars = paste0("out[", i,",","4", "]"), inc_warmup = FALSE)),c(0.05,0.25,0.5,0.75,0.95))
}

##Rstan posterior estimations of proportions of pregnant women in 5 compartments 
z00<-y00<-y10<-y11<-y01<-matrix(0,nrow = 36,ncol = 5)
for (i in 1:36) {
  y00[i,]<-quantile(unlist(extract(fit_model, pars = paste0("y_hat[", i,",","1", "]"), inc_warmup = FALSE)),c(0.05,0.25,0.5,0.75,0.95))
  y10[i,]<-quantile(unlist(extract(fit_model, pars = paste0("y_hat[", i,",","2", "]"), inc_warmup = FALSE)),c(0.05,0.25,0.5,0.75,0.95))
  y11[i,]<-quantile(unlist(extract(fit_model, pars = paste0("y_hat[", i,",","3", "]"), inc_warmup = FALSE)),c(0.05,0.25,0.5,0.75,0.95))
  y01[i,]<-quantile(unlist(extract(fit_model, pars = paste0("y_hat[", i,",","4", "]"), inc_warmup = FALSE)),c(0.05,0.25,0.5,0.75,0.95))
  z00[i,]<-quantile(unlist(extract(fit_model, pars = paste0("y_hat[", i,",","5", "]"), inc_warmup = FALSE)),c(0.05,0.25,0.5,0.75,0.95))
}

date<-as.Date("2020-04-20")+c(0,7*seq(35))
data_x00 = data.frame(t=date, 
                      median = res00[,3], 
                      lower = res00[,1], 
                      upper = res00[,5],
                      lower1 = res00[,2], 
                      upper1 = res00[,4])


data1_x00 = data.frame( t=date,med= c(Result_median[1:5,1],rep(NA,8),Result_median[6:13,1],rep(NA,2),Result_median[14:26,1]),upp=c(Result_upp[1:5,1],rep(NA,8),Result_upp[6:13,1],rep(NA,2),Result_upp[14:26,1]),low=c(Result_low[1:5,1],rep(NA,8),Result_low[6:13,1],rep(NA,2),Result_low[14:26,1]))

x00p1<-ggplot(data_x00, aes(x=t, y = median)) +
  ylab("Proportion")+xlab("")+
  geom_line(size = 1) +
  ggtitle("")+
  scale_x_continuous(breaks = as.Date(date)[c(1,4,8,12,16,20,24,28,32,36)])+
  geom_pointrange(data=data1_x00, aes(x=t,y=med,ymin=low, ymax=upp), inherit.aes = FALSE, shape = 21, size=0.5, colour = "black", fill = colors_Dark[2])+
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "#D55E00", alpha=0.3, colour = NA)+
  geom_ribbon(aes(ymin=lower1, ymax=upper1), fill = "#D55E00", alpha=0.5, colour = NA)+
  # geom_ribbon(aes(ymin=lower1, ymax=upper1), alpha=0.5, colour = NA)+
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  theme(
    text = element_text(size=18 ),
    plot.title = element_text(face = "bold", size = 18 ,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=18 ),
    axis.title.x = element_text(size=18 ),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )


data_x10 = data.frame(t=date, 
                      median = res10[,3], 
                      lower = res10[,1], 
                      upper = res10[,5],
                      lower1 = res10[,2], 
                      upper1 = res10[,4])
data1_x10 = data.frame( t=date,med= c(Result_median[1:5,2],rep(NA,8),Result_median[6:13,2],rep(NA,2),Result_median[14:26,2]),upp=c(Result_upp[1:5,2],rep(NA,8),Result_upp[6:13,2],rep(NA,2),Result_upp[14:26,2]),low=c(Result_low[1:5,2],rep(NA,8),Result_low[6:13,2],rep(NA,2),Result_low[14:26,2]))

x10p1<-ggplot(data_x10, aes(x=t, y = median)) +
  ylab("Proportion")+
  scale_x_continuous(breaks = as.Date(date)[c(1,4,8,12,16,20,24,28,32,36)])+
  geom_line(size = 1) +
  ggtitle("")+xlab("")+
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "#D55E00", alpha=0.3, colour = NA)+
  geom_ribbon(aes(ymin=lower1, ymax=upper1), fill = "#D55E00", alpha=0.5, colour = NA)+
  geom_pointrange(data=data1_x10, aes(x=t,y=med,ymin=low, ymax=upp), inherit.aes = FALSE, shape = 21, size=0.5, colour = "black", fill = colors_Dark[2])+
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  theme(
    text = element_text(size=18),
    plot.title = element_text(face = "bold", size = 18,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=18),
    axis.title.x = element_text(size=18),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )


data_x11 = data.frame(t=date, 
                      median = res11[,3], 
                      lower = res11[,1], 
                      upper = res11[,5],
                      lower1 = res11[,2], 
                      upper1 = res11[,4])
data1_x11 = data.frame( t=date,med= c(Result_median[1:5,3],rep(NA,8),Result_median[6:13,3],rep(NA,2),Result_median[14:26,3]),upp=c(Result_upp[1:5,3],rep(NA,8),Result_upp[6:13,3],rep(NA,2),Result_upp[14:26,3]),low=c(Result_low[1:5,3],rep(NA,8),Result_low[6:13,3],rep(NA,2),Result_low[14:26,3]))
data_x11$lower[data_x11$lower<=0]<-0

x11p1<-ggplot(data_x11, aes(x=t, y = median))+
  ylab("Proportion")+
  scale_x_continuous(breaks = as.Date(date)[c(1,4,8,12,16,20,24,28,32,36)])+
  scale_y_continuous(breaks = c(0.00,0.05,0.10,0.15,0.20),limits = c(0,0.2))+
  xlab("")+
  geom_line(size = 1) +  
  ggtitle("")+
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "#D55E00", alpha=0.3, colour = NA)+
  geom_ribbon(aes(ymin=lower1, ymax=upper1), fill = "#D55E00", alpha=0.5, colour = NA)+
  geom_pointrange(data=data1_x11, aes(x=t,y=med,ymin=low, ymax=upp), inherit.aes = FALSE, shape = 21, size=0.5, colour = "black", fill = colors_Dark[2])+
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  theme(
    text = element_text(size=18),
    plot.title = element_text(face = "bold", size = 18,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=18),
    axis.title.x = element_text(size=18),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

data_x01 = data.frame(t=date, 
                      median = res01[,3], 
                      lower = res01[,1], 
                      upper = res01[,5],
                      lower1 = res01[,2], 
                      upper1 = res01[,4])
data1_x01 = data.frame( t=date,med= c(Result_median[1:5,4],rep(NA,8),Result_median[6:13,4],rep(NA,2),Result_median[14:26,4]),upp=c(Result_upp[1:5,4],rep(NA,8),Result_upp[6:13,4],rep(NA,2),Result_upp[14:26,4]),low=c(Result_low[1:5,4],rep(NA,8),Result_low[6:13,4],rep(NA,2),Result_low[14:26,4]))

x01p1<-ggplot(data_x01, aes(x=t, y = median)) +
  ylab("Proportion")+
  scale_x_continuous(breaks = as.Date(date)[c(1,4,8,12,16,20,24,28,32,36)])+
  scale_y_continuous(breaks = c(0.00,0.05,0.10,0.15,0.20),limits = c(0,0.2))+
  xlab("")+
  geom_line(size = 1) +  
  ggtitle("")+
  geom_ribbon(aes(ymin=lower, ymax=upper), fill = "#D55E00", alpha=0.3, colour = NA)+
  geom_ribbon(aes(ymin=lower1, ymax=upper1), fill = "#D55E00", alpha=0.5, colour = NA)+
  geom_pointrange(data=data1_x01, aes(x=t,y=med,ymin=low, ymax=upp), inherit.aes = FALSE, shape = 21, size=0.5, colour = "black", fill = colors_Dark[2])+
  # scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  theme(
    text = element_text(size=18),
    plot.title = element_text(face = "bold", size = 18,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=18),
    axis.title.x = element_text(size=18),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

ggarrange(x00p1, x10p1, x11p1,x01p1,
          labels = c("(A)", "(B)", "(C)","(D)"),
          ncol = 2, nrow = 2)

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

saveRDS(data_x01, file=paste(folder,"/x01_M2.rds", sep = ""))
saveRDS(data_x11, file=paste(folder,"/x11_M2.rds", sep = ""))
saveRDS(data_x10, file=paste(folder,"/x10_M2.rds", sep = ""))
saveRDS(data_x00, file=paste(folder,"/x00_M2.rds", sep = ""))


###expouser-seroprevalence
data_seroexpo = data.frame(output = c(rep("Exposure", 36), rep("Seroprevalence", 36)), 
                           t=c(date,date), 
                           median = c((1-y00)[,3], 
                                      (res11+res01)[,3]), 
                           lower1 = c((1-y00)[,5], 
                                      (res11+res01)[,1] ), 
                           upper1 = c((1-y00)[,1],
                                      (res11+res01)[,5] ),
                           lower2 = c((1-y00)[,4], 
                                      (res11+res01)[,2] ),
                           upper2 = c((1-y00)[,2], 
                                      (res11+res01)[,4] ))
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

saveRDS(data_seroexpo, file=paste(folder,"/seroexpo_M2.rds", sep = ""))

Result_median<-readRDS(paste(folder,"/Result_median_2.rds", sep = ""))
Result_upp<-readRDS(paste(folder,"/Result_upp_2.rds", sep = ""))
Result_low<-readRDS(paste(folder,"/Result_low_2.rds", sep = ""))

data2<-data.frame(date = data1_x11[1],
                  med = c(Result_median[1:5],rep(NA,8),Result_median[6:13],rep(NA,2),Result_median[14:26]),
                  upp = c(Result_upp[1:5],rep(NA,8),Result_upp[6:13],rep(NA,2),Result_upp[14:26]),
                  low = c(Result_low[1:5],rep(NA,8),Result_low[6:13],rep(NA,2),Result_low[14:26]))
data2<-na.omit(data2)
data_point = data.frame( t=data2$t, value= data2$med, upper= data2$upp,lower = data2$low )

p<-ggplot(data_seroexpo, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  xlab("")+ylab("")+
  geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.1, colour = NA)+
  geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.3, colour = NA)+
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6), limit = c(0, 0.6))+
  scale_x_continuous(breaks = as.Date(date)[c(1,4,8,12,16,20,24,28,32,36)])+
  geom_pointrange(data=data_point, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])


p+scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  ggtitle("Seroprevalence Vs exposure estimated by pregnant women data")+
  theme(
    text = element_text(size=18),
    plot.title = element_text(face = "bold", size = 18,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=18),
    axis.title.x = element_text(size=18),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )


