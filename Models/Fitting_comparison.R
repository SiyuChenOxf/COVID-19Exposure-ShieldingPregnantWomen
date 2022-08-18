#This code is to plot 4 & 2-compartment model fitting comparison 
#For Fig 2 and Fig 3
rm(list = ls())
set.seed(100)

library(ggplot2)
library(rstan)
library(deSolve)
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

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

Result_median<-readRDS(file=paste(folder,"/Result_median_4.rds", sep = ""))
Result_upp<-readRDS(file=paste(folder,"/Result_upp_4.rds", sep = ""))
Result_low<-readRDS(file=paste(folder,"/Result_low_4.rds", sep = ""))


#MODEL PREDICTION MODEL 1
x00_hat_1<-readRDS(file=paste(folder,"/x00_M1.rds", sep = ""))
x10_hat_1<-readRDS(file=paste(folder,"/x10_M1.rds", sep = ""))
x11_hat_1<-readRDS(file=paste(folder,"/x11_M1.rds", sep = ""))
x01_hat_1<-readRDS(file=paste(folder,"/x01_M1.rds", sep = ""))

#MODEL 2
x00_hat_2<-readRDS(file=paste(folder,"/x00_M2.rds", sep = ""))
x10_hat_2<-readRDS(file=paste(folder,"/x10_M2.rds", sep = ""))
x11_hat_2<-readRDS(file=paste(folder,"/x11_M2.rds", sep = ""))
x01_hat_2<-readRDS(file=paste(folder,"/x01_M2.rds", sep = ""))


###MODEL 3
x00_hat_3<-readRDS(file=paste(folder,"/x00_M3.rds", sep = ""))
x10_hat_3<-readRDS(file=paste(folder,"/x10_M3.rds", sep = ""))
x11_hat_3<-readRDS(file=paste(folder,"/x11_M3.rds", sep = ""))
x01_hat_3<-readRDS(file=paste(folder,"/x01_M3.rds", sep = ""))


#MODEL 4
x00_hat_4<-readRDS(file=paste(folder,"/x00_M4.rds", sep = ""))
x10_hat_4<-readRDS(file=paste(folder,"/x10_M4.rds", sep = ""))
x11_hat_4<-readRDS(file=paste(folder,"/x11_M4.rds", sep = ""))
x01_hat_4<-readRDS(file=paste(folder,"/x01_M4.rds", sep = ""))

date<-as.Date("2020-04-20")+c(0,7*seq(35))

data_x00 = data.frame(output=c(rep("Model 1",length(date)),rep("Model 2",length(date)),rep("Model 3",length(date)),rep("Model 4",length(date))),
                      t=rep(date,4), 
                      median=c(
                        median_1 = x00_hat_1$median, 
                        median_2 = x00_hat_2$median, 
                        median_3 = x00_hat_3$median, 
                        median_4 = x00_hat_4$median), 
                      
                      low_25=c(
                        lower_1 = x00_hat_1$lower1, 
                        lower_2 = x00_hat_2$lower1, 
                        lower_3 = x00_hat_3$lower1, 
                        lower_4 = x00_hat_4$lower1), 
                      
                      upp_75=c(
                        upper_1 = x00_hat_1$upper1, 
                        upper_2 = x00_hat_2$upper1, 
                        upper_3 = x00_hat_3$upper1, 
                        upper_4 = x00_hat_4$upper1),
                      
                      low05=c(
                        lower1_1 = x00_hat_1$lower, 
                        lower2_1 = x00_hat_2$lower, 
                        lower3_1 = x00_hat_3$lower, 
                        lower4_1 = x00_hat_4$lower), 
                      
                      upp95=c(
                        upper1_1 = x00_hat_1$upper, 
                        upper2_1 = x00_hat_2$upper, 
                        upper3_1 = x00_hat_3$upper, 
                        upper4_1 = x00_hat_4$upper))

data1_x00 = data.frame( t=date,med= c(Result_median[1:5,1],rep(NA,8),Result_median[6:13,1],rep(NA,2),Result_median[14:26,1]),upp=c(Result_upp[1:5,1],rep(NA,8),Result_upp[6:13,1],rep(NA,2),Result_upp[14:26,1]),low=c(Result_low[1:5,1],rep(NA,8),Result_low[6:13,1],rep(NA,2),Result_low[14:26,1]))

x00p1<-ggplot(data_x00, aes(x=t, y = median, group = output, colour = output)) +
  ylab("Proportion of population")+xlab("")+
  geom_line(size = 1) +
  ggtitle("")+
  scale_x_date(breaks = c(as.Date(c("2020-04-01","2020-06-01","2020-08-01","2020-10-01","2020-12-01"))),
               labels = c("Apr","Jun","Aug","Oct","Dec"),
               limits = as.Date(c("2020-04-01","2020-12-31")))+
  # scale_y_continuous(breaks = c(0.00,0.05,0.10,0.15,0.20),limits = c(0,0.2))+
   geom_pointrange(data=data1_x00, aes(x=t,y=med,ymin=low, ymax=upp), inherit.aes = FALSE, shape = 21, size=0.5, colour = "black", fill = colors_Dark[2])+
  geom_ribbon(aes(ymin=low05, ymax=upp95, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower1, ymax=upper1), alpha=0.5, colour = NA)+
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  theme(
    text = element_text(size=20 ),
    plot.title = element_text(face = "bold", size = 20 ,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=20 ),
    axis.title.x = element_text(size=20 ),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

data_x10 = data.frame(output=c(rep("Model 1",length(date)),rep("Model 2",length(date)),rep("Model 3",length(date)),rep("Model 4",length(date))),
                      t=rep(date,4), 
                      median=c(
                        median_1 = x10_hat_1$median, 
                        median_2 = x10_hat_2$median, 
                        median_3 = x10_hat_3$median, 
                        median_4 = x10_hat_4$median), 
                      
                      low_25=c(
                        lower_1 = x10_hat_1$lower1, 
                        lower_2 = x10_hat_2$lower1, 
                        lower_3 = x10_hat_3$lower1, 
                        lower_4 = x10_hat_4$lower1), 
                      
                      upp_75=c(
                        upper_1 = x10_hat_1$upper1, 
                        upper_2 = x10_hat_2$upper1, 
                        upper_3 = x10_hat_3$upper1, 
                        upper_4 = x10_hat_4$upper1),
                      
                      low05=c(
                        lower1_1 = x10_hat_1$lower, 
                        lower2_1 = x10_hat_2$lower, 
                        lower3_1 = x10_hat_3$lower, 
                        lower4_1 = x10_hat_4$lower), 
                      
                      upp95=c(
                        upper1_1 = x10_hat_1$upper, 
                        upper2_1 = x10_hat_2$upper, 
                        upper3_1 = x10_hat_3$upper, 
                        upper4_1 = x10_hat_4$upper))
data1_x10 = data.frame( t=date,med= c(Result_median[1:5,2],rep(NA,8),Result_median[6:13,2],rep(NA,2),Result_median[14:26,2]),upp=c(Result_upp[1:5,2],rep(NA,8),Result_upp[6:13,2],rep(NA,2),Result_upp[14:26,2]),low=c(Result_low[1:5,2],rep(NA,8),Result_low[6:13,2],rep(NA,2),Result_low[14:26,2]))

x10p1<-ggplot(data_x10, aes(x=t, y = median, group = output, colour = output)) +
  ylab("Proportion of population")+
  scale_x_date(breaks = c(as.Date(c("2020-04-01","2020-06-01","2020-08-01","2020-10-01","2020-12-01"))),
               labels = c("Apr","Jun","Aug","Oct","Dec"),
               limits = as.Date(c("2020-04-01","2020-12-31")))+
  # scale_y_continuous(breaks = c(0.00,0.05,0.10,0.15,0.20),limits = c(0,0.2))+
  geom_line(size = 1) +
  ggtitle("")+xlab("")+
  geom_ribbon(aes(ymin=low05, ymax=upp95, fill = output), alpha=0.2, colour = NA)+
  # geom_point(data=data1_x10, aes(x=t,y=value))+
  geom_pointrange(data=data1_x10, aes(x=t,y=med,ymin=low, ymax=upp), inherit.aes = FALSE, shape = 21, size=0.5, colour = "black", fill = colors_Dark[2])+
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  theme(
    text = element_text(size=20),
    plot.title = element_text(face = "bold", size = 20,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(0.2,1),
    legend.title = element_blank(),
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

data_x11 = data.frame(output=c(rep("Model 1",length(date)),rep("Model 2",length(date)),rep("Model 3",length(date)),rep("Model 4",length(date))),
                      t=rep(date,4), 
                      median=c(
                        median_1 = x11_hat_1$median, 
                        median_2 = x11_hat_2$median, 
                        median_3 = x11_hat_3$median, 
                        median_4 = x11_hat_4$median), 
                      
                      low_25=c(
                        lower_1 = x11_hat_1$lower1, 
                        lower_2 = x11_hat_2$lower1, 
                        lower_3 = x11_hat_3$lower1, 
                        lower_4 = x11_hat_4$lower1), 
                      
                      upp_75=c(
                        upper_1 = x11_hat_1$upper1, 
                        upper_2 = x11_hat_2$upper1, 
                        upper_3 = x11_hat_3$upper1, 
                        upper_4 = x11_hat_4$upper1),
                      
                      low05=c(
                        lower1_1 = x11_hat_1$lower, 
                        lower2_1 = x11_hat_2$lower, 
                        lower3_1 = x11_hat_3$lower, 
                        lower4_1 = x11_hat_4$lower), 
                      
                      upp95=c(
                        upper1_1 = x11_hat_1$upper, 
                        upper2_1 = x11_hat_2$upper, 
                        upper3_1 = x11_hat_3$upper, 
                        upper4_1 = x11_hat_4$upper))
data1_x11 = data.frame( t=date,med= c(Result_median[1:5,3],rep(NA,8),Result_median[6:13,3],rep(NA,2),Result_median[14:26,3]),upp=c(Result_upp[1:5,3],rep(NA,8),Result_upp[6:13,3],rep(NA,2),Result_upp[14:26,3]),low=c(Result_low[1:5,3],rep(NA,8),Result_low[6:13,3],rep(NA,2),Result_low[14:26,3]))

x11p1<-ggplot(data_x11, aes(x=t, y = median, group = output, colour = output))+
  ylab("Proportion of population")+
  scale_x_date(breaks = c(as.Date(c("2020-04-01","2020-06-01","2020-08-01","2020-10-01","2020-12-01"))),
               labels = c("Apr","Jun","Aug","Oct","Dec"),
               limits = as.Date(c("2020-04-01","2020-12-31")))+
  scale_y_continuous(breaks = c(0.00,0.05,0.10,0.15,0.20),limits = c(0,0.2))+
  xlab("")+
  geom_line(size = 1) +  
  ggtitle(" ")+
  geom_ribbon(aes(ymin=low05, ymax=upp95, fill = output), alpha=0.2, colour = NA)+
  geom_pointrange(data=data1_x11, aes(x=t,y=med,ymin=low, ymax=upp), inherit.aes = FALSE, shape = 21, size=0.5, colour = "black", fill = colors_Dark[2])+
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  theme(
    text = element_text(size=20),
    plot.title = element_text(face = "bold", size = 20,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

data_x01 = data.frame(output=c(rep("Model 1",length(date)),rep("Model 2",length(date)),rep("Model 3",length(date)),rep("Model 4",length(date))),
                      t=rep(date,4), 
                      median=c(
                        median_1 = x01_hat_1$median, 
                        median_2 = x01_hat_2$median, 
                        median_3 = x01_hat_3$median, 
                        median_4 = x01_hat_4$median), 
                      
                      low_25=c(
                        lower_1 = x01_hat_1$lower1, 
                        lower_2 = x01_hat_2$lower1, 
                        lower_3 = x01_hat_3$lower1, 
                        lower_4 = x01_hat_4$lower1), 
                      
                      upp_75=c(
                        upper_1 = x01_hat_1$upper1, 
                        upper_2 = x01_hat_2$upper1, 
                        upper_3 = x01_hat_3$upper1, 
                        upper_4 = x01_hat_4$upper1),
                      
                      low05=c(
                        lower1_1 = x01_hat_1$lower, 
                        lower2_1 = x01_hat_2$lower, 
                        lower3_1 = x01_hat_3$lower, 
                        lower4_1 = x01_hat_4$lower), 
                      
                      upp95=c(
                        upper1_1 = x01_hat_1$upper, 
                        upper2_1 = x01_hat_2$upper, 
                        upper3_1 = x01_hat_3$upper, 
                        upper4_1 = x01_hat_4$upper))
data1_x01 = data.frame( t=date,med= c(Result_median[1:5,4],rep(NA,8),Result_median[6:13,4],rep(NA,2),Result_median[14:26,4]),upp=c(Result_upp[1:5,4],rep(NA,8),Result_upp[6:13,4],rep(NA,2),Result_upp[14:26,4]),low=c(Result_low[1:5,4],rep(NA,8),Result_low[6:13,4],rep(NA,2),Result_low[14:26,4]))

x01p1<-ggplot(data_x01, aes(x=t, y = median, group = output, colour = output)) +
  ylab("Proportion of population")+
  scale_x_date(breaks = c(as.Date(c("2020-04-01","2020-06-01","2020-08-01","2020-10-01","2020-12-01"))),
               labels = c("Apr","Jun","Aug","Oct","Dec"),
               limits = as.Date(c("2020-04-01","2020-12-31")))+
  scale_y_continuous(breaks = c(0.00,0.05,0.10,0.15,0.20,0.25),limits = c(0,0.25))+
  xlab("")+
  geom_line(size = 1) +  
  ggtitle("")+
  geom_ribbon(aes(ymin=low05, ymax=upp95, fill = output), alpha=0.2, colour = NA)+
  geom_pointrange(data=data1_x01, aes(x=t,y=med,ymin=low, ymax=upp), inherit.aes = FALSE, shape = 21, size=0.5, colour = "black", fill = colors_Dark[2])+
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +
  theme(
    text = element_text(size=20),
    plot.title = element_text(face = "bold", size = 20,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

ggarrange(x00p1, x10p1, x11p1,x01p1,
          labels = c("(A)", "(B)", "(C)","(D)"),
          ncol = 2, nrow = 2)

#MODEL 1
seroexpo_M1<-readRDS(file=paste(folder,"/seroexpo_M1.rds", sep = ""))
#MODEL 2
seroexpo_M2<-readRDS(file=paste(folder,"/seroexpo_M2.rds", sep = ""))
###MODEL 3
seroexpo_M3<-readRDS(file=paste(folder,"/seroexpo_M3.rds", sep = ""))
#MODEL 4
seroexpo_M4<-readRDS(file=paste(folder,"/seroexpo_M4.rds", sep = ""))
 

###To plot exposure-seroprevalence among different models 
data_seroexpo = data.frame(output = c(rep("Exposure Model 1",length(date)), rep("Seroprevalence Model 1", length(date)),
                                      rep("Exposure Model 2",length(date)), rep("Seroprevalence Model 2", length(date)),
                                      rep("Exposure Model 3",length(date)), rep("Seroprevalence Model 3", length(date)),
                                      rep("Exposure Model 4",length(date)), rep("Seroprevalence Model 4", length(date))), 
                           t=c(date,date), 
                           median = c(seroexpo_M1[seroexpo_M1$output=="Exposure","median"], 
                                      seroexpo_M1[seroexpo_M1$output=="Seroprevalence","median"],
                                      
                                      seroexpo_M2[seroexpo_M2$output=="Exposure","median"], 
                                      seroexpo_M2[seroexpo_M2$output=="Seroprevalence","median"],
                                      
                                      seroexpo_M3[seroexpo_M3$output=="Exposure","median"], 
                                      seroexpo_M3[seroexpo_M3$output=="Seroprevalence","median"],
                                      
                                      seroexpo_M4[seroexpo_M4$output=="Exposure","median"], 
                                      seroexpo_M4[seroexpo_M4$output=="Seroprevalence","median"]), 
                           
                           lower1 = c(seroexpo_M1[seroexpo_M1$output=="Exposure","lower1"], 
                                      seroexpo_M1[seroexpo_M1$output=="Seroprevalence","lower1"],
                                      
                                      seroexpo_M2[seroexpo_M2$output=="Exposure","lower1"], 
                                      seroexpo_M2[seroexpo_M2$output=="Seroprevalence","lower1"],
                                      
                                      seroexpo_M3[seroexpo_M3$output=="Exposure","lower1"], 
                                      seroexpo_M3[seroexpo_M3$output=="Seroprevalence","lower1"],
                                      
                                      seroexpo_M4[seroexpo_M4$output=="Exposure","lower1"], 
                                      seroexpo_M4[seroexpo_M4$output=="Seroprevalence","lower1"]),
                           
                           upper1 = c(seroexpo_M1[seroexpo_M1$output=="Exposure","upper1"], 
                                      seroexpo_M1[seroexpo_M1$output=="Seroprevalence","upper1"],
                                      
                                      seroexpo_M2[seroexpo_M2$output=="Exposure","upper1"], 
                                      seroexpo_M2[seroexpo_M2$output=="Seroprevalence","upper1"],
                                      
                                      seroexpo_M3[seroexpo_M3$output=="Exposure","upper1"], 
                                      seroexpo_M3[seroexpo_M3$output=="Seroprevalence","upper1"],
                                      
                                      seroexpo_M4[seroexpo_M4$output=="Exposure","median"], 
                                      seroexpo_M4[seroexpo_M4$output=="Seroprevalence","median"]),
                           
                           lower2 = c(seroexpo_M1[seroexpo_M1$output=="Exposure","lower2"], 
                                      seroexpo_M1[seroexpo_M1$output=="Seroprevalence","lower2"],
                                      
                                      seroexpo_M2[seroexpo_M2$output=="Exposure","lower2"], 
                                      seroexpo_M2[seroexpo_M2$output=="Seroprevalence","lower2"],
                                      
                                      seroexpo_M3[seroexpo_M3$output=="Exposure","lower2"], 
                                      seroexpo_M3[seroexpo_M3$output=="Seroprevalence","lower2"],
                                      
                                      seroexpo_M4[seroexpo_M4$output=="Exposure","median"], 
                                      seroexpo_M4[seroexpo_M4$output=="Seroprevalence","median"]),
                           
                           upper2 = c(seroexpo_M1[seroexpo_M1$output=="Exposure","upper2"], 
                                      seroexpo_M1[seroexpo_M1$output=="Seroprevalence","upper2"],
                                      
                                      seroexpo_M2[seroexpo_M2$output=="Exposure","upper2"], 
                                      seroexpo_M2[seroexpo_M2$output=="Seroprevalence","upper2"],
                                      
                                      seroexpo_M3[seroexpo_M3$output=="Exposure","upper2"], 
                                      seroexpo_M3[seroexpo_M3$output=="Seroprevalence","upper2"],
                                      
                                      seroexpo_M4[seroexpo_M4$output=="Exposure","upper2"], 
                                      seroexpo_M4[seroexpo_M4$output=="Seroprevalence","upper2"]))

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

Result_median <- readRDS(file=paste(folder,"/Result_median_2.rds", sep = ""))
Result_upp <- readRDS(file=paste(folder,"/Result_upp_2.rds", sep = ""))
Result_low <- readRDS(file=paste(folder,"/Result_low_2.rds", sep = ""))

data_point = data.frame( t=date,value= c(Result_median[1:5],rep(NA,8),Result_median[6:13],rep(NA,2),Result_median[14:26]),upper=c(Result_upp[1:5],rep(NA,8),Result_upp[6:13],rep(NA,2),Result_upp[14:26]),lower=c(Result_low[1:5],rep(NA,8),Result_low[6:13],rep(NA,2),Result_low[14:26]))

p<-ggplot(data_seroexpo, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = lwd) +  xlab("")+ylab("")+
  # geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.2, colour = NA)+
  geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.2, colour = NA)+
  scale_x_continuous(breaks = as.Date(date)[c(1,4,8,12,16,20,24,28,32,36)])+
  scale_x_date(breaks = c(as.Date(c("2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01"))),
               labels = c("Apr 2020","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan 2021"),
               limits = as.Date(c("2020-04-01","2021-01-01")))+
  geom_pointrange(data=data_point, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=pt_size, colour = "black", fill = colors_Dark[2])+
  annotate("text", x = as.Date("2021-01-01"), y = 0.40, label = "Exposure",size = 5)+
  annotate("text", x = as.Date("2020-12-25"), y = 0.06, label = "Seroprevalence",size = 5)

p5<-p+scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  theme_minimal() +xlab("")+ylab("Proportion of population")+
  ggtitle("Trajectory of exposure and seroprevalence estimated in pregnant women")+
  theme(
    text = element_text(size=20),
    plot.title = element_text(face = "bold", size = 20,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(0.20,1.0),
    legend.title = element_blank(),
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    # axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )

###To plot exposure-seroprevalence for general population 
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

data<-read.csv(file=paste(folder,"/MultiSeroModelResults.csv", sep = ""))

data1 = data.frame(output = c(rep("Exposure", length(data$Date)), rep("Seroprevalence", length(data$Date))), 
                   t=c(as.Date(as.Date("2020-01-02")+seq(from=0,to=364,by=1)),as.Date(as.Date("2020-01-02")+seq(from=0,to=364,by=1))), 
                   median = c(gene_expo_inf_med=data$Expo_NYC.median,gene_sero_inf_med=data$Sero_NYC.median),
                   lower  = c(gene_expo_inf_med=data$Expo_NYC.95.lower,gene_sero_inf_med=data$Sero_NYC.95.lower), 
                   upper  = c(gene_expo_inf_med=data$Expo_NYC.95.upper,gene_sero_inf_med=data$Sero_NYC.95.upper))
data2<-data.frame(t=as.Date(as.Date("2020-01-02")+seq(from=0,to=364,by=1)),
                  value=data$adj_sero_NYC_median,
                  upper=data$adj_sero_NYC_upper,
                  lower=data$adj_sero_NYC_lower)
data3<-data.frame(t=as.Date(as.Date("2020-01-02")+seq(from=0,to=364,by=1)),
                  value=data$preg.obs,
                  upper=data$preg.obs,
                  lower=data$preg.obs)


p6<-ggplot(data1, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("")+xlab("")+ylab("Proportion")+
  # geom_point(data=data3,aes(x=t, y = value))+
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limit = c(0, 0.6))+
  scale_x_continuous(breaks = c(as.Date(c("2020-01-01","2020-01-27","2020-02-17","2020-03-09","2020-03-30")),as.Date(date)[c(seq(1,length(date),by=3))],as.Date("2020-12-31")))+
  geom_pointrange(data=data2, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 21, size=0.8, colour = "black", fill = colors_Dark[2])+
  # geom_pointrange(data=data3, aes(x=t,y=value,ymin=lower, ymax=upper), inherit.aes = FALSE, shape = 19, size=0.8, colour = "black", fill = colors_Dark[1])+
  theme_minimal() +
  theme(
    text = element_text(size=20),
    plot.title = element_text(face = "bold", size = 20,hjust = 0.5),
    legend.background = element_rect(fill = "white", size = 1.25, colour = "white"),
    legend.justification = c(0, 1),
    legend.position = c(0.05,1.0),
    legend.title = element_blank(),
    axis.title.y = element_text(size=20),
    axis.title.x = element_text(size=20),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  )
p6
ggarrange(p5, p6,
          labels = c("A", "B"),
          ncol = 1, nrow = 2)

