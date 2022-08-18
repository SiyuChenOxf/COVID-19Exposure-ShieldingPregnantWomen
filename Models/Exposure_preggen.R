#This code is to plot Fig 4
rm(list = ls())
set.seed(100)

library(ggplot2)
library(rstan)
library(deSolve)
library(ggpubr)
library(ggplot2)
library(bayesplot)
library(dplyr)

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

data<-read.csv(paste(folder,"/MultiSeroModelResults.csv", sep = ""))

## common values
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


data1 = data.frame(output = c(rep("Exposure estimated by general population", length(data$Date)), rep("Seroprevalence estimated by general population", length(data$Date))), 
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


ggplot(data1, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("")+xlab("")+ylab("Proportion")+
  # geom_point(data=data3,aes(x=t, y = value))+
  geom_ribbon(aes(ymin=lower, ymax=upper, fill = output), alpha=0.2, colour = NA)+
  # geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.5, colour = NA)+
  # scale_y_continuous(breaks = c(0,5,10,15,20,25), limit = c(0, ymax_exposure))+
  scale_x_date(breaks = as.Date(c("2020-01-01","2020-02-01","2020-03-01","2020-04-01","2020-05-01","2020-06-01","2020-07-01","2020-08-01","2020-09-01","2020-10-01","2020-11-01","2020-12-01","2021-01-01")), labels=c("Jan 2020","Feb", "Mar", "Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Jan 2021"), limit = as.Date(c("2020-01-02","2020-12-31")))+
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
#---------u=4--------#
seroexpo_M1 <- readRDS(paste(folder,"/seroexpo_M4.rds", sep = ""))
seroexpo_M1<-seroexpo_M1[(seroexpo_M1$output=="Exposure"),-1]
seroexpo_M1$lower<-seroexpo_M1$lower1
seroexpo_M1$upper<-seroexpo_M1$upper1
data_preg<-seroexpo_M1
data_preg$Date<-as.Date(data_preg$t)

data_tot<-read.csv(paste(folder,"/MultiSeroModelResults.csv", sep = ""))
data_tot$Date<-as.Date(as.Date("2020-01-02")+seq(from=0,to=length(data_tot$Date)-1))
data<-merge(data_tot,data_preg,by.x =  "Date",by.y = "Date",all=TRUE)
data$median<-as.numeric(data$median)

x0<-data$Date[1:(which(data$median>0)[1]-1)]
x1<-data$Expo_NYC.median[1:(which(data$median>0)[1]-1)]
x2<-data$Expo_NYC.95.lower[1:(which(data$median>0)[1]-1)]
x3<-data$Expo_NYC.95.upper[1:(which(data$median>0)[1]-1)]

dataNEW<-data[which(data$median>0),]

xx0<-c(x0,dataNEW$Date)
xx1<-c(x1,dataNEW$Expo_NYC.median)
xx2<-c(x2,dataNEW$Expo_NYC.95.lower)
xx3<-c(x3,dataNEW$Expo_NYC.95.upper)

yy1_4<-c(rep(NA,length(x1)),dataNEW$median)
yy2_4<-c(rep(NA,length(x2)),dataNEW$lower)
yy3_4<-c(rep(NA,length(x3)),dataNEW$upper)
yy4_4<-c(rep(NA,length(x2)),dataNEW$lower1)
yy5_4<-c(rep(NA,length(x3)),dataNEW$upper1)
yy6_4<-c(rep(NA,length(x2)),dataNEW$lower2)
yy7_4<-c(rep(NA,length(x3)),dataNEW$upper2)

#---------u=3--------#
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

seroexpo_M1 <- readRDS(paste(folder,"/seroexpo_M3.rds", sep = ""))
seroexpo_M1<-seroexpo_M1[(seroexpo_M1$output=="Exposure"),-1]
seroexpo_M1$lower<-seroexpo_M1$lower1
seroexpo_M1$upper<-seroexpo_M1$upper1
data_preg<-seroexpo_M1
data_preg$Date<-as.Date(data_preg$t)

data_tot<-read.csv(paste(folder,"/MultiSeroModelResults.csv", sep = ""))
data_tot$Date<-as.Date(as.Date("2020-01-02")+seq(from=0,to=length(data_tot$Date)-1))
data<-merge(data_tot,data_preg,by.x =  "Date",by.y = "Date",all=TRUE)
data$median<-as.numeric(data$median)

x0<-data$Date[1:(which(data$median>0)[1]-1)]
x1<-data$Expo_NYC.median[1:(which(data$median>0)[1]-1)]
x2<-data$Expo_NYC.95.lower[1:(which(data$median>0)[1]-1)]
x3<-data$Expo_NYC.95.upper[1:(which(data$median>0)[1]-1)]

dataNEW<-data[which(data$median>0),]

xx0<-c(x0,dataNEW$Date)
xx1<-c(x1,dataNEW$Expo_NYC.median)
xx2<-c(x2,dataNEW$Expo_NYC.95.lower)
xx3<-c(x3,dataNEW$Expo_NYC.95.upper)

yy1_3<-c(rep(NA,length(x1)),dataNEW$median)
yy2_3<-c(rep(NA,length(x2)),dataNEW$lower)
yy3_3<-c(rep(NA,length(x3)),dataNEW$upper)
yy4_3<-c(rep(NA,length(x2)),dataNEW$lower1)
yy5_3<-c(rep(NA,length(x3)),dataNEW$upper1)
yy6_3<-c(rep(NA,length(x2)),dataNEW$lower2)
yy7_3<-c(rep(NA,length(x3)),dataNEW$upper2)

#---------u=2--------#
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

seroexpo_M1 <- readRDS(paste(folder,"/seroexpo_M2.rds", sep = ""))
seroexpo_M1<-seroexpo_M1[(seroexpo_M1$output=="Exposure"),-1]
seroexpo_M1$lower<-seroexpo_M1$lower1
seroexpo_M1$upper<-seroexpo_M1$upper1
data_preg<-seroexpo_M1
data_preg$Date<-as.Date(data_preg$t)

data_tot<-read.csv(paste(folder,"/MultiSeroModelResults.csv", sep = ""))
data_tot$Date<-as.Date(as.Date("2020-01-02")+seq(from=0,to=length(data_tot$Date)-1))
data<-merge(data_tot,data_preg,by.x =  "Date",by.y = "Date",all=TRUE)
data$median<-as.numeric(data$median)

x0<-data$Date[1:(which(data$median>0)[1]-1)]
x1<-data$Expo_NYC.median[1:(which(data$median>0)[1]-1)]
x2<-data$Expo_NYC.95.lower[1:(which(data$median>0)[1]-1)]
x3<-data$Expo_NYC.95.upper[1:(which(data$median>0)[1]-1)]

dataNEW<-data[which(data$median>0),]

xx0<-c(x0,dataNEW$Date)
xx1<-c(x1,dataNEW$Expo_NYC.median)
xx2<-c(x2,dataNEW$Expo_NYC.95.lower)
xx3<-c(x3,dataNEW$Expo_NYC.95.upper)

yy1_2<-c(rep(NA,length(x1)),dataNEW$median)
yy2_2<-c(rep(NA,length(x2)),dataNEW$lower)
yy3_2<-c(rep(NA,length(x3)),dataNEW$upper)
yy4_2<-c(rep(NA,length(x2)),dataNEW$lower1)
yy5_2<-c(rep(NA,length(x3)),dataNEW$upper1)
yy6_2<-c(rep(NA,length(x2)),dataNEW$lower2)
yy7_2<-c(rep(NA,length(x3)),dataNEW$upper2)

#---------u=1--------#
folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")

seroexpo_M1 <- readRDS(paste(folder,"/seroexpo_M1.rds", sep = ""))
seroexpo_M1<-seroexpo_M1[(seroexpo_M1$output=="Exposure"),-1]
seroexpo_M1$lower<-seroexpo_M1$lower1
seroexpo_M1$upper<-seroexpo_M1$upper1
data_preg<-seroexpo_M1
data_preg$Date<-as.Date(data_preg$t)

data_tot<-read.csv(paste(folder,"/MultiSeroModelResults.csv", sep = ""))
data_tot$Date<-as.Date(as.Date("2020-01-02")+seq(from=0,to=length(data_tot$Date)-1))
data<-merge(data_tot,data_preg,by.x =  "Date",by.y = "Date",all=TRUE)
data$median<-as.numeric(data$median)

x0<-data$Date[1:(which(data$median>0)[1]-1)]
x1<-data$Expo_NYC.median[1:(which(data$median>0)[1]-1)]
x2<-data$Expo_NYC.95.lower[1:(which(data$median>0)[1]-1)]
x3<-data$Expo_NYC.95.upper[1:(which(data$median>0)[1]-1)]

dataNEW<-data[which(data$median>0),]

xx0<-c(x0,dataNEW$Date)
xx1<-c(x1,dataNEW$Expo_NYC.median)
xx2<-c(x2,dataNEW$Expo_NYC.95.lower)
xx3<-c(x3,dataNEW$Expo_NYC.95.upper)

yy1_1<-c(rep(NA,length(x1)),dataNEW$median)
yy2_1<-c(rep(NA,length(x2)),dataNEW$lower)
yy3_1<-c(rep(NA,length(x3)),dataNEW$upper)
yy4_1<-c(rep(NA,length(x2)),dataNEW$lower1)
yy5_1<-c(rep(NA,length(x3)),dataNEW$upper1)
yy6_1<-c(rep(NA,length(x2)),dataNEW$lower2)
yy7_1<-c(rep(NA,length(x3)),dataNEW$upper2)

date<-as.Date("2020-04-20")+c(0,7*seq(35))
data_result = data.frame(output = c(rep("Exposure estimated by general population", length(xx0)), 
                                    rep("Exposure estimated by pregnant women in Model 4", length(xx0)),
                                    rep("Exposure estimated by pregnant women in Model 3", length(xx0)),
                                    rep("Exposure estimated by pregnant women in Model 2", length(xx0)),
                                    rep("Exposure estimated by pregnant women in Model 1", length(xx0))), 
                         t=c(as.Date(xx0),as.Date(xx0),as.Date(xx0),as.Date(xx0),as.Date(xx0)), 
                         median = c(xx1,yy1_4,yy1_3,yy1_2,yy1_1),
                         lower  = c(xx2,yy2_4,yy2_3,yy2_2,yy2_1), 
                         upper  = c(xx3,yy3_4,yy3_3,yy2_3,yy2_1),
                         lower1  = c(xx2,yy4_4,yy4_3,yy4_2,yy4_1), 
                         upper1  = c(xx3,yy5_4,yy5_3,yy5_2,yy5_1),
                         lower2  = c(xx2,yy6_4,yy6_3,yy6_2,yy6_1), 
                         upper2  = c(xx3,yy7_4,yy7_3,yy7_2,yy7_1))


p4<-ggplot(data_result, aes(x=t, y = median, group = output, colour = output)) +
  geom_line(size = 1) +  ggtitle("")+xlab("")+ylab("Proportion of population")+
  # geom_ribbon(aes(ymin=lower, ymax=upper, fill = output), alpha=0.2, colour = NA)+
  #geom_ribbon(aes(ymin=lower1, ymax=upper1, fill = output), alpha=0.1, colour = NA)+
  geom_ribbon(aes(ymin=lower2, ymax=upper2, fill = output), alpha=0.2, colour = NA)+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6), limit = c(0, 0.6))+
  scale_x_date(breaks = c(as.Date(c("2020-01-01","2020-03-01","2020-05-01","2020-07-01","2020-09-01","2020-11-01","2021-01-01"))),
               labels = c("Jan 2020","Mar","May","Jul","Sep","Nov","Jan 2021"),
               limits = as.Date(c("2020-01-01","2021-01-01")))+
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
    # axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black",size = 1),
    plot.margin = margin(t=0, r=right_margin, b=0, l=0, "cm")
  ) 

p4






