#This code is to clean and describe raw data before modelling

rm(list = ls())

folder_strings = unlist(strsplit(getwd(), '/'))
folder_strings[length(folder_strings)] = "Data"
folder = paste(folder_strings, sep = "", collapse = "/")
# setwd("C:/Users/schen/Desktop/MfPW/Rcodes/NEWWWW")

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

#Read raw data
data<-read.csv(paste(folder,"/Prevalence_PCR_SER_encrypted.csv", sep = ""),header = TRUE) 
data$Date_Birth<-as.Date(data$Date_Birth)
dataNEW<-data[order(data$Date_Birth),]              #sort the data by birth date
dataNEW<-dataNEW[dataNEW$Date_Birth<="2021-02-15",] #truncated before high vaccination rate

#Categorizing pcr & sero results for every pregnant women
data6<-data.frame(Date_Birth=dataNEW$Date_Birth, sero=dataNEW$BirthSER_M,pcr=dataNEW$BirthPCR_M,res=rep(0,length(dataNEW$Date_Birth)))
for (i in 1:length(data6$Date_Birth)) {
  if(data6$pcr[i]=="Negative" & data6$sero[i]=="Negative"){
    data6$res[i]="pcr Negative & sero Negative"
  }
  else if(data6$pcr[i]=="Negative" & data6$sero[i]=="Never_tested"){
    data6$res[i]="pcr Negative & sero Missing"
  }
  else if(data6$pcr[i]=="Negative" & data6$sero[i]=="No_Birth_test"){
    data6$res[i]="pcr Negative & sero Missing"
  }
  else if(data6$pcr[i]=="Negative" & data6$sero[i]=="Positive"){
    data6$res[i]="pcr Negative & sero Positive"
  }
  else if(data6$pcr[i]=="Not_tested" & data6$sero[i]=="Negative"){
    data6$res[i]="pcr Missing & sero Negative"
  }
  else if(data6$pcr[i]=="Not_tested" & data6$sero[i]=="Never_tested"){
    data6$res[i]="pcr Missing & sero Missing"
  }
  else if(data6$pcr[i]=="Not_tested" & data6$sero[i]=="No_Birth_test"){
    data6$res[i]="pcr Missing & sero Missing"
  }
  else if(data6$pcr[i]=="Not_tested" & data6$sero[i]=="Positive"){
    data6$res[i]="pcr Missing & sero Positive"
  }
  else if(data6$pcr[i]=="Positive" & data6$sero[i]=="Negative"){
    data6$res[i]="pcr Positive & sero Negative"
  }
  else if(data6$pcr[i]=="Positive" & data6$sero[i]=="Never_tested"){
    data6$res[i]="pcr Positive & sero Missing"
  }
  else if(data6$pcr[i]=="Positive" & data6$sero[i]=="No_Birth_test"){
    data6$res[i]="pcr Positive & sero Missing"
  }
  else if(data6$pcr[i]=="Positive" & data6$sero[i]=="Positive"){
    data6$res[i]="pcr Positive & sero Positive"
  }
}
#By date pcr & sero results
data7<-data.frame(Date_Birth=unique(dataNEW$Date_Birth),A=rep(0,length(unique(dataNEW$Date_Birth))),B=rep(0,length(unique(dataNEW$Date_Birth))),C=rep(0,length(unique(dataNEW$Date_Birth))),D=rep(0,length(unique(dataNEW$Date_Birth))),E=rep(0,length(unique(dataNEW$Date_Birth))),F=rep(0,length(unique(dataNEW$Date_Birth))),G=rep(0,length(unique(dataNEW$Date_Birth))),H=rep(0,length(unique(dataNEW$Date_Birth))),I=rep(0,length(unique(dataNEW$Date_Birth))))
for (i in 1:length(unique(dataNEW$Date_Birth))) {
  c<-which(data6$Date_Birth==unique(dataNEW$Date_Birth)[i])
  for (j in 1:length(c)) {
    if(data6$res[c][j]=="pcr Negative & sero Missing"){
      data7[i,2]=data7[i,2]+1
    }
    else if(data6$res[c][j]=="pcr Missing & sero Missing"){
      data7[i,3]=data7[i,3]+1
    }
    else if(data6$res[c][j]=="pcr Positive & sero Missing"){
      data7[i,4]=data7[i,4]+1
    }
    else if(data6$res[c][j]=="pcr Negative & sero Positive"){
      data7[i,5]=data7[i,5]+1
    }
    else if(data6$res[c][j]=="pcr Missing & sero Positive"){
      data7[i,6]=data7[i,6]+1
    }
    else if(data6$res[c][j]=="pcr Positive & sero Positive"){
      data7[i,7]=data7[i,7]+1
    }
    else if(data6$res[c][j]=="pcr Negative & sero Negative"){
      data7[i,8]=data7[i,8]+1
    }
    else if(data6$res[c][j]=="pcr Missing & sero Negative"){
      data7[i,9]=data7[i,9]+1
    }
    else if(data6$res[c][j]=="pcr Positive & sero Negative"){
      data7[i,10]=data7[i,10]+1
    }
  }
}
data_trun<-data7[4:188,]
sum(data_trun[,-1],na.rm = TRUE) #2682

x<-data.frame(Date_Birth=as.Date(as.Date("2020-4-16")+seq(1,365)),
              x1=rep(0,365))

b<-merge(x,data7,by="Date_Birth",all=T)
sum(b[4:255,-c(1:2)],na.rm = TRUE) #2682

#Missing any record in that date
b[which(b$Date_Birth==as.Date("2020-07-27")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-08-08")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-08-09")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-08-12")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-08-22")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-08-23")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-08-28")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-09-10")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-09-16")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-09-26")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-09-27")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-10-06")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-10-15")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-10-16")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-10-26")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-10-30")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-11-01")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-11-04")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-11-10")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-11-11")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-11-21")),3:11]<-0
b[which(b$Date_Birth==as.Date("2020-11-22")),3:11]<-0
sum(b[4:255,-c(1:2)],na.rm = TRUE) #2682
# write.csv(b,'wer.csv')

#Transfer daily count data into calendar week (2020) proportion 
index<-c(4,4+seq(35)*7)
Sumweekly<-matrix(0,nrow = length(index),ncol = 9)
for (i in 1:length(index)) {
  Sumweekly[i,]<-apply(b[index[i]:(index[i]+6),3:11],2,sum)
}
sum(Sumweekly,na.rm = TRUE)

date<-as.Date(b[index,1])
xx<-data.frame(date=date,A=Sumweekly[,1],B=Sumweekly[,2],C=Sumweekly[,3],D=Sumweekly[,4],E=Sumweekly[,5],F=Sumweekly[,6],G=Sumweekly[,7],H=Sumweekly[,8],I=Sumweekly[,9])
xx[c(12,13),2:10]<-NA                                  #row 12 and 13 A column is too much which means serology data has too many missing points at that date!
xx<-xx[,c("date","D","F","G","I")]
colnames(xx)<-c("Date","PnSp","PpSp","PnSn","PpSn")
prop_xx<-data.frame(Date=as.Date(xx[,"Date"]), prop_PnSp=xx[,2]/apply(xx[2:5],1,sum),prop_PpSp=xx[,3]/apply(xx[2:5],1,sum),prop_PnSn=xx[,4]/apply(xx[2:5],1,sum), prop_PpSn=xx[,5]/apply(xx[2:5],1,sum))
y<-xx[c(1:5,14:21,24:36),]
y<-cbind(y[,"PnSn"],y[,"PpSn"],y[,"PpSp"],y[,"PnSp"])   
saveRDS(y,"data.rds")

n_week2<-length(y[,1]) #number of calender week of data

#To save and plot 4-categorized data and its associated credible interval
median<-upp<-low<-matrix(0,nrow = n_week2,4)
for (i in 1:n_week2) {
 
  data_mn <-  list(N_trials = sum(y[i,]),ans = y[i,])
  fit_mn<-stan(file=paste(folder,"/DriMultinorm.stan", sep = ""),
                     data = data_mn,
                     iter = 10000,
                     chains = 4,
                     seed = 4,
                     control = list(adapt_delta = 0.99,max_treedepth = 15))
  pos_theta <-  rstan::extract(fit_mn)$theta

  median[i,]<-apply(pos_theta, 2, median) 
  upp[i,]<-apply(pos_theta, 2, function(x) quantile(x, probs = 0.95)) 
  low[i,]<-apply(pos_theta, 2, function(x) quantile(x, probs = 0.05)) 

  print(i)
}

data00<-data.frame(t=as.Date(xx$Date),
                         med=c(median[1:5,1],rep(NA,8),median[6:13,1],rep(NA,2),median[14:26,1]),
                         lower=c(low[1:5,1],rep(NA,8),low[6:13,1],rep(NA,2),low[14:26,1]),
                         upper=c(upp[1:5,1],rep(NA,8),upp[6:13,1],rep(NA,2),upp[14:26,1]))
p00<-ggplot(data00,aes(x=t,y=med))+
  geom_point(size = 0.1) +  ggtitle("PCR-Sero-")+ ylab(" Percentage (%)  ") +  xlab("  ")+
  geom_pointrange(aes(ymin=lower, ymax=upper),colour = "black", fill = colors_Dark[2])+
  scale_x_date(breaks = as.Date(data00$t), labels=as.Date(data00$t), limit = as.Date(c("2020-04-20","2020-12-21")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))


data10<-data.frame(t=as.Date(xx$Date),
                   med=c(median[1:5,2],rep(NA,8),median[6:13,2],rep(NA,2),median[14:26,2]),
                   lower=c(low[1:5,2],rep(NA,8),low[6:13,2],rep(NA,2),low[14:26,2]),
                   upper=c(upp[1:5,2],rep(NA,8),upp[6:13,2],rep(NA,2),upp[14:26,2]))
p10<-ggplot(data10,aes(x=t,y=med))+
  geom_point(size = 0.1) +  ggtitle("PCR+Sero-")+ ylab(" Percentage (%)  ") +  xlab("  ")+
  geom_pointrange(aes(ymin=lower, ymax=upper),colour = "black", fill = colors_Dark[2])+
  scale_x_date(breaks = as.Date(data00$t), labels=as.Date(data00$t), limit = as.Date(c("2020-04-20","2020-12-21")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

data11<-data.frame(t=as.Date(xx$Date),
                   med=c(median[1:5,3],rep(NA,8),median[6:13,3],rep(NA,2),median[14:26,3]),
                   lower=c(low[1:5,3],rep(NA,8),low[6:13,3],rep(NA,2),low[14:26,3]),
                   upper=c(upp[1:5,3],rep(NA,8),upp[6:13,3],rep(NA,2),upp[14:26,3]))
p11<-ggplot(data11,aes(x=t,y=med))+
  geom_point(size = 0.1) +  ggtitle("PCR+Sero+")+ ylab(" Percentage (%)  ") +  xlab("  ")+
  geom_pointrange(aes(ymin=lower, ymax=upper),colour = "black", fill = colors_Dark[2])+
  scale_x_date(breaks = as.Date(data00$t), labels=as.Date(data00$t), limit = as.Date(c("2020-04-20","2020-12-21")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

data01<-data.frame(t=as.Date(xx$Date),
                   med=c(median[1:5,4],rep(NA,8),median[6:13,4],rep(NA,2),median[14:26,4]),
                   lower=c(low[1:5,4],rep(NA,8),low[6:13,4],rep(NA,2),low[14:26,4]),
                   upper=c(upp[1:5,4],rep(NA,8),upp[6:13,4],rep(NA,2),upp[14:26,4]))
p01<-ggplot(data01,aes(x=t,y=med))+
  geom_point(size = 0.1) +  ggtitle("PCR-Sero+")+ ylab(" Percentage (%)  ") +  xlab("  ")+
  geom_pointrange(aes(ymin=lower, ymax=upper),colour = "black", fill = colors_Dark[2])+
  scale_x_date(breaks = as.Date(data00$t), labels=as.Date(data00$t), limit = as.Date(c("2020-04-20","2020-12-21")))+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))

ggarrange(p00, p10, p11,p01,
          ncol = 2, nrow = 2)

saveRDS(median,"Result_median_4.rds")
saveRDS(upp,"Result_upp_4.rds")
saveRDS(low,"Result_low_4.rds")

####save and plot seroprevalence data and its associated credible interval
median<-upp<-low<-matrix(0,nrow = n_week2,1)
for (i in 1:n_week2) {
  
  data_list <- list(N = apply(y, 1,sum)[i],  n = apply(y[,c(3,4)],1,sum)[i])
  
  fit_mn<-stan(file=paste(folder,"/BinomialUniform.stan", sep = ""),
                   data = data_list,
                   iter = 10000,
                   chains = 4,
                   seed = 4,
                   control = list(adapt_delta = 0.99,max_treedepth = 15))
  pos_theta <-  rstan::extract(fit_mn)$theta
  
  median[i,]<-quantile(pos_theta, 0.5)
  upp[i,]<-quantile(pos_theta, 0.95)
  low[i,]<-quantile(pos_theta, 0.05)
  
  print(i)
}

A<-c(median[1:5,1],rep(NA,8),median[6:13,1],rep(NA,2),median[14:26,1])
B<-c(low[1:5,1],rep(NA,8),low[6:13,1],rep(NA,2),low[14:26,1])
C<-c(upp[1:5,1],rep(NA,8),upp[6:13,1],rep(NA,2),upp[14:26,1])

data_seropre_another<- data.frame(med=A,low=B,upp=C)
data_seropre_another$time<-as.Date(data11$t)
ggplot(data_seropre_another,aes(x=time,y=med))+
  geom_point(size = 0.1) +  ggtitle("seroprevalence")+ ylab(" Percentage (%)  ") +  xlab("  ")+
  geom_pointrange(aes(ymin=low, ymax=upp),colour = "black", fill = colors_Dark[2])

saveRDS(median,"Result_median_2.rds")
saveRDS(upp,"Result_upp_2.rds")
saveRDS(low,"Result_low_2.rds")
