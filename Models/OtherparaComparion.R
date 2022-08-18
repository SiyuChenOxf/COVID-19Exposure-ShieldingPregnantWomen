#This file is to plot comparisons of estimations of parameters among four models 
#For Fig S4

library(ggpubr)
library(ggplot2)

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

fit_1para <- readRDS(paste(folder,"/fit_1para.rds", sep = ""))
k1_M1<-rstan::extract(fit_1para)$k1
k11_M1<-rstan::extract(fit_1para)$k11
k10_M1<-rstan::extract(fit_1para)$k10
k01_M1<-rstan::extract(fit_1para)$k01

fit_2para <- readRDS(paste(folder,"/fit_2para.rds", sep = ""))
k1_M2<-rstan::extract(fit_2para)$k1
k2_M2<-rstan::extract(fit_2para)$k2
k11_M2<-rstan::extract(fit_2para)$k11
k10_M2<-rstan::extract(fit_2para)$k10
k01_M2<-rstan::extract(fit_2para)$k01

fit_3para <- readRDS(paste(folder,"/fit_3para.rds", sep = ""))
k1_M3<-rstan::extract(fit_3para)$k1
k2_M3<-rstan::extract(fit_3para)$k2
k3_M3<-rstan::extract(fit_3para)$k3
k11_M3<-rstan::extract(fit_1para)$k11
k10_M3<-rstan::extract(fit_1para)$k10
k01_M3<-rstan::extract(fit_1para)$k01

fit_4para <- readRDS(paste(folder,"/fit_4para.rds", sep = ""))
k1_M4<-rstan::extract(fit_4para)$k1
k2_M4<-rstan::extract(fit_4para)$k2
k3_M4<-rstan::extract(fit_4para)$k3
k4_M4<-rstan::extract(fit_4para)$k4
k11_M4<-rstan::extract(fit_4para)$k11
k10_M4<-rstan::extract(fit_4para)$k10
k01_M4<-rstan::extract(fit_4para)$k01

a<-data.frame(k11_M1=k11_M1,
              k11_M2=k11_M2,
              k11_M3=k11_M3,
              k11_M4=k11_M4,
              k10_M1=k10_M1,
              k10_M2=k10_M2,
              k10_M3=k10_M3,
              k10_M4=k10_M4,
              k01_M1=k01_M1,
              k01_M2=k01_M2,
              k01_M3=k01_M3,
              k01_M4=k01_M4 )

data<-data.frame(group=rep(c("Model 1","Model 2","Model 3","Model 4"),3),
                 labels = rep(c("k11", "k10","k01"),each=4),
                 med=c(apply(a, 2, function(x) quantile(x, probs = 0.5))),
                 upp=c(apply(a, 2, function(x) quantile(x, probs = 0.75))),
                 low=c(apply(a, 2, function(x) quantile(x, probs = 0.25))))

data$group <- as.factor(data$group)
data$labels <- as.factor(data$labels)

data <- data %>%
  group_by(labels, group) %>%
  summarise(
    med = med,
    upp = upp,
    low = low
  )
ggplot(data , aes(labels, med)) +ylim(0,1)+xlab("")+ylab("")+
  geom_pointrange(
    aes(ymin = low, ymax = upp, color = group),
    position = position_dodge(0.5)
  )+
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


####
d=data.frame(Model2=c("k1"), mean=c(median(k1_M1)), lower=c(quantile(k1_M1,0.05)), upper=c(quantile(k1_M1,0.95)))
p1<-ggplot() + ylim(c(0,0.025))+xlab("")+ylab("")+
  geom_pointrange(data=d, mapping=aes(x=Model2, y=mean, ymin=upper, ymax=lower))  +
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


d=data.frame(Model2=c("k1","k2"), mean=c(median(k1_M2),median(k2_M2)), lower=c(quantile(k1_M2,0.05),quantile(k2_M2,0.05)), upper=c(quantile(k1_M2,0.95),quantile(k2_M2,0.95)))
p2<-ggplot() + ylim(c(0,0.025))+xlab("")+ylab("")+
  geom_pointrange(data=d, mapping=aes(x=Model2, y=mean, ymin=upper, ymax=lower))  +
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

d=data.frame(Model2=c("k1","k2","k3"), mean=c(median(k1_M3),median(k2_M3),median(k3_M3)), lower=c(quantile(k1_M3,0.05),quantile(k2_M3,0.05),quantile(k3_M3,0.05)), upper=c(quantile(k1_M3,0.95),quantile(k2_M3,0.95),quantile(k3_M3,0.95)))
p3<-ggplot() + ylim(c(0,0.025))+xlab("")+ylab("")+
  geom_pointrange(data=d, mapping=aes(x=Model2, y=mean, ymin=upper, ymax=lower))  +
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

d=data.frame(Model2=c("k1","k2","k3","k4"), mean=c(median(k1_M4),median(k2_M4),median(k3_M4),median(k4_M4)), lower=c(quantile(k1_M4,0.05),quantile(k2_M4,0.05),quantile(k3_M4,0.05),quantile(k4_M4,0.05)), upper=c(quantile(k1_M4,0.95),quantile(k2_M4,0.95),quantile(k3_M4,0.95),quantile(k4_M4,0.95)))
p4<-ggplot() + ylim(c(0,0.025))+xlab("")+ylab("")+
  geom_pointrange(data=d, mapping=aes(x=Model2, y=mean, ymin=upper, ymax=lower))  +
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
ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2,labels = c("A", "B", "B","D"))
