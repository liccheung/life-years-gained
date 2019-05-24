rm(list=ls(all=TRUE))
library("haven")
library(lcrisks)
library(ggplot2)
library(plotrix)
require(gridExtra)
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/anlst_mod.RData")

#last death is 7.16/7.17 years for CT/X-ray arms.  Cut-off death at 7 years
anlst$death_7yr <- ifelse(anlst$death==1 & anlst$death.years>=7,0,anlst$death)
anlst$death.years_7yr <- ifelse(anlst$death.years>=7,7,anlst$death.years)
ct <- subset(anlst,screen_group=="CT")
xray <- subset(anlst,screen_group=="X-ray")
kmm_ct <- survfit(Surv(ct$death.years_7yr,ct$death_7yr)~1)
kmm_xray <- survfit(Surv(xray$death.years_7yr,xray$death_7yr)~1)

compare_ctsurv <- rep(NA,length(kmm_xray$time))
for(i in 1:length(kmm_xray$time)){
  compare_ctsurv[i] <- min(kmm_ct$surv[kmm_ct$time<=kmm_xray$time[i]]) 
}
which(compare_ctsurv<kmm_xray$surv)
switch_time <- kmm_xray$time[365]

#Average Life-Years Gained
obs.lyg <- function(ct,xray){
  kmm_ct <- survfit(Surv(ct$death.years_7yr,ct$death_7yr)~1)
  kmm_xray <- survfit(Surv(xray$death.years_7yr,xray$death_7yr)~1)
  
  overall_ly_ct <- sum(diff(c(0,kmm_ct$time))*(kmm_ct$surv))
  overall_ly_xray <- sum(diff(c(0,kmm_xray$time))*(kmm_xray$surv))
  overall_lyg <- overall_ly_ct-overall_ly_xray  #overall life years gained from screening
  return(c(overall_ly_ct,overall_ly_xray,overall_lyg,mean(rbind(ct,xray)$lyg)))
}

overall_ly_ct_neg <- sum(diff(c(0,kmm_ct$time[kmm_ct$time<=switch_time]))*(kmm_ct$surv[kmm_ct$time<=switch_time]))
overall_ly_xray_neg <- sum(diff(c(0,kmm_xray$time[kmm_xray$time<=switch_time]))*(kmm_xray$surv[kmm_xray$time<=switch_time]))
overall_lyg_neg <- overall_ly_ct_neg-overall_ly_xray_neg  #overall life years gained from screening
365.25*overall_lyg_neg
overall_ly_ct_pos <- sum(diff(kmm_ct$time[sum(kmm_ct$time<=switch_time):length(kmm_ct$time)])*(kmm_ct$surv[(min(which(kmm_ct$time>switch_time))):length(kmm_ct$time)]))
overall_ly_xray_pos <- sum(diff(kmm_xray$time[sum(kmm_xray$time<=switch_time):length(kmm_xray$time)])*(kmm_xray$surv[(min(which(kmm_xray$time>switch_time))):length(kmm_xray$time)]))
overall_lyg_pos <- overall_ly_ct_pos-overall_ly_xray_pos  #overall life years gained from screening
365.25*overall_lyg_pos
365.25*obs.lyg(ct,xray)

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/kmm_diff_v2.jpeg",width=11.5,height=8,units='in',res=900)
plot(kmm_ct$time,kmm_ct$surv,type="l",lwd=3,ylim=c(0.9,1),xlab="Years",ylab="Survival",main="Life gained from CT screening in the NLST")
lines(kmm_xray$time,kmm_xray$surv,lwd=3,col="blue")
axis.break(2,.899,style="slash")
xx <- c(kmm_ct$time[kmm_ct$time>switch_time], rev(kmm_xray$time[kmm_xray$time>switch_time]))
yy <- c(kmm_ct$surv[kmm_ct$time>switch_time],rev(kmm_xray$surv[kmm_xray$time>switch_time]))
polygon(xx,yy,col="green")
xx <- c(kmm_ct$time[kmm_ct$time<=switch_time], rev(kmm_xray$time[kmm_xray$time<=switch_time]))
yy <- c(kmm_ct$surv[kmm_ct$time<=switch_time],rev(kmm_xray$surv[kmm_xray$time<=switch_time]))
polygon(xx,yy,col="red")
arrows(x0=3.88, y0=0.9634095, x1 = 3.38, y1 = 0.951643, length = 0.1, angle = 30,
       code = 1, col = "blue")
text(x=3.38,y=0.95,"X-ray arm", col="blue")
arrows(x0=4, y0=min(kmm_ct$surv[kmm_ct$time<=4])+.001, x1 = 4.5, y1 = 0.978, length = 0.1, angle = 30,
       code = 1, col = "black")
text(x=4.5,y=0.98,"CT arm", col="black")
arrows(x0=0.75, y0=0.9963, x1 = 0.5, y1 = 0.982, length = 0.1, angle = 30,
       code = 1)
text(x=0.5,0.98,"-1.5 days",col="red")
arrows(x0=5, y0=0.951, x1 = 5.25, y1 = 0.965, length = 0.1, angle = 30,
       code = 1)
text(x=5.25,0.967,"+6.8 days",col="green")
text(x=2,y=0.92,"+5.3 days of life gained from CT screening after 7 years",col="green")
text(x=7,y=0.9275,"1,892\ndeaths")
text(x=7,y=0.9075,"2,005\ndeaths", col="blue")
dev.off()

load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/anlst_mod2.RData")
ct <- subset(anlst,screen_group=="CT")
xray <- subset(anlst,screen_group=="X-ray")
times <- seq(0.125,25.625,by=.125)
ctsurv <- colMeans(surv[anlst$screen_group=="CT",]-
  0.204*cbind(anlst[anlst$screen_group=="CT",99:138],matrix(rep(anlst[anlst$screen_group=="CT",138],165),nrow=sum(anlst$screen_group=="CT"))))
xraysurv <- colMeans(surv[anlst$screen_group=="X-ray",])
plot(times,ctsurv,type="l")
lines(times,xraysurv,col="blue")
lines(kmm_ct$time,kmm_ct$surv,col="red")
lines(kmm_xray$time,kmm_xray$surv,col="green")
