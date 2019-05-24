rm(list=ls(all=TRUE))
library("haven")
library(lcrisks)
library(ggplot2)
library(plotrix)
require(gridExtra)
nlst_new <- read_sas("~/Desktop/Lung cancer/lrisk/sasfile/nlst/nlst364_prsn_20171023.sas7bdat")
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/anlst_mod.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/mortality.model.v4.RData")
load("~/Desktop/Lung cancer/lrisk/other/nlst/mortalityrisk_pop_012715.Rdata")


anlst$t0 <- ifelse(is.na(anlst$truefalse_scrnres_ly0)|anlst$truefalse_scrnres_ly0==0,
             0,anlst$truefalse_scrnres_ly0)  # EXCLUDE NO/INADEQUATE
anlst$t0 <- ifelse(anlst$t0==0,"No Screen",
             ifelse(anlst$t0<=2,"True Positive",
                    ifelse(anlst$t0==3,"False Positive","Negative")))

anlst$t1 <- ifelse(is.na(anlst$truefalse_scrnres_ly1)|anlst$truefalse_scrnres_ly1==0,
             0,anlst$truefalse_scrnres_ly1)  # EXCLUDE NO/INADEQUATE
anlst$t1 <- ifelse(anlst$t1==0,"No Screen",
             ifelse(anlst$t1<=2,"True Positive",
                    ifelse(anlst$t1==3,"False Positive","Negative")))

anlst$t2 <- ifelse(is.na(anlst$truefalse_scrnres_ly2)|anlst$truefalse_scrnres_ly2==0,
             0,anlst$truefalse_scrnres_ly2)  # EXCLUDE NO/INADEQUATE
anlst$t2 <- ifelse(anlst$t2==0,"No Screen",
             ifelse(anlst$t2<=2,"True Positive",
                    ifelse(anlst$t2==3,"False Positive","Negative")))

# CUMULATIVE FALSE-POSITIVE
# 0= no screen
# 1= true positive
# 2= false positive
# 3= negative
anlst$false_pos <- ifelse(anlst$t0=="True Positive"|anlst$t1=="True Positive"|anlst$t2=="True Positive",1,
                    ifelse(anlst$t0=="False Positive"|anlst$t1=="False Positive"|anlst$t2=="False Positive",2,
                           ifelse(anlst$t0=="No Screen"&anlst$t1=="No Screen"&anlst$t2=="No Screen",0,3)))

anlst$totalfalsepos <- as.numeric(anlst$t0=="False Positive") + as.numeric(anlst$t1=="False Positive") + as.numeric(anlst$t2=="False Positive")

nlst_new2 <- nlst_new[c("pid","rand_month","rand_year")]
anlst <- merge(anlst,nlst_new2,by="pid")
anlst$days_rand_studyend <- 365*(2008-anlst$rand_year)+365/12*(13-anlst$rand_month)
nlst_old <- mortalityrisk_pop_012715[c("pid","educat")]
anlst <- merge(anlst,nlst_old,by="pid")
anlst$edu6<-ifelse((anlst$educat==1)|(anlst$educat==2),1,
             ifelse((anlst$educat==3)|(anlst$educat==8)|(anlst$educat==98)|(anlst$educat==99),2,
                    ifelse(anlst$educat==4,3,
                           ifelse(anlst$educat==5,4,
                                  ifelse((anlst$educat==6),5,
                                         ifelse(anlst$educat==7,6,0))))))


anlst <- subset(anlst,death.years>=0) #remove those with negative death years
persons <- data.frame(age=anlst$age,
                      female=anlst$female,
                      smkyears=anlst$smkyears,
                      qtyears=anlst$qtyears,
                      cpd=anlst$cpd,
                      race=as.numeric(as.character(anlst$race)),
                      emp=anlst$emp,
                      fam.lung.trend=anlst$fam.lung.trend,
                      bmi=anlst$bmi,
                      edu6=anlst$edu6)

missval <- ifelse(!is.na(persons$age) & !is.na(persons$female) & !is.na(persons$smkyears) & !is.na(persons$qtyears) &
                    !is.na(persons$cpd) & !is.na(persons$race) & !is.na(persons$emp) &
                    !is.na(persons$fam.lung.trend) & !is.na(persons$bmi) & !is.na(persons$edu6),0,1)

anlst$lcrat[missval==0] <- lcrisk(persons[missval==0,],5)[,5]/1000

quant.cp <- quantile(anlst$lcrat,c(.2,.4,.6,.8))
anlst$lcrat.cat <- ifelse(anlst$lcrat<quant.cp[1],1,
                          ifelse(anlst$lcrat<quant.cp[2],2,
                                 ifelse(anlst$lcrat<quant.cp[3],3,
                                        ifelse(anlst$lcrat<quant.cp[4],4,5))))

#Average Life-Years Gained
obs.lyg <- function(ct,xray){
  kmm_ct <- survfit(Surv(ct$death.years_6yr,ct$death_6yr)~1)
  kmm_xray <- survfit(Surv(xray$death.years_6yr,xray$death_6yr)~1)

   overall_ly_ct <- sum(diff(kmm_ct$time)*kmm_ct$surv[2:length(kmm_ct$surv)])
   overall_ly_xray <- sum(diff(kmm_xray$time)*kmm_xray$surv[2:length(kmm_xray$surv)])
   overall_lyg <- overall_ly_ct-overall_ly_xray  #overall life years gained from screening
   return(c(overall_ly_ct,overall_ly_xray,overall_lyg,mean(rbind(ct,xray)$lyg)))
}

#last death is between 6.1-6.2 years for CT/X-ray arms.  Cut-off death at 6 years
anlst$death_6yr <- ifelse(anlst$death==1 & anlst$death.years>=6,0,anlst$death)
anlst$death.years_6yr <- ifelse(anlst$death.years>=6,6,anlst$death.years)
ct <- subset(anlst,screen_group=="CT")
xray <- subset(anlst,screen_group=="X-ray")
kmm_ct <- survfit(Surv(ct$death.years_6yr,ct$death_6yr)~1)
kmm_xray <- survfit(Surv(xray$death.years_6yr,xray$death_6yr)~1)

overall_ly_ct_neg <- sum(diff(kmm_ct$time[kmm_ct$time<=1.711156])*(kmm_ct$surv[2:sum(kmm_ct$time<=1.711156)]))
overall_ly_xray_neg <- sum(diff(kmm_xray$time[kmm_xray$time<=1.711156])*(kmm_xray$surv[2:sum(kmm_xray$time<=1.711156)]))
overall_lyg_neg <- overall_ly_ct_neg-overall_ly_xray_neg  #overall life years gained from screening
365.25*overall_lyg_neg
overall_ly_ct_pos <- sum(diff(kmm_ct$time[sum(kmm_ct$time<=1.711156):length(kmm_ct$time)])*(kmm_ct$surv[(min(which(kmm_ct$time>1.711156))):length(kmm_ct$time)]))
overall_ly_xray_pos <- sum(diff(kmm_xray$time[sum(kmm_xray$time<=1.711156):length(kmm_xray$time)])*(kmm_xray$surv[(min(which(kmm_xray$time>1.711156))):length(kmm_xray$time)]))
overall_lyg_pos <- overall_ly_ct_pos-overall_ly_xray_pos  #overall life years gained from screening
365.25*overall_lyg_pos
365.25*obs.lyg(ct,xray)

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/kmm_diff.jpeg",width=11.5,height=8,units='in',res=900)
plot(kmm_ct$time,kmm_ct$surv,type="l",lwd=3,ylim=c(0.9,1),xlab="Years",ylab="Survival",main="Life gained from CT screening in the NLST")
lines(kmm_xray$time,kmm_xray$surv,lwd=3,col="blue")
axis.break(2,.899,style="slash")
xx <- c(kmm_ct$time[kmm_ct$time>1.711156], rev(kmm_xray$time[kmm_xray$time>1.711156]))
yy <- c(kmm_ct$surv[kmm_ct$time>1.711156],rev(kmm_xray$surv[kmm_xray$time>1.711156]))
polygon(xx,yy,col="green")
xx <- c(kmm_ct$time[kmm_ct$time<=1.711156], rev(kmm_xray$time[kmm_xray$time<=1.711156]))
yy <- c(kmm_ct$surv[kmm_ct$time<=1.711156],rev(kmm_xray$surv[kmm_xray$time<=1.711156]))
polygon(xx,yy,col="red")
arrows(x0=3.88, y0=0.9634095, x1 = 3.38, y1 = 0.951643, length = 0.1, angle = 30,
       code = 1, col = "blue")
text(x=3.38,y=0.95,"X-ray arm", col="blue")
arrows(x0=4, y0=min(kmm_ct$surv[kmm_ct$time<=4])+.001, x1 = 4.5, y1 = 0.978, length = 0.1, angle = 30,
       code = 1, col = "black")
text(x=4.5,y=0.98,"CT arm", col="black")
arrows(x0=0.75, y0=0.9963, x1 = 0.5, y1 = 0.982, length = 0.1, angle = 30,
       code = 1)
text(x=0.5,0.98,"-0.3 days",col="red")
arrows(x0=5, y0=0.943, x1 = 5.25, y1 = 0.958, length = 0.1, angle = 30,
       code = 1)
text(x=5.25,0.96,"+4.7 days",col="green")
text(x=2,y=0.92,"+4.4 days of life gained from CT screening after 6 years",col="green")
text(x=6,y=0.9275,"1,871\ndeaths")
text(x=6,y=0.9075,"1,992\ndeaths")
dev.off()

fit <- coxph(Surv(death.years,death)~as.factor(screen_group), data=anlst)
survCT <- survfit(fit,newdata=data.frame(screen_group=as.factor("CT"))) 
survXray <- survfit(fit,newdata=data.frame(screen_group=as.factor("X-ray"))) 

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/cox_diff.jpeg",width=11.5,height=8,units='in',res=900)
plot(survCT$time,survCT$surv,type="l",lwd=3,ylim=c(0.9,1),xlab="Years",ylab="Survival")
lines(survXray$time,survXray$surv,lwd=3,col="blue")
xx <- c(survCT$time, rev(survXray$time))
yy <- c(survCT$surv,rev(survXray$surv))
polygon(xx,yy,col="green")
dev.off()

overall_lyg <- sum(diff(survCT$time)*((survCT$surv-survXray$surv)[2:length(survCT$surv)]))
365.25*overall_lyg

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/compare_kmm_cox.jpeg",width=11.5,height=8,units='in',res=900)
plot(kmm_ct$time,kmm_ct$surv,type="l",lwd=3,ylim=c(0.9,1),xlab="Years",ylab="Survival",main="Life gained from CT screening in the NLST")
lines(kmm_xray$time,kmm_xray$surv,lwd=3,col="blue")
lines(survCT$time,survCT$surv,lwd=3,col="orange")
lines(survXray$time,survXray$surv,lwd=3,col="green")
dev.off()













ctq1 <- subset(anlst,screen_group=="CT" & lcrat.cat==1)
xrayq1 <- subset(anlst,screen_group=="X-ray" & lcrat.cat==1)
ctq2 <- subset(anlst,screen_group=="CT" & lcrat.cat==2)
xrayq2 <- subset(anlst,screen_group=="X-ray" & lcrat.cat==2)
ctq3 <- subset(anlst,screen_group=="CT" & lcrat.cat==3)
xrayq3 <- subset(anlst,screen_group=="X-ray" & lcrat.cat==3)
ctq4 <- subset(anlst,screen_group=="CT" & lcrat.cat==4)
xrayq4 <- subset(anlst,screen_group=="X-ray" & lcrat.cat==4)
ctq5 <- subset(anlst,screen_group=="CT" & lcrat.cat==5)
xrayq5 <- subset(anlst,screen_group=="X-ray" & lcrat.cat==5)
obslyg <- rbind(obs.lyg(ct,xray),obs.lyg(ctq1,xrayq1),obs.lyg(ctq2,xrayq2),obs.lyg(ctq3,xrayq3),obs.lyg(ctq4,xrayq4),obs.lyg(ctq5,xrayq5))
rownames(obslyg) <- c("Overall","Q1","Q2","Q3","Q4","Q5")
colnames(obslyg) <- c("LY CT","LY X-ray","LYG","Est. Mean LYG")

p1 <- ggplot(ct, aes(x=lcrat.cat,y=365.25*lyg)) +
  geom_boxplot(outlier.size=.5,aes(group=cut_width(lcrat.cat,1))) +
  geom_hline(yintercept=19.8,col="green") + 
  ggtitle("2b. Estimated life gained from screening by risk quintiles in the CT arm") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_text(x=1.25,y=23,label="Average expected life gained=19.8 days") + 
  labs(x = "Risk Quintiles (lowest to highest)", y = "Days of life gained")

trans_new(risk_tran, transform, inverse, breaks = extended_breaks(),
          format = format_format(), domain = c(-Inf, Inf))

p2 <- ggplot(anlst,aes(x=100*lcrat,y=365*lyg), log10="x") + geom_point(alpha = 0.05) + 
      scale_x_log10(breaks=c(0.5,1,2,3,10),limits=c(0.25,10)) +
      scale_y_continuous(limits=c(0,75)) +
      labs(x="% lung cancer risk") + 
      labs(y="Days of life gained") + 
      geom_hline(yintercept=19.8,col="green") +
      geom_vline(xintercept=1.16,col="blue") +
      geom_vline(xintercept=1.77,col="blue") +
      geom_vline(xintercept=2.56,col="blue") +  
      geom_vline(xintercept=4.02,col="blue") +  
      annotate("text",x=0.5,y=21,label="Average expected life gained=19.8 days") +
      annotate("text",x=0.5,y=70,label="Q1") +
      annotate("text",x=1.5,y=70,label="Q2") +
      annotate("text",x=2.25,y=70,label="Q3") +
      annotate("text",x=3.25,y=70,label="Q4") +
      annotate("text",x=7,y=70,label="Q5") +
      ggtitle("2a. Estimated life gained from screening and 5-year lung cancer risks") +
      theme(plot.title = element_text(hjust=0.5))

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_lyg_and_risks.jpeg",width=11.5,height=8,units='in',res=900)
plot(p2)
dev.off()

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_lyg_and_risks2.jpeg",width=11.5,height=8,units='in',res=900)
plot(p1)
dev.off()




quant.cp <- quantile(anlst$lyg,c(.2,.4,.6,.8))
anlst$lyg.cat <- ifelse(anlst$lyg<quant.cp[1],1,
                          ifelse(anlst$lyg<quant.cp[2],2,
                                 ifelse(anlst$lyg<quant.cp[3],3,
                                        ifelse(anlst$lyg<quant.cp[4],4,5))))
fp_risk <- table(anlst$false_pos,anlst$lcrat.cat)
numfp_risk <- table(anlst$totalfalsepos,anlst$lcrat.cat)
fp_lg <- table(anlst$false_pos,anlst$lyg.cat)
numfp_lg <- table(anlst$totalfalsepos,anlst$lyg.cat)

fp <- rbind(fp_risk[3,],fp_lg[3,])
colnames(fp) <- c("Q1","Q2","Q3","Q4","Q5")
rownames(fp) <- c("Risk","Life Gained")
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/persons_falsepos.jpeg",width=11.5,height=8,units='in',res=900)
barCenters <- barplot(height=fp,
                      beside = TRUE, las = 1,
                      ylim=c(0,4000), 
                      ylab = "Persons with false-positive results",
                      border = "black", axes = TRUE,
                      legend.text = TRUE, main="Number of people with false-positive results",
                      args.legend = list(title = "Quintiles", x = "top",cex = 1))
dev.off()

numfp <- rbind(numfp_risk[2,]+2*numfp_risk[3,]+3*numfp_risk[4,],
               numfp_lg[2,]+2*numfp_lg[3,]+3*numfp_lg[4,])
colnames(numfp) <- c("Q1","Q2","Q3","Q4","Q5")
rownames(numfp) <- c("Risk","Life Gained")
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/number_falsepos.jpeg",width=11.5,height=8,units='in',res=900)
barCenters <- barplot(height=numfp,
                      beside = TRUE, las = 1,
                      ylim=c(0,6000), 
                      ylab = "Number of false-positive results",
                      border = "black", axes = TRUE,
                      legend.text = TRUE, main="Number of false-positive results", 
                      args.legend = list(title = "Quintiles", x = "top",cex = 1))
dev.off()

false_true_ratio <- rbind(fp_risk[3,]/fp_risk[2,],fp_lg[3,]/fp_lg[2,])