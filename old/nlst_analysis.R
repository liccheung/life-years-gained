rm(list=ls(all=TRUE))
library("haven")
library(lcrisks)
library(ggplot2)
require(gridExtra)
nlst_new <- read_sas("~/Desktop/Lung cancer/lrisk/sasfile/nlst/nlst364_prsn_20171023.sas7bdat")
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/anlst_mod.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/mortality.model.v4.RData")
load("~/Desktop/Lung cancer/lrisk/other/nlst/mortalityrisk_pop_012715.Rdata")

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
  kmm_ct <- survfit(Surv(ct$death.years,ct$death)~1)
  kmm_xray <- survfit(Surv(xray$death.years,xray$death)~1)

   overall_ly_ct <- sum(diff(kmm_ct$time)*kmm_ct$surv[2:length(kmm_ct$surv)])
   overall_ly_xray <- sum(diff(kmm_xray$time)*kmm_xray$surv[2:length(kmm_xray$surv)])
   overall_lyg <- overall_ly_ct-overall_ly_xray  #overall life years gained from screening
   return(c(overall_ly_ct,overall_ly_xray,overall_lyg,mean(rbind(ct,xray)$lyg)))
}

ct <- subset(anlst,screen_group=="CT")
xray <- subset(anlst,screen_group=="X-ray")
kmm_ct <- survfit(Surv(ct$death.years,ct$death)~1)
kmm_xray <- survfit(Surv(xray$death.years,xray$death)~1)
plot(kmm_ct$time,kmm_ct$surv,type="l",ylim=c(0.9,1))
lines(kmm_xray$time,kmm_xray$surv,col="red")

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
  geom_hline(yintercept=16.8) + 
  ggtitle("Estimated Life Gained from Screening in CT arm") +
  theme(plot.title = element_text(hjust=0.5)) +
  geom_text(x=1.25,y=20,label="Average gain in NLST=16.8 days") + 
  labs(x = "Risk Quintiles (lowest to highest)", y = "Days of life gained")

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_estlyg_dist.jpeg",width=11.5,height=8,units='in',res=900)
plot(p1)
dev.off()

ct55 <- subset(anlst,screen_group=="CT" & age<=55)
xray55 <- subset(anlst,screen_group=="X-ray" & age<=55)
ct56 <- subset(anlst,screen_group=="CT" & age==56)
xray56 <- subset(anlst,screen_group=="X-ray" & age==56)
ct57 <- subset(anlst,screen_group=="CT" & age==57)
xray57 <- subset(anlst,screen_group=="X-ray" & age==57)
ct58 <- subset(anlst,screen_group=="CT" & age==58)
xray58 <- subset(anlst,screen_group=="X-ray" & age==58)
ct59 <- subset(anlst,screen_group=="CT" & age==59)
xray59 <- subset(anlst,screen_group=="X-ray" & age==59)
ct60 <- subset(anlst,screen_group=="CT" & age==60)
xray60 <- subset(anlst,screen_group=="X-ray" & age==60)
ct61 <- subset(anlst,screen_group=="CT" & age==61)
xray61 <- subset(anlst,screen_group=="X-ray" & age==61)
ct62 <- subset(anlst,screen_group=="CT" & age==62)
xray62 <- subset(anlst,screen_group=="X-ray" & age==62)
ct63 <- subset(anlst,screen_group=="CT" & age==63)
xray63 <- subset(anlst,screen_group=="X-ray" & age==63)
ct64 <- subset(anlst,screen_group=="CT" & age==64)
xray64 <- subset(anlst,screen_group=="X-ray" & age==64)
ct65 <- subset(anlst,screen_group=="CT" & age==65)
xray65 <- subset(anlst,screen_group=="X-ray" & age==65)
ct66 <- subset(anlst,screen_group=="CT" & age==66)
xray66 <- subset(anlst,screen_group=="X-ray" & age==66)
ct67 <- subset(anlst,screen_group=="CT" & age==67)
xray67 <- subset(anlst,screen_group=="X-ray" & age==67)
ct68 <- subset(anlst,screen_group=="CT" & age==68)
xray68 <- subset(anlst,screen_group=="X-ray" & age==68)
ct69 <- subset(anlst,screen_group=="CT" & age==69)
xray69 <- subset(anlst,screen_group=="X-ray" & age==69)
ct70 <- subset(anlst,screen_group=="CT" & age==70)
xray70 <- subset(anlst,screen_group=="X-ray" & age==70)
ct71 <- subset(anlst,screen_group=="CT" & age==71)
xray71 <- subset(anlst,screen_group=="X-ray" & age==71)
ct72 <- subset(anlst,screen_group=="CT" & age==72)
xray72 <- subset(anlst,screen_group=="X-ray" & age==72)
ct73 <- subset(anlst,screen_group=="CT" & age==73)
xray73 <- subset(anlst,screen_group=="X-ray" & age==73)
ct74 <- subset(anlst,screen_group=="CT" & age>=74)
xray74 <- subset(anlst,screen_group=="X-ray" & age>=74)
obslyg <- rbind(obs.lyg(ct55,xray55),obs.lyg(ct56,xray56),obs.lyg(ct57,xray57),obs.lyg(ct58,xray58),obs.lyg(ct59,xray59),
                obs.lyg(ct60,xray60),obs.lyg(ct61,xray61),obs.lyg(ct62,xray62),obs.lyg(ct63,xray63),obs.lyg(ct64,xray64),
                obs.lyg(ct65,xray65),obs.lyg(ct66,xray66),obs.lyg(ct67,xray67),obs.lyg(ct68,xray68),obs.lyg(ct69,xray69),
                obs.lyg(ct70,xray70),obs.lyg(ct71,xray71),obs.lyg(ct72,xray72),obs.lyg(ct73,xray73),obs.lyg(ct74,xray74))

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_lyg_age.jpeg",width=11.5,height=8,units='in',res=900)
plot(seq(55,74,1),365*obslyg[,4],ylim=c(-150,150),type="l",xlab="Age at first CT screen",ylab="Days of life gained",main="Life Gained from Screening")
lines(seq(55,74,1),365*obslyg[,3],col="red")
dev.off()

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_estlyg_age.jpeg",width=11.5,height=8,units='in',res=900)
plot(seq(55,74,1),365*obslyg[,4],ylim=c(0,35),type="l",xlab="Age at first CT screen",ylab="Days of life gained",main="Estimated Life Gained from Screening")
dev.off()

#ct55_59 <- subset(anlst,screen_group=="CT" & age<=59)
#xray55_59 <- subset(anlst,screen_group=="X-ray" & age<=59)
#ct60_64 <- subset(anlst,screen_group=="CT" & age>=60 & age<=64)
#xray60_64 <- subset(anlst,screen_group=="X-ray" & age>=60 & age<=64)
#ct65_69 <- subset(anlst,screen_group=="CT" & age>=65 & age<=69)
#xray65_69 <- subset(anlst,screen_group=="X-ray" & age>=65 & age<=69)
#ct70_74 <- subset(anlst,screen_group=="CT" & age>=70 & age<=74)
#xray70_74 <- subset(anlst,screen_group=="X-ray" & age>=70 & age<=74)
#obslyg <- rbind(obs.lyg(ct55_59,xray55_59),obs.lyg(ct60_64,xray60_64),obs.lyg(ct65_69,xray65_69),obs.lyg(ct70_74,xray70_74))

#jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_lyg_age2.jpeg",width=11.5,height=8,units='in',res=900)
#plot(seq(57,72,5),365*obslyg[,4],ylim=c(-50,65),type="l",xlab="Age at first CT screen",ylab="Days of life gained",main="Life Gained from Screening")
#lines(seq(57,72,5),365*obslyg[,3],col="red")
#dev.off()

#highlight 60 year old versus 74 year old in NLST (both has 3% risk of lung cancer after 5 years)
#ct[ct$pid==114068,]  3.25% risk of lung cancer, 74 year old female with prior cancer and heart attack, bmi of 43.5 - 13.4 days of life gained
#former smoker, quit 3 years ago, 1 pack per day, 57 years of smoking, no family history of lung cancer
#ct[ct$pid==203093,] 3.17% risk of lung cancer, 60 year old male with no comorbidities, bmi of 25.8 - 31.7 days of life gained
#current smoker, 1 pack per day, 40 years of smoking, with family history of lung cancer
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_illustration.jpeg",width=8,height=11.5,units='in',res=900)
par(mfrow=c(1,2))
c <- survfit(morat,newdata=ct[ct$pid==114068,],start.time=ct$age[ct$pid==114068])
c_alt <- survfit(morat,newdata=ct[ct$pid==114068,],start.time=ct$age[ct$pid==114068]+5,stop.time=100)
surv_5 <- c$surv[(c$time-ct$age[ct$pid==114068])==5]+((1-lcrisk(persons,5)[3]/1000)**.8-(1-lcrisk(persons,5)[3]/1000))

times <- seq(1,19,1)
ctsurv <- function(i){
  persons <- data.frame(age=ct$age[ct$pid==114068],
                        female=ct$female[ct$pid==114068],
                        smkyears=ct$smkyears[ct$pid==114068],
                        qtyears=ct$qtyears[ct$pid==114068],
                        cpd=ct$cpd[ct$pid==114068],
                        race=as.numeric(as.character(ct$race[ct$pid==114068])),
                        emp=ct$emp[ct$pid==114068],
                        fam.lung.trend=ct$fam.lung.trend[ct$pid==114068],
                        bmi=ct$bmi[ct$pid==114068],
                        edu6=ct$edu6[ct$pid==114068])
  if (i<=5){
     return(c$surv[(c$time-ct$age[ct$pid==114068])==i]+((1-lcrisk(persons,i)[3]/1000)**.8-(1-lcrisk(persons,i)[3]/1000)))
  } else {
    return(max(c$surv[(c$time-ct$age[ct$pid==114068])==i],as.numeric(surv_5*c_alt$surv[(c_alt$time-ct$age[ct$pid==114068])==i])))
  }      
}
ct_surv <- unlist(sapply(times,ctsurv))
c$surv[(c$time-ct$age[ct$pid==114068]) %in% times]



plot(times+ct$age[ct$pid==114068],c$surv[c$time %in% (times+ct$age[ct$pid==114068])],type="l",xlim=c(74,100),col="red",xlab="age",ylab="Survival")
lines(times+ct$age[ct$pid==114068],ct_surv,type="l")

sum(diff(times)*(ct_surv-c$surv[(c$time-ct$age[ct$pid==114068]) %in% times])[2:length(times)])


c <- survfit(morat,newdata=ct[ct$pid==203093,],start.time=ct$age[ct$pid==203093])
c_alt <- survfit(morat,newdata=ct[ct$pid==203093,],start.time=ct$age[ct$pid==203093]+5,stop.time=100)
surv_5 <- c$surv[(c$time-ct$age[ct$pid==203093])==5]+((1-lcrisk(persons,5)[3]/1000)**.8-(1-lcrisk(persons,5)[3]/1000))

times <- seq(1,33,1)
ctsurv <- function(i){
  persons <- data.frame(age=ct$age[ct$pid==203093],
                        female=ct$female[ct$pid==203093],
                        smkyears=ct$smkyears[ct$pid==203093],
                        qtyears=ct$qtyears[ct$pid==203093],
                        cpd=ct$cpd[ct$pid==203093],
                        race=as.numeric(as.character(ct$race[ct$pid==203093])),
                        emp=ct$emp[ct$pid==203093],
                        fam.lung.trend=ct$fam.lung.trend[ct$pid==203093],
                        bmi=ct$bmi[ct$pid==203093],
                        edu6=ct$edu6[ct$pid==203093])
  if (i<=5){
    return(c$surv[(c$time-ct$age[ct$pid==203093])==i]+((1-lcrisk(persons,i)[3]/1000)**.8-(1-lcrisk(persons,i)[3]/1000)))
  } else {
    return(max(c$surv[(c$time-ct$age[ct$pid==203093])==i],as.numeric(surv_5*c_alt$surv[(c_alt$time-ct$age[ct$pid==203093])==i])))
  }      
}
ct_surv <- unlist(sapply(times,ctsurv))
names(ct_surv) <- NULL

c$surv[(c$time-ct$age[ct$pid==203093]) %in% times]


plot(times+ct$age[ct$pid==203093],c$surv[c$time %in% (times+ct$age[ct$pid==203093])],type="l",xlim=c(60,100),col="red",xlab="age",ylab="Survival")
lines(times+ct$age[ct$pid==203093],ct_surv,type="l")

sum(diff(times)*(ct_surv-c$surv[(c$time-ct$age[ct$pid==203093]) %in% times])[2:length(times)])
dev.off()

anlst$comorbidities <- anlst$emp+anlst$hypertension+anlst$chd+anlst$angina+anlst$heartattack+
                       anlst$heartdisease+anlst$stroke+anlst$diab+anlst$bron+anlst$kidney+
                       anlst$liver+anlst$prior.cancer+anlst$speceq
anlst$comorbidities[anlst$comorbidities>=2] <- 2
xray0 <- subset(anlst,screen_group=="X-ray" & comorbidities==0)
ct0 <- subset(anlst,screen_group=="CT" & comorbidities==0)
xray1 <- subset(anlst,screen_group=="X-ray" & comorbidities==1)
ct1 <- subset(anlst,screen_group=="CT" & comorbidities==1)
xray2 <- subset(anlst,screen_group=="X-ray" & comorbidities>=2)
ct2 <- subset(anlst,screen_group=="CT" & comorbidities>=2)
obslyg <- rbind(obs.lyg(ct0,xray0),obs.lyg(ct1,xray1),obs.lyg(ct2,xray2))
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_lyg_comorbidities.jpeg",width=11.5,height=8,units='in',res=900)
p <- boxplot(365*lyg~comorbidities,data=anlst,ylim=c(-100,150),xlab="Number of comorbidities at first CT screen",ylab="Days of life gained",main="Life Gained from Screening")
lines(seq(1,3,1),365*obslyg[,3],col="red")
dev.off()

table(xray1$lcrat.cat)
table(ct1$lcrat.cat)

xray0_rq1 <- subset(anlst,screen_group=="X-ray" & comorbidities==0 & lcrat.cat==1)
ct0_rq1 <- subset(anlst,screen_group=="CT" & comorbidities==0 & lcrat.cat==1)
xray0_rq2 <- subset(anlst,screen_group=="X-ray" & comorbidities==0 & lcrat.cat==2)
ct0_rq2 <- subset(anlst,screen_group=="CT" & comorbidities==0 & lcrat.cat==2)
xray0_rq3 <- subset(anlst,screen_group=="X-ray" & comorbidities==0 & lcrat.cat==3)
ct0_rq3 <- subset(anlst,screen_group=="CT" & comorbidities==0 & lcrat.cat==3)
xray0_rq4 <- subset(anlst,screen_group=="X-ray" & comorbidities==0 & lcrat.cat==4)
ct0_rq4 <- subset(anlst,screen_group=="CT" & comorbidities==0 & lcrat.cat==4)
xray0_rq5 <- subset(anlst,screen_group=="X-ray" & comorbidities==0 & lcrat.cat==5)
ct0_rq5 <- subset(anlst,screen_group=="CT" & comorbidities==0 & lcrat.cat==5)
xray1_rq1 <- subset(anlst,screen_group=="X-ray" & comorbidities==1 & lcrat.cat==1)
ct1_rq1 <- subset(anlst,screen_group=="CT" & comorbidities==1 & lcrat.cat==1)
xray1_rq2 <- subset(anlst,screen_group=="X-ray" & comorbidities==1 & lcrat.cat==2)
ct1_rq2 <- subset(anlst,screen_group=="CT" & comorbidities==1 & lcrat.cat==2)
xray1_rq3 <- subset(anlst,screen_group=="X-ray" & comorbidities==1 & lcrat.cat==3)
ct1_rq3 <- subset(anlst,screen_group=="CT" & comorbidities==1 & lcrat.cat==3)
xray1_rq4 <- subset(anlst,screen_group=="X-ray" & comorbidities==1 & lcrat.cat==4)
ct1_rq4 <- subset(anlst,screen_group=="CT" & comorbidities==1 & lcrat.cat==4)
xray1_rq5 <- subset(anlst,screen_group=="X-ray" & comorbidities==1 & lcrat.cat==5)
ct1_rq5 <- subset(anlst,screen_group=="CT" & comorbidities==1 & lcrat.cat==5)
xray2_rq1 <- subset(anlst,screen_group=="X-ray" & comorbidities==2 & lcrat.cat==1)
ct2_rq1 <- subset(anlst,screen_group=="CT" & comorbidities==2 & lcrat.cat==1)
xray2_rq2 <- subset(anlst,screen_group=="X-ray" & comorbidities==2 & lcrat.cat==2)
ct2_rq2 <- subset(anlst,screen_group=="CT" & comorbidities==2 & lcrat.cat==2)
xray2_rq3 <- subset(anlst,screen_group=="X-ray" & comorbidities==2 & lcrat.cat==3)
ct2_rq3 <- subset(anlst,screen_group=="CT" & comorbidities==2 & lcrat.cat==3)
xray2_rq4 <- subset(anlst,screen_group=="X-ray" & comorbidities==2 & lcrat.cat==4)
ct2_rq4 <- subset(anlst,screen_group=="CT" & comorbidities==2 & lcrat.cat==4)
xray2_rq5 <- subset(anlst,screen_group=="X-ray" & comorbidities==2 & lcrat.cat==5)
ct2_rq5 <- subset(anlst,screen_group=="CT" & comorbidities==2 & lcrat.cat==5)
obslyg <- rbind(obs.lyg(ct0_rq1,xray0_rq1),obs.lyg(ct0_rq2,xray0_rq2),obs.lyg(ct0_rq3,xray0_rq3),obs.lyg(ct0_rq4,xray0_rq4),obs.lyg(ct0_rq5,xray0_rq5),
                obs.lyg(ct1_rq1,xray1_rq1),obs.lyg(ct1_rq2,xray1_rq2),obs.lyg(ct1_rq3,xray1_rq3),obs.lyg(ct1_rq4,xray1_rq4),obs.lyg(ct1_rq5,xray1_rq5),
                obs.lyg(ct2_rq1,xray2_rq1),obs.lyg(ct2_rq2,xray2_rq2),obs.lyg(ct2_rq3,xray2_rq3),obs.lyg(ct2_rq4,xray2_rq4),obs.lyg(ct2_rq5,xray2_rq5))

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_lyg_comorbidities_risk.jpeg",width=11.5,height=8,units='in',res=900)
plot(seq(1,5,1),365*obslyg[1:5,3],ylim=c(-150,150),ylab="life days gained",xlab="risk quintile",main="Life gained by risk and comorbidities")
points(seq(1,5,1),365*obslyg[6:10,3],col="blue")
points(seq(1,5,1),365*obslyg[11:15,3],col="red")
dev.off()

quant.cp <- quantile(anlst$lcrat,c(1/3,2/3))
anlst$lcrat.tert <- ifelse(anlst$lcrat<quant.cp[1],1,
                          ifelse(anlst$lcrat<quant.cp[2],2,3))


xray0_rq1 <- subset(anlst,screen_group=="X-ray" & comorbidities==0 & lcrat.tert==1)
ct0_rq1 <- subset(anlst,screen_group=="CT" & comorbidities==0 & lcrat.tert==1)
xray0_rq2 <- subset(anlst,screen_group=="X-ray" & comorbidities==0 & lcrat.tert==2)
ct0_rq2 <- subset(anlst,screen_group=="CT" & comorbidities==0 & lcrat.tert==2)
xray0_rq3 <- subset(anlst,screen_group=="X-ray" & comorbidities==0 & lcrat.tert==3)
ct0_rq3 <- subset(anlst,screen_group=="CT" & comorbidities==0 & lcrat.tert==3)
xray1_rq1 <- subset(anlst,screen_group=="X-ray" & comorbidities==1 & lcrat.tert==1)
ct1_rq1 <- subset(anlst,screen_group=="CT" & comorbidities==1 & lcrat.tert==1)
xray1_rq2 <- subset(anlst,screen_group=="X-ray" & comorbidities==1 & lcrat.tert==2)
ct1_rq2 <- subset(anlst,screen_group=="CT" & comorbidities==1 & lcrat.tert==2)
xray1_rq3 <- subset(anlst,screen_group=="X-ray" & comorbidities==1 & lcrat.tert==3)
ct1_rq3 <- subset(anlst,screen_group=="CT" & comorbidities==1 & lcrat.tert==3)
xray2_rq1 <- subset(anlst,screen_group=="X-ray" & comorbidities==2 & lcrat.tert==1)
ct2_rq1 <- subset(anlst,screen_group=="CT" & comorbidities==2 & lcrat.tert==1)
xray2_rq2 <- subset(anlst,screen_group=="X-ray" & comorbidities==2 & lcrat.tert==2)
ct2_rq2 <- subset(anlst,screen_group=="CT" & comorbidities==2 & lcrat.tert==2)
xray2_rq3 <- subset(anlst,screen_group=="X-ray" & comorbidities==2 & lcrat.tert==3)
ct2_rq3 <- subset(anlst,screen_group=="CT" & comorbidities==2 & lcrat.tert==3)
obslyg <- rbind(obs.lyg(ct0_rq1,xray0_rq1),obs.lyg(ct0_rq2,xray0_rq2),obs.lyg(ct0_rq3,xray0_rq3),
                obs.lyg(ct1_rq1,xray1_rq1),obs.lyg(ct1_rq2,xray1_rq2),obs.lyg(ct1_rq3,xray1_rq3),
                obs.lyg(ct2_rq1,xray2_rq1),obs.lyg(ct2_rq2,xray2_rq2),obs.lyg(ct2_rq3,xray2_rq3))

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_lyg_comorbidities_risktert.jpeg",width=11.5,height=8,units='in',res=900)
plot(seq(1,3,1),365*obslyg[1:3,3],ylim=c(-150,150),ylab="life days gained",xlab="risk tertile",main="Life gained by risk and comorbidities")
points(seq(1,3,1),365*obslyg[4:6,3],col="blue")
points(seq(1,3,1),365*obslyg[7:9,3],col="red")
dev.off()


anlst$comorb_score <- 0.48287*anlst$emp+0.13541*anlst$hypertension+0.18167*anlst$chd-0.10730*anlst$angina+0.27952*anlst$heartattack+
  0.15432*anlst$heartdisease+0.25903*anlst$stroke+0.45245*anlst$diab+0.11276*anlst$bron+0.40338*anlst$kidney+
  0.46008*anlst$liver+0.29890*anlst$prior.cancer+0.55969*anlst$speceq

quant.cp <- quantile(anlst$comorb_score[anlst$comorb_score>0],c(1/3,2/3))
anlst$comorb.cat <- ifelse(anlst$comorb_score<=0,0,
                           ifelse(anlst$comorb_score<=quant.cp[1],1,
                                  ifelse(anlst$comorb_score<=quant.cp[2],2,3)))

xray0 <- subset(anlst,screen_group=="X-ray" & comorb.cat==0)
ct0 <- subset(anlst,screen_group=="CT" & comorb.cat==0)
xray1 <- subset(anlst,screen_group=="X-ray" & comorb.cat==1)
ct1 <- subset(anlst,screen_group=="CT" & comorb.cat==1)
xray2 <- subset(anlst,screen_group=="X-ray" & comorb.cat==2)
ct2 <- subset(anlst,screen_group=="CT" & comorb.cat==2)
xray3 <- subset(anlst,screen_group=="X-ray" & comorb.cat==3)
ct3 <- subset(anlst,screen_group=="CT" & comorb.cat==3)
obslyg <- rbind(obs.lyg(ct0,xray0),obs.lyg(ct1,xray1),obs.lyg(ct2,xray2),obs.lyg(ct3,xray3))
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_lyg_comorb_score.jpeg",width=11.5,height=8,units='in',res=900)
p <- boxplot(365*lyg~comorb.cat,data=anlst,ylim=c(-100,150),xlab="Comorbidity score (categorical) at first CT screen",ylab="Days of life gained",main="Life Gained from Screening")
lines(seq(1,4,1),365*obslyg[,3],col="red")
dev.off()

xrayhypertension <- subset(anlst,screen_group=="X-ray" & hypertension==1)
cthypertension <- subset(anlst,screen_group=="CT" & hypertension==1)
365.25*obs.lyg(cthypertension,xrayhypertension)

anlst$lungcomorbidities <- anlst$emp+anlst$bron
anlst$comorbidities[anlst$lungcomorbidities>=3] <- 3
xray0 <- subset(anlst,screen_group=="X-ray" & lungcomorbidities==0)
ct0 <- subset(anlst,screen_group=="CT" & lungcomorbidities==0)
xray1 <- subset(anlst,screen_group=="X-ray" & lungcomorbidities==1)
ct1 <- subset(anlst,screen_group=="CT" & lungcomorbidities==1)
xray2 <- subset(anlst,screen_group=="X-ray" & lungcomorbidities>=2)
ct2 <- subset(anlst,screen_group=="CT" & lungcomorbidities>=2)
obslyg <- rbind(obs.lyg(ct0,xray0),obs.lyg(ct1,xray1),obs.lyg(ct2,xray2))
jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_lyg_lungcomorbidities.jpeg",width=11.5,height=8,units='in',res=900)
p <- boxplot(365*lyg~lungcomorbidities,data=anlst,ylim=c(-100,150),xlab="Number of lung comorbidities at first CT screen",ylab="Days of life gained",main="Life Gained from Screening")
lines(seq(1,3,1),365*obslyg[,3],col="red")
dev.off()

kmm_ct1 <- survfit(Surv(ct1$death.years,ct1$death)~1)
kmm_xray1 <- survfit(Surv(xray1$death.years,xray1$death)~1)

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nlst_lungcomorbidities_surv.jpeg",width=11.5,height=8,units='in',res=900)
plot(kmm_ct1,ylim=c(.8,1))
lines(kmm_xray1,col="red")
dev.off()

morat <- coxph(Surv(age, age+death.years, death) ~ 
                    female + race + edu6 + rand_year + 
                    I(bmi <= 18.5) + I(bmi > 18.5 & bmi <= 20) + I(bmi > 25 & bmi <= 30) + I(bmi > 30 & bmi <= 35) + I(bmi > 35) + 
                    emp + hypertension + heartattack + stroke + diab + bron + prior.cancer + 
                    I(log(qtyears + 1)) + I(log(cpd + 1)) + I(sqrt(pkyr.cat)),
                    data = xray, model = TRUE, y = TRUE)

ely_est <- function(x){
  c <- survfit(morat,newdata=z[x,],start.time=z$age[x])
  ely <- sum(diff(c(z$age[x],c$time))*c$surv)  
  return(ely)
}

z <- ct
ely_res <- sapply(1:nrow(z),ely_est)

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/scatter_risks_lyg.jpeg",width=11.5,height=8,units='in',res=900)
ggplot(anlst,aes(x=100*lcrat,y=365*lyg)) + geom_point(alpha = 0.05) + 
  labs(x="% lung cancer risk") + 
  labs(y="Days of life gained") + 
  coord_cartesian(xlim = c(0, 10),ylim = c(0,100)) +
  geom_hline(yintercept=16.8) +
  annotate("text",x=9,y=20,label="Average gain in NLST=16.8 days")
dev.off()
ct$ely_res <- ely_res
save(ct,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/ct.Rdata")
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/ct.Rdata")
ct$expecteddayslived <- pmin(365.25*ct$ely_res,ct$days_rand_studyend)
ct$dayslived_studyend <- ifelse(ct$death==1,365.25*ct$death.years,ct$days_rand_studyend)
ct$actual_life_gained <- ct$dayslived_studyend-ct$expecteddayslived

ctlived <- subset(ct,death==0)
ctdeath <- subset(ct,death==1)

p2 <- ggplot(ctlived, aes(x=lcrat.cat,y=actual_life_gained)) +
  geom_boxplot(outlier.size=.5,aes(group=cut_width(lcrat.cat,1))) +
  ggtitle("Life Gained from Screening in CT arm (Alive by end of study)") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(x = "Risk Quintiles (lowest to highest)", y = "Days of life gained")

p3 <- ggplot(ctdeath, aes(x=lcrat.cat,y=actual_life_gained)) +
  geom_boxplot(outlier.size=.5,aes(group=cut_width(lcrat.cat,1))) +
  ggtitle("Life Gained from Screening in CT arm (Dead by end of study)") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(x = "Risk Quintiles (lowest to highest)", y = "Days of life gained")

p4 <- ggplot(ct, aes(x=lcrat.cat,y=actual_life_gained)) +
  geom_boxplot(outlier.size=.5,aes(group=cut_width(lcrat.cat,1))) +
  ggtitle("Life Gained from Screening in CT arm") +
  theme(plot.title = element_text(hjust=0.5)) +
  labs(x = "Risk Quintiles (lowest to highest)", y = "Days of life gained")

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/model_assisted_obs_lyg_alive.jpeg",width=11.5,height=8,units='in',res=900)
plot(p2)
dev.off()

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/model_assisted_obs_lyg_dead.jpeg",width=11.5,height=8,units='in',res=900)
plot(p3)
dev.off()

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/model_assisted_obs_lyg_all.jpeg",width=11.5,height=8,units='in',res=900)
plot(p4)
dev.off()


ct55_59 <- subset(ct,age<=59)
ct60_64 <- subset(ct,age>=60 & age<=64)
ct65_69 <- subset(ct,age>=65 & age<=69)
ct70_74 <- subset(ct,age>=70 & age<=74)
mean(ct55_59$actual_life_gained)
mean(ct60_64$actual_life_gained)
mean(ct65_69$actual_life_gained)
mean(ct70_74$actual_life_gained)
#model-assisted actual life gained in -34, -46, -68, -109 for ages 55-59, 60-64, 65-69, 70-74
