rm(list=ls())
library(survival)
library(survey)
library(lcrisks)
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis.imputed")


nhis <- subset(nhis,deathage>age)
master  <- svydesign(id=~psu, strata=~strata, weights=~wt_mort5, data=nhis, nest=TRUE)
design <- subset(master, age>=40)
morat <- svycoxph(Surv(age, deathage, died) ~ 
                  female+race+edu6+I(bmi <= 18.5)+cpd+pkyr.cat+bmi+qtyears+smkyears+birthyear+
                  emp+hypertension+chd+angina+heartattack+heartdisease+stroke+diab+bron+kidney+liver+speceq+prior.cancer,
                  design=design, data = nhis, model = TRUE, y = TRUE) 


load(file="~/Desktop/Lung cancer/lrisk/other/nhis2010_12/nhis_impute_1.v3.RData")
persons <- data.frame(age=nhis$age,
                      female=nhis$female,
                      smkyears=nhis$smkyears,
                      qtyears=nhis$qtyears,
                      cpd=nhis$cpd,
                      race=as.numeric(as.character(nhis$race)),
                      emp=nhis$emp,
                      fam.lung.trend=nhis$fam.lung.trend,
                      bmi=nhis$bmi,
                      edu6=nhis$edu6)
nhis$lcd1 <- lcrisk(persons,1)[,3]/1000
nhis$lcd2 <- lcrisk(persons,2)[,3]/1000   
nhis$lcd3 <- lcrisk(persons,3)[,3]/1000   
nhis$lcd4 <- lcrisk(persons,4)[,3]/1000   
nhis$lcd5 <- lcrisk(persons,5)[,3]/1000   

nhis$birthyear <- nhis$year-nhis$age

lyg <- function(x){
  print(x)
  c <- survfit(morat,newdata=nhis[x,],start.time=nhis$age[x])
  ely <- sum(diff(c(nhis$age[x],c$time))*c$surv)  
  c <- survfit(morat,newdata=nhis[x,],start.time=nhis$age[x]+5,stop.time=100)
  lyg <- (1-nhis$lcd1[x])**.8-(1-nhis$lcd1[x]) + 
         (1-nhis$lcd2[x])**.8-(1-nhis$lcd2[x]) +
         (1-nhis$lcd3[x])**.8-(1-nhis$lcd3[x]) +
         (1-nhis$lcd4[x])**.8-(1-nhis$lcd4[x]) +
         ((1-nhis$lcd5[x])**.8-(1-nhis$lcd5[x]))*(1+sum(diff(c(nhis$age[x]+5,c$time))*c$surv))
  res <- c(ely,ely+lyg,lyg)
  return(res)
}

lyg_res <- sapply(1:nrow(nhis),lyg)
lyg_res <- t(lyg_res)
nhis$ely <- lyg_res[,1]
nhis$ely_ct <- lyg_res[,2]
nhis$lyg <- lyg_res[,3]
save(nhis,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj")

load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_proj")
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
#overall
quantiles <- c(0,.25,.5,.75,1)
svyquantile(~lyg,master,quantiles)*365.25
svyquantile(~ely,master,quantiles)

# stratify everything by lung cancer risks
quintiles <- c(.2,.4,.6,.8)
svyquantile(~lcd5,master,quintiles)
design <- subset(master,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,lcd5>0.00421&lcd5<=0.0110)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,lcd5>0.0110)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)


subcohort <- subset(master,female==0)
svyquantile(~lyg,subcohort,quantiles)*365.25
svyquantile(~ely,subcohort,quantiles)
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,female==1)
svyquantile(~lyg,subcohort,quantiles)*365.25
svyquantile(~ely,subcohort,quantiles)
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)


svytable(~current,master)/62027579
design <- subset(master,current==1)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,current==1)
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

svytable(~former,master)/62027579
design <- subset(master,former==1)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,former==1)  #need to change everything from here on to former
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)






cumsum(svytable(~age,master)/62027579)
design <- subset(master,age<50)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,age>=50&age<60)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,age>=60&age<70)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,age>=70&age<80)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,age>=80)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,age<50)  
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,age>=50 & age<60) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,age>=60 & age<70) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,age>=70 & age<80) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,age>=80) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

quintiles <- c(0,.2,.4,.6,.8,1)
svyquantile(~ely,master,quintiles)
design <- subset(master,ely<=12)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,ely>12&ely<=17)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,ely>17&ely<=21.5)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,ely>21.5&ely<=26)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,ely>26)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,ely<=12) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,ely>12 & ely<=17) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,ely>17 & ely<=21.5) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,ely>21.5 & ely<=26) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,ely>26) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
                      nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
                      nhis$liver+nhis$prior.cancer

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
svytable(~comorbidities,master)/62027579
design <- subset(master,comorbidities==0)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,comorbidities==1)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,comorbidities==2)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(master,comorbidities>=3)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,comorbidities==0) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,comorbidities==1) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,comorbidities==2) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)

subcohort <- subset(master,comorbidities>=3) 
design <- subset(subcohort,lcd5<=0.0007)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0007&lcd5<=0.0018)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0018&lcd5<=0.0042)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.0042&lcd5<=0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)
design <- subset(subcohort,lcd5>0.011)
svyquantile(~lyg,design,quantiles)*365.25
svyquantile(~ely,design,quantiles)