rm(list=ls(all=TRUE))
library(survival)
library(survey)
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/mortality.model.v2.RData")

drest <- function(x){
  print(x)
  c <- survfit(morat,newdata=z[x,],start.time=z$age[x])
  yr_flwup <- 2012-z$intyear[x]  #death data until end of 2011
  if (max(c$surv[c$time>=(z$age[x]+yr_flwup)])==min(c$surv[c$time<=(z$age[x]+yr_flwup)])){
    death_risk <- 1-c$surv[c$time==(z$age[x]+yr_flwup)]
  } else {
    #impute midpoint between estimated risks
    m <- (max(c$surv[c$time>=(z$age[x]+yr_flwup)])-min(c$surv[c$time<=(z$age[x]+yr_flwup)]))/(max(c$time[c$time>=(z$age[x]+yr_flwup)])-min(c$time[c$time<=(z$age[x]+yr_flwup)]))
    b <- max(c$surv[c$time>=(z$age[x]+yr_flwup)])-m*max(c$time[c$time>=(z$age[x]+yr_flwup)])
    death_risk <- 1-(m*(z$age[x]+yr_flwup)+b)
  }
  return(death_risk)
}

load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_1.RData")
nhis$birthyear <- nhis$year-nhis$age
notanaly <- subset(nhis,analypop==0|age<40|age>80)
z <- subset(nhis,analypop==1&age>=40&age<=80)  #impute variables for analysis population
z$death_risk <- sapply(1:nrow(z),drest)
nhis <- rbind(z,notanaly)
save(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_1_mod.RData")

load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_2.RData")
nhis$birthyear <- nhis$year-nhis$age
z <- nhis
z$death_risk <- sapply(1:nrow(z),drest)
save(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_2_mod.RData")

load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_3.RData")
nhis$birthyear <- nhis$year-nhis$age
z <- nhis
z$death_risk <- sapply(1:nrow(z),drest)
save(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_3_mod.RData")

load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_4.RData")
nhis$birthyear <- nhis$year-nhis$age
z <- nhis
z$death_risk <- sapply(1:nrow(z),drest)
save(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_4_mod.RData")

load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_5.RData")
nhis$birthyear <- nhis$year-nhis$age
z <- nhis
z$death_risk <- sapply(1:nrow(z),drest)
save(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_5_mod.RData")