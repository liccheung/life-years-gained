rm(list=ls(all=TRUE))
library(survival)
library(survey)
load(file="/home/cheunglc/lyg/mortality.model.v2.RData")

drest <- function(x){
  c <- survfit(morat,newdata=z[x,],start.time=z$age[x])
  yr_flwup <- 2012-z$intyear[x]  #death data until end of 2011
  if (max(c$surv[c$time>=(z$age[1]+yr_flwup)])==min(c$surv[c$time<=(z$age[1]+yr_flwup)])){
    death_risk <- 1-c$surv[c$time==(z$age[1]+yr_flwup)]
  } else {
    m <- (max(c$surv[c$time>=(z$age[1]+yr_flwup)])-min(c$surv[c$time<=(z$age[1]+yr_flwup)]))/(max(c$time[c$time>=(z$age[1]+yr_flwup)])-min(c$time[c$time<=(z$age[1]+yr_flwup)]))
    b <- max(c$surv[c$time>=(z$age[1]+yr_flwup)])-m*max(c$time[c$time>=(z$age[1]+yr_flwup)])
    death_risk <- 1-(m*(z$age[1]+yr_flwup)+b)
  }
  return(death_risk)
}

load(file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_1.RData")
nhis$birthyear <- nhis$year-nhis$age
z <- nhis
z$death_risk <- sapply(1:nrow(z),drest)
save(file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_1_mod.RData")

load(file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_2.RData")
nhis$birthyear <- nhis$year-nhis$age
z <- nhis
z$death_risk <- sapply(1:nrow(z),drest)
save(file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_2_mod.RData")

load(file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_3.RData")
nhis$birthyear <- nhis$year-nhis$age
z <- nhis
z$death_risk <- sapply(1:nrow(z),drest)
save(file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_3_mod.RData")

load(file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_4.RData")
nhis$birthyear <- nhis$year-nhis$age
z <- nhis
z$death_risk <- sapply(1:nrow(z),drest)
save(file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_4_mod.RData")

load(file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_5.RData")
nhis$birthyear <- nhis$year-nhis$age
z <- nhis
z$death_risk <- sapply(1:nrow(z),drest)
save(file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_5_mod.RData")
