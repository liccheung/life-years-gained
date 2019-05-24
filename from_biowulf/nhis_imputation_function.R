#impute smoking variables - only in analysis population
# EXCLUSIONS OF NEVER/UNKNOWN SMOKERS
# PRIOR LUNG CANCER
# MISSING/NEGATIVE SMOKE AGE START
# ADULTS NOT 40+ YEARS OF AGE
nhis$analypop <- ifelse(nhis$never==0 & nhis$unknown==0 & nhis$lung.cancer.before==0 & 
                          !(is.na(nhis$smoke.age.start)|nhis$smoke.age.start<0),1,0)

notanaly <- subset(nhis,analypop==0)
nhis <- subset(nhis,analypop==1)  #impute variables for analysis population

# c(-1,rep(0,9),1) creates numeric repeating string -1 0 0 0 0 0 0 0 0 0 1
# seq(0,1,by=.1)   creates numeric string 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
nhis$age.cat      <- cut(nhis$age, breaks=quantile(nhis$age,seq(0,1,by=.1))+c(-1,rep(0,9),1))

nhis$age.quartile <- cut(nhis$age, breaks=quantile(nhis$age,seq(0,1,by=.25))+c(-1,rep(0,3),1),labels=1:4)

nhis$group <- interaction(nhis$age.cat, nhis$female,drop=TRUE) # AGE AND GENDER-SPECIFIC GROUP

impute.function <- function(strata, variable){
	values <- nhis[,variable]
	values <- values[!is.na(values)&strata==nhis$group&nhis$former==1]
sample(values,size=1)
}

qtyears.impute <- sapply(nhis$group[is.na(nhis$qtyears)],impute.function,variable="qtyears")
nhis$qtyears[is.na(nhis$qtyears)] <- qtyears.impute
nhis$agequit[is.na(nhis$agequit)] <- nhis$age[is.na(nhis$agequit)]-qtyears.impute


nhis$group <- interaction(nhis$group, nhis$former)

impute.function <- function(strata, variable){
	values <- nhis[,variable]
	values <- values[!is.na(values)&strata==nhis$group&nhis$year==2005]
sample(values,size=1)
}

cpd.impute <- sapply(nhis$group[is.na(nhis$cpd)],impute.function,variable="cpd")
nhis$cpd[is.na(nhis$cpd)] <- cpd.impute

drops <- c("age.cat","age.quartile","group")
nhis <- nhis[ , !(names(nhis) %in% drops)]
nhis <- rbind(nhis,notanaly)


#Impute on entire population
# c(-1,rep(0,9),1) creates numeric repeating string -1 0 0 0 0 0 0 0 0 0 1
# seq(0,1,by=.1)   creates numeric string 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
nhis$age.cat      <- cut(nhis$age, breaks=quantile(nhis$age,seq(0,1,by=.1))+c(-1,rep(0,9),1))

nhis$age.quartile <- cut(nhis$age, breaks=quantile(nhis$age,seq(0,1,by=.25))+c(-1,rep(0,3),1),labels=1:4)

nhis$group <- interaction(nhis$age.cat, nhis$female,drop=TRUE) # AGE AND GENDER-SPECIFIC GROUP

impute.function <- function(strata, variable){
  values <- nhis[,variable]
  values <- values[!is.na(values)&strata==nhis$group]
  sample(values,size=1)
}

race.impute <- sapply(nhis$group[is.na(nhis$race)],impute.function,variable="race")
nhis$race[is.na(nhis$race)] <- race.impute

bmi.impute <- sapply(nhis$group[is.na(nhis$bmi)],impute.function,variable="bmi")
nhis$bmi[is.na(nhis$bmi)] <- bmi.impute

fam.lung.trend.impute <- sapply(nhis$group[nhis$year!=2005],impute.function,variable="fam.lung.trend")
nhis$fam.lung.trend[nhis$year!=2005] <- fam.lung.trend.impute 

# Set smoking cancer to lung cancer
nhis$fam.smoke.cancer <- ifelse(nhis$fam.lung.trend>=1,1,0)
# Assume family lung cancer are all LATE onset. 0=None, 1=early, 2=late
nhis$fam.cancer.onset <- ifelse(nhis$fam.lung.trend>=1,2,0)

nhis$group <- interaction(nhis$group, nhis$fam.smoke.cancer)
fam.cancer.impute <- sapply(nhis$group[nhis$year!=2005],impute.function,variable="fam.cancer")
nhis$fam.cancer[nhis$year!=2005] <- fam.cancer.impute

nhis$copd[is.na(nhis$copd)] <- ifelse(nhis$emp[is.na(nhis$copd)]==1 | nhis$bron[is.na(nhis$copd)]==1,1,0)

# PREPARE VARIABLES
nhis$bmi[nhis$bmi>40] <- 40
nhis$bmi[nhis$bmi<15] <- 15

nhis$black <- ifelse(nhis$race=="Non-Hispanic Black",1,0)
nhis$hispanic <- ifelse(nhis$race=="Hispanic",1,0)
nhis$asian <- ifelse(nhis$race=="Asian",1,0)
nhis$islander <- ifelse(nhis$race=="Pacific Islander",1,0)
nhis$indian <- ifelse(nhis$race=="American Indian/Alaskan Native",1,0)
nhis$qt.trend <- ifelse(nhis$qtyears<1,0,ifelse(nhis$qtyears<=5,1,2))
nhis$race <- factor(ifelse(nhis$black==1,1,
						   ifelse(nhis$hispanic==1,2,
							ifelse(nhis$race=="Non-Hispanic White",0,3)))	)						
							
# PACKYEARS
nhis$smkyears[is.na(nhis$smkyears)] <- (nhis$age-nhis$smoke.age.start-nhis$qtyears)[is.na(nhis$smkyears)] 
nhis$smkyears[nhis$smkyears<=0] <- 1
nhis$packyears <- nhis$smkyears*(nhis$cpd/20)
nhis$packyears[nhis$packyears>90] <- 90
nhis$pkyr.cat <- nhis$packyears

nhis$medicare.eligible <- nhis$age>=55&nhis$age<75&nhis$packyears>=30
nhis$medicare.eligible[nhis$former==1] <- nhis$qtyears[nhis$former==1]<=15&nhis$medicare.eligible[nhis$former==1] 

nhis$uspstf.eligible <- nhis$age>=55&nhis$age<=80&nhis$packyears>=30
nhis$uspstf.eligible[nhis$former==1] <- nhis$qtyears[nhis$former==1]<=15&nhis$uspstf.eligible[nhis$former==1]