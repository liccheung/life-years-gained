#impute smoking variables - only in analysis population
# EXCLUSIONS OF NEVER/UNKNOWN SMOKERS
# PRIOR LUNG CANCER
# MISSING/NEGATIVE SMOKE AGE START
# ADULTS NOT 40+ YEARS OF AGE
nhis$analypop <- ifelse(nhis$never==0 & nhis$unknown==0 & 
                          !(is.na(nhis$smoke.age.start)|nhis$smoke.age.start<0),1,0)

notanaly <- subset(nhis,analypop==0)
nhis <- subset(nhis,analypop==1)  #impute variables for analysis population

nhis$age.cat <- cut(nhis$age, breaks=quantile(nhis$age,seq(0,1,by=.1))+c(-1,rep(0,9),1)) #age deciles
nhis$age.quartile <- cut(nhis$age, breaks=quantile(nhis$age,seq(0,1,by=.25))+c(-1,rep(0,3),1),labels=1:4) #age quartiles

#death probability is different for each nhis year as earlier years have more follow-up
nhis$group <- interaction(nhis$age.cat, nhis$female, nhis$year, nhis$died, drop=TRUE) 

impute.function <- function(strata, variable){
  values <- nhis[,variable]
  values <- values[!is.na(values)&strata==nhis$group]
  sample(values,size=1)
}

#impute 832 with missing race - random assigment with same proportions as non-missing
race.impute <- sapply(nhis$group[is.na(nhis$race)],impute.function,variable="race")
nhis$race[is.na(nhis$race)] <- race.impute

#death probability is different for each nhis year as earlier years have more follow-up
nhis$group <- interaction(nhis$age.cat, nhis$female, nhis$year, nhis$died, drop=TRUE) 
cpd.factor <- function(female, age){
  cpddiff$diff[cpddiff$age==age&cpddiff$female==female]
}

impute.function <- function(strata, case, year, variable){
  values <- nhis[,variable]
  values <- values[!is.na(values)&strata==nhis$group]
  
  #if no one is in that group, then sample from case group from that year
  if(length(values)==0){
    values <- nhis[,variable]
    values <- values[!is.na(values)&case==nhis$died & !is.na(values)&year==nhis$year]
  }
  sample(values,size=1)
}

#impute 33013 with missing cig per day
cpd.impute <- mapply(impute.function,strata=nhis$group[is.na(nhis$cpd)],
                     case = nhis$died[is.na(nhis$cpd)],
                     year = nhis$year[is.na(nhis$cpd)],
                     MoreArgs=list(variable="cpd"))

cpd.impute.add <- mapply(cpd.factor, age=nhis$age.quartile[is.na(nhis$cpd)], female=nhis$female[is.na(nhis$cpd)])

nhis$cpd[is.na(nhis$cpd)] <- ifelse(nhis$former[is.na(nhis$cpd)]==1,cpd.impute+cpd.impute.add,cpd.impute)
nhis$cpd[nhis$cpd<5] <- 5

impute.function <- function(strata, variable){
  values <- nhis[,variable]
  values <- values[!is.na(values)&strata==nhis$group&nhis$former==1]
  sample(values,size=1)
}

qtyears.impute <- sapply(nhis$group[is.na(nhis$qtyears)],impute.function,variable="qtyears")
nhis$qtyears[is.na(nhis$qtyears)] <- qtyears.impute

drops <- c("age.cat","age.quartile","group")
nhis <- nhis[ , !(names(nhis) %in% drops)]
nhis <- rbind(nhis,notanaly)
nhis$race <- as.factor(nhis$race)


nhis$age.cat <- cut(nhis$age, breaks=quantile(nhis$age,seq(0,1,by=.1))+c(-1,rep(0,9),1)) #age deciles
nhis$age.quartile <- cut(nhis$age, breaks=quantile(nhis$age,seq(0,1,by=.25))+c(-1,rep(0,3),1),labels=1:4) #age quartiles
#death probability is different for each nhis year as earlier years have more follow-up
nhis$group <- interaction(nhis$age.cat, nhis$female, nhis$year, nhis$died, drop=TRUE) 

impute.function <- function(strata, variable){
  values <- nhis[,variable]
  values <- values[!is.na(values)&strata==nhis$group]
  sample(values,size=1)
}

#1521 missing bmi - random assigment with same proportions as non-missing (by age cat, gender, and mortality)
bmi.impute <- sapply(nhis$group[is.na(nhis$bmi)],impute.function,variable="bmi")
nhis$bmi[is.na(nhis$bmi)] <- bmi.impute

#impute 287 missing education - random assigment with same proportions as non-missing (by age cat and gender)
nhis$edu6 <- as.numeric(nhis$edu6)
nhis$edu6[nhis$edu6==max(nhis$edu6)] <- NA
edu6.impute <- sapply(nhis$group[is.na(nhis$edu6)],impute.function,variable="edu6")
nhis$edu6[is.na(nhis$edu6)] <- edu6.impute

# PREPARE VARIABLES
nhis$bmi[nhis$bmi>40] <- 40
nhis$bmi[nhis$bmi<15] <- 15

nhis$black <- ifelse(nhis$race=="Non-Hispanic Black",1,0)
nhis$hispanic <- ifelse(nhis$race=="Hispanic",1,0)
nhis$asian <- ifelse(nhis$race=="Asian",1,0)
nhis$islander <- ifelse(nhis$race=="Pacific Islander",1,0)
nhis$indian <- ifelse(nhis$race=="American Indian/Alaskan Native",1,0)
nhis$qt.trend <- ifelse(nhis$qtyears<=1,0,ifelse(nhis$qtyears<=5,1,2))
nhis$race <- factor(ifelse(nhis$black==1,1,
						   ifelse(nhis$hispanic==1,2,
							ifelse(nhis$race=="Non-Hispanic White",0,3)))	)						
							

# PACKYEARS
nhis$smkyears[is.na(nhis$smkyears)] <- (nhis$age-nhis$smoke.age.start-nhis$qtyears)[is.na(nhis$smkyears)] 
nhis$smkyears[nhis$smkyears<=0] <- 1
nhis$packyears <- nhis$smkyears*(nhis$cpd/20)
nhis$packyears[nhis$packyears>90] <- 90
nhis$pkyr.cat <- nhis$packyears

nhis$heart <- ifelse(nhis$heartattack==1|nhis$heartdisease==1|nhis$hypertension==1|nhis$chd==1|nhis$stroke==1|nhis$angina==1,1,0)
nhis$resp <- ifelse(nhis$emp==1|nhis$bron==1,1,0)