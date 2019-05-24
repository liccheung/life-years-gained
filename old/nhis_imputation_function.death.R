nhis$age.cat <- cut(nhis$age, breaks=quantile(nhis$age,seq(0,1,by=.1))+c(-1,rep(0,9),1)) #age deciles
nhis$age.quartile <- cut(nhis$age, breaks=quantile(nhis$age,seq(0,1,by=.25))+c(-1,rep(0,3),1),labels=1:4) #age quartiles
nhis$group <- interaction(nhis$age.cat, nhis$female, nhis$lung.cancer.death,drop=TRUE)

impute.function <- function(strata, variable){
	values <- nhis[,variable]
	values <- values[!is.na(values)&strata==nhis$group]
sample(values,size=1)
}

cpd.factor <- function(female, age){
	cpddiff$diff[cpddiff$age==age&cpddiff$female==female]
}

#impute 850 with missing race - random assigment with same proportions
nhis$race[nhis$race=="Missing"] <- sample(nhis$race[nhis$race!="Missing"],size=sum(nhis$race=="Missing"))
race.impute <- sapply(nhis$group[is.na(nhis$race)],impute.function,variable="race")

#1605 missing bmi
bmi.impute <- sapply(nhis$group[is.na(nhis$bmi)],impute.function,variable="bmi")
nhis$bmi[is.na(nhis$bmi)] <- bmi.impute

#impute 306 missing education
nhis$edu6 <- as.numeric(nhis$edu6)
nhis$edu6[nhis$edu6==max(nhis$edu6)] <- NA
edu6.impute <- sapply(nhis$group[is.na(nhis$edu6)],impute.function,variable="edu6")
nhis$edu6[is.na(nhis$edu6)] <- edu6.impute

impute.function <- function(strata, case, variable){
	values <- nhis[,variable]
	values <- values[!is.na(values)&strata==nhis$group]
	if(length(values)==0){
		values <- nhis[,variable]
		values <- values[!is.na(values)&case==nhis$lung.cancer.death]
		}
sample(values,size=1)
}

#impute 34674 with missing cig per day
cpd.impute <- mapply(impute.function,strata=nhis$group[is.na(nhis$cpd)],
											case = nhis$lung.cancer.death[is.na(nhis$cpd)],
											MoreArgs=list(variable="cpd"))

cpd.impute.add <- mapply(cpd.factor, age=nhis$age.quartile[is.na(nhis$cpd)], female=nhis$female[is.na(nhis$cpd)])
											
nhis$cpd[is.na(nhis$cpd)] <- cpd.impute+cpd.impute.add
nhis$cpd[nhis$cpd<5] <- 5

impute.function <- function(strata, variable){
	values <- nhis[,variable]
	values <- values[!is.na(values)&strata==nhis$group&nhis$former==1]
sample(values,size=1)
}

qtyears.impute <- sapply(nhis$group[is.na(nhis$qtyears)],impute.function,variable="qtyears")
nhis$qtyears[is.na(nhis$qtyears)] <- qtyears.impute

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
							
nhis$pneu <- 0

# PACKYEARS
nhis$smkyears[is.na(nhis$smkyears)] <- (nhis$age-nhis$smoke.age.start-nhis$qtyears)[is.na(nhis$smkyears)] 
nhis$smkyears[nhis$smkyears<=0] <- 1
nhis$packyears <- nhis$smkyears*(nhis$cpd/20)
nhis$packyears[nhis$packyears>90] <- 90
nhis$pkyr.cat <- nhis$packyears
nhis$age.stopped <- ifelse(nhis$former==1,nhis$age-nhis$qtyears,0)
nhis$age.stopped[nhis$smkyears==0] <- nhis$smoke.age.start[nhis$smkyears==0] 
nhis$uspstf.eligible <- nhis$age>=55&nhis$age<=80&nhis$packyears>=30
nhis$uspstf.eligible[nhis$former==1] <- nhis$qtyears[nhis$former==1]<=15&nhis$uspstf.eligible[nhis$former==1]
nhis$medicare.eligible <- nhis$age>=55&nhis$age<75&nhis$packyears>=30
nhis$medicare.eligible[nhis$former==1] <- nhis$qtyears[nhis$former==1]<=15&nhis$medicare.eligible[nhis$former==1] 



