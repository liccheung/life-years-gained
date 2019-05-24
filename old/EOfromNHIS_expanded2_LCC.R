# PERFORMANCE
library(survey)
library(gdata)
setwd("P:/bb/lrisk/other/lifeyearsgained/nhis")


#calculate observed, expected, E/O, Var(E/O) for the following populations:
#overall; USPSTF eligible; Ineligible, 50-80; All Ineligible; All 50-80; Other Ineligible
popEO <- function(design){

   #observed, expected, E/O, Var(E/O) for all subjects
	O <- svytotal(~five.year, design)  #observed
	E <- svytotal(~risk, design)  #expected
	EO.hat <- svyratio(num=~risk, denom=~five.year, design)  #E/O
	overall <- cbind(coef(O), coef(E), coef(EO.hat), vcov(EO.hat))  #summary statistics
	
   #observed, expected, E/O, Var(E/O) for each group
	eligible.design <- subset(design, eligible.group==1)
	O <- svytotal(~five.year, eligible.design)
	E <- svytotal(~risk, eligible.design)
	EO.hat <- svyratio(num=~risk, denom=~five.year, eligible.design)
	eligible1 <- cbind(coef(O), coef(E), coef(EO.hat), vcov(EO.hat))

	ineligible.design <- subset(design, eligible.group==2)
	O <- svytotal(~five.year, ineligible.design)
	E <- svytotal(~risk, ineligible.design)
	EO.hat <- svyratio(num=~risk, denom=~five.year, ineligible.design)
	eligible2 <- cbind(coef(O), coef(E), coef(EO.hat), vcov(EO.hat))

	ineligible.design <- subset(design, eligible.group!=1)
	O <- svytotal(~five.year, ineligible.design)
	E <- svytotal(~risk, ineligible.design)
	EO.hat <- svyratio(num=~risk, denom=~five.year, ineligible.design)
	eligible3 <- cbind(coef(O), coef(E), coef(EO.hat), vcov(EO.hat))

	ineligible.design <- subset(design, eligible.group!=3)
	O <- svytotal(~five.year, ineligible.design)
	E <- svytotal(~risk, ineligible.design)
	EO.hat <- svyratio(num=~risk, denom=~five.year, ineligible.design)
	eligible4 <- cbind(coef(O), coef(E), coef(EO.hat), vcov(EO.hat))

	ineligible.design <- subset(design, eligible.group==3)
	O <- svytotal(~five.year, ineligible.design)
	E <- svytotal(~risk, ineligible.design)
	EO.hat <- svyratio(num=~risk, denom=~five.year, ineligible.design)
	eligible5 <- cbind(coef(O), coef(E), coef(EO.hat), vcov(EO.hat))


result <- rbind(overall, eligible1, eligible2, eligible3, eligible4, eligible5)

row.names(result) <- c("All Smokers",
								"USPSTF Eligible",
								"Ineligible, 50-80",
								"All Ineligible",
								"All 50-80",
								"Other Ineligible")

result
}

#For each categorical grouping, calculate expected, observed, E/O, Var(E/O)
popriskEO <- function(factor, design){

	f <- function(group){
		
			keep <- design$variables[,factor]==group
			subset.design <- subset(design, keep)
			O <- svytotal(~five.year, subset.design)
			E <- svytotal(~risk, subset.design)
			EO.hat <- svyratio(num=~risk, denom=~five.year, subset.design)
			
		cbind(coef(E), coef(O), coef(EO.hat), vcov(EO.hat))
		}

	design$variables[,factor] <- drop.levels(design$variables[,factor], reorder=F)  #drop any unused levels
	thelevels <- levels(design$variables[,factor])

lapply(thelevels, f)
}

#output weighted quantile of Y
pop.quantile <- function(y, wt, breaks){

   #calculate empirical cdf of Y 
	freq <- tapply(wt, y, sum)
	freq <- cumsum(freq/sum(freq))	

	index <- sapply(breaks, function(x) {
			if(x==0)
				1
			else if (x==1)
				length(freq)	
			else
			    min(which(freq>=x))	
			})
	
as.numeric(names(freq)[index])
}

Labels <- function(riskgroup){
		the_levels <- levels(riskgroup)
		the_levels <- gsub("(\\(|\\])","",the_levels)
		the_levels <- strsplit(the_levels, split=",")
	    first <- as.numeric(sapply(the_levels, function(x) x[1]))
	    second <- as.numeric(sapply(the_levels, function(x) x[2]))
	    first <- format(round(first*100,2)+rep(c(0,0.1,0),c(2,length(the_levels)-3,1)),ns=2,trim=TRUE)
	    second <- format(round(second*100,2), ns=2, trim=TRUE)
	    labels <- paste(first,"-",second,"%",sep="")
	    labels[1] <- paste("<",second[1],"%",sep="")
	    labels[length(labels)] <- paste(">",second[length(labels)-1],"%",sep="")
labels	    
}

by <- 0.2

#######################################################################################
#FOR IMPUTED DATSET 1
#######################################################################################
load(file="nhis_risk1.RData")

#1 - eligible group, 2 - not eligible group, but between 50<=age<=80, 3 - not eligible group outside age range;
nhis$eligible.group0 <- ifelse(nhis$medicare.eligible,1,ifelse(nhis$age>=50&nhis$age<=80,2,3))
nhis$eligible.group <- ifelse(nhis$uspstf.eligible,1,ifelse(nhis$age>=50&nhis$age<=80,2,3))

nhis$age.cat <- factor(nhis$age>65, lab=c("65-",">65"))
nhis$female.cat <- factor(nhis$female, lab=c("Male","Female"))
nhis$smoke.stat.cat <- factor(nhis$current, lab=c("Former","Current"))
nhis$packyears.cat <- factor(cut(nhis$packyears, breaks=c(0,20,40,60,max(nhis$packyears,na.rm=T)+1),right=F), exclude=NULL)
levels(nhis$age.cat)

#######################################################################################
#create quartile/quintile categories based on different populations
#######################################################################################
range <- range(nhis$risk)

#categorizes nhis$risk into quintiles based on Medicare-eligible group
breaks <- pop.quantile(nhis$risk[nhis$eligible.group0==1], nhis$wt[nhis$eligible.group0==1], seq(0,1,by=by))
breaks[1] <- range[1]-1
breaks[length(breaks)] <- range[2]+2
nhis$riskgroup1 <- cut(nhis$risk, breaks=breaks)


#######################################################################################
#calculate observed, expected, E/O, Var(E/O) for the following populations:
#overall; uspstf eligible; Ineligible, 50-80; All Ineligible; All 50-58; Other Ineligible
#######################################################################################
master  <- svydesign(id=~psu, strata=~strata, weights=~wt, data=nhis, nest=TRUE)
overall1 <- popEO(master)

design <- subset(master, eligible.group==1)
pop.age11 <- popriskEO("age.cat", design )
pop.female11 <- popriskEO("female.cat", design)
pop.race11 <- popriskEO("race", design )
pop.smoke.stat11 <- popriskEO("smoke.stat.cat", design )
pop.packyears11 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==2)
pop.age21 <- popriskEO("age.cat", design )
pop.female21 <- popriskEO("female.cat", design)
pop.race21 <- popriskEO("race", design )
pop.smoke.stat21 <- popriskEO("smoke.stat.cat", design )
pop.packyears21 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==3)
pop.age31 <- popriskEO("age.cat", design )
pop.female31 <- popriskEO("female.cat", design)
pop.race31 <- popriskEO("race", design )
pop.smoke.stat31 <- popriskEO("smoke.stat.cat", design )
pop.packyears31 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group!=3)
pop.age41 <- popriskEO("age.cat", design )
pop.female41 <- popriskEO("female.cat", design)
pop.race41 <- popriskEO("race", design )
pop.smoke.stat41 <- popriskEO("smoke.stat.cat", design )
pop.packyears41 <- popriskEO("packyears.cat", design )

#######################################################################################
#calculate observed, expected, E/O, Var(E/O) of quartiles/quintiles for the following populations:
#Other Ineligible overall; uspstf eligible; Ineligible, 50-80; All 50-58; 
#######################################################################################

design <- subset(master, eligible.group==3)
risk1 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group==1)
risk11 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group==2)
risk21 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group!=3)
risk31 <- popriskEO("riskgroup1", design) 


#######################################################################################
#FOR IMPUTED DATSET 2
#######################################################################################
load(file="nhis_risk2.RData")

nhis$eligible.group0 <- ifelse(nhis$medicare.eligible,1,ifelse(nhis$age>=50&nhis$age<=80,2,3))
nhis$eligible.group <- ifelse(nhis$uspstf.eligible,1,ifelse(nhis$age>=50&nhis$age<=80,2,3))

nhis$age.cat <- factor(nhis$age>65, lab=c("65-",">65"))
nhis$female.cat <- factor(nhis$female, lab=c("Male","Female"))
nhis$smoke.stat.cat <- factor(nhis$current, lab=c("Former","Current"))
nhis$packyears.cat <- factor(cut(nhis$packyears, breaks=c(0,20,40,60,max(nhis$packyears,na.rm=T)+1),right=F), exclude=NULL)

range <- range(nhis$risk)

breaks <- pop.quantile(nhis$risk[nhis$eligible.group0==1], nhis$wt[nhis$eligible.group0==1], seq(0,1,by=by))
breaks[1] <- range[1]-1
breaks[length(breaks)] <- range[2]+2
nhis$riskgroup1 <- cut(nhis$risk, breaks=breaks)


master  <- svydesign(id=~psu, strata=~strata, weights=~wt, data=nhis, nest=TRUE)
overall2 <- popEO(master)

design <- subset(master, eligible.group==1)
pop.age12 <- popriskEO("age.cat", design )
pop.female12 <- popriskEO("female.cat", design)
pop.race12 <- popriskEO("race", design )
pop.smoke.stat12 <- popriskEO("smoke.stat.cat", design )
pop.packyears12 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==2)
pop.age22 <- popriskEO("age.cat", design )
pop.female22 <- popriskEO("female.cat", design)
pop.race22 <- popriskEO("race", design )
pop.smoke.stat22 <- popriskEO("smoke.stat.cat", design )
pop.packyears22 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==3)
pop.age32 <- popriskEO("age.cat", design )
pop.female32 <- popriskEO("female.cat", design)
pop.race32 <- popriskEO("race", design )
pop.smoke.stat32 <- popriskEO("smoke.stat.cat", design )
pop.packyears32 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group!=3)
pop.age42 <- popriskEO("age.cat", design )
pop.female42 <- popriskEO("female.cat", design)
pop.race42 <- popriskEO("race", design )
pop.smoke.stat42 <- popriskEO("smoke.stat.cat", design )
pop.packyears42 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==3)
risk2 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group==1)
risk12 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group==2)
risk22 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group!=3)
risk32 <- popriskEO("riskgroup1", design) 


#######################################################################################
#FOR IMPUTED DATSET 3
#######################################################################################
load(file="nhis_risk3.RData")

nhis$eligible.group0 <- ifelse(nhis$medicare.eligible,1,ifelse(nhis$age>=50&nhis$age<=80,2,3))
nhis$eligible.group <- ifelse(nhis$uspstf.eligible,1,ifelse(nhis$age>=50&nhis$age<=80,2,3))

nhis$age.cat <- factor(nhis$age>65, lab=c("65-",">65"))
nhis$female.cat <- factor(nhis$female, lab=c("Male","Female"))
nhis$smoke.stat.cat <- factor(nhis$current, lab=c("Former","Current"))
nhis$packyears.cat <- factor(cut(nhis$packyears, breaks=c(0,20,40,60,max(nhis$packyears,na.rm=T)+1),right=F), exclude=NULL)

range <- range(nhis$risk)

breaks <- pop.quantile(nhis$risk[nhis$eligible.group0==1], nhis$wt[nhis$eligible.group0==1], seq(0,1,by=by))
breaks[1] <- range[1]-1
breaks[length(breaks)] <- range[2]+2
nhis$riskgroup1 <- cut(nhis$risk, breaks=breaks)

master  <- svydesign(id=~psu, strata=~strata, weights=~wt, data=nhis, nest=TRUE)
overall3 <- popEO(master)


design <- subset(master, eligible.group==1)
pop.age13 <- popriskEO("age.cat", design )
pop.female13 <- popriskEO("female.cat", design)
pop.race13 <- popriskEO("race", design )
pop.smoke.stat13 <- popriskEO("smoke.stat.cat", design )
pop.packyears13 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==2)
pop.age23 <- popriskEO("age.cat", design )
pop.female23 <- popriskEO("female.cat", design)
pop.race23 <- popriskEO("race", design )
pop.smoke.stat23 <- popriskEO("smoke.stat.cat", design )
pop.packyears23 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==3)
pop.age33 <- popriskEO("age.cat", design )
pop.female33 <- popriskEO("female.cat", design)
pop.race33 <- popriskEO("race", design )
pop.smoke.stat33 <- popriskEO("smoke.stat.cat", design )
pop.packyears33 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group!=3)
pop.age43 <- popriskEO("age.cat", design )
pop.female43 <- popriskEO("female.cat", design)
pop.race43 <- popriskEO("race", design )
pop.smoke.stat43 <- popriskEO("smoke.stat.cat", design )
pop.packyears43 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==3)
risk3 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group==1)
risk13 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group==2)
risk23 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group!=3)
risk33 <- popriskEO("riskgroup1", design) 


#######################################################################################
#FOR IMPUTED DATSET 4
#######################################################################################
load(file="nhis_risk4.RData")

nhis$eligible.group0 <- ifelse(nhis$medicare.eligible,1,ifelse(nhis$age>=50&nhis$age<=80,2,3))
nhis$eligible.group <- ifelse(nhis$uspstf.eligible,1,ifelse(nhis$age>=50&nhis$age<=80,2,3))

nhis$age.cat <- factor(nhis$age>65, lab=c("65-",">65"))
nhis$female.cat <- factor(nhis$female, lab=c("Male","Female"))
nhis$smoke.stat.cat <- factor(nhis$current, lab=c("Former","Current"))
nhis$packyears.cat <- factor(cut(nhis$packyears, breaks=c(0,20,40,60,max(nhis$packyears,na.rm=T)+1),right=F), exclude=NULL)

range <- range(nhis$risk)

breaks <- pop.quantile(nhis$risk[nhis$eligible.group0==1], nhis$wt[nhis$eligible.group0==1], seq(0,1,by=by))
breaks[1] <- range[1]-1
breaks[length(breaks)] <- range[2]+2
nhis$riskgroup1 <- cut(nhis$risk, breaks=breaks)

master  <- svydesign(id=~psu, strata=~strata, weights=~wt, data=nhis, nest=TRUE)
overall4 <- popEO(master)

design <- subset(master, eligible.group==1)
pop.age14 <- popriskEO("age.cat", design )
pop.female14 <- popriskEO("female.cat", design)
pop.race14 <- popriskEO("race", design )
pop.smoke.stat14 <- popriskEO("smoke.stat.cat", design )
pop.packyears14 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==2)
pop.age24 <- popriskEO("age.cat", design )
pop.female24 <- popriskEO("female.cat", design)
pop.race24 <- popriskEO("race", design )
pop.smoke.stat24 <- popriskEO("smoke.stat.cat", design )
pop.packyears24 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==3)
pop.age34 <- popriskEO("age.cat", design )
pop.female34 <- popriskEO("female.cat", design)
pop.race34 <- popriskEO("race", design )
pop.smoke.stat34 <- popriskEO("smoke.stat.cat", design )
pop.packyears34 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group!=3)
pop.age44 <- popriskEO("age.cat", design )
pop.female44 <- popriskEO("female.cat", design)
pop.race44 <- popriskEO("race", design )
pop.smoke.stat44 <- popriskEO("smoke.stat.cat", design )
pop.packyears44 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==3)
risk4 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group==1)
risk14 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group==2)
risk24 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group!=3)
risk34 <- popriskEO("riskgroup1", design) 


#######################################################################################
#FOR IMPUTED DATSET 5
#######################################################################################
load(file="nhis_risk5.RData")

nhis$eligible.group0 <- ifelse(nhis$medicare.eligible,1,ifelse(nhis$age>=50&nhis$age<=80,2,3))
nhis$eligible.group <- ifelse(nhis$uspstf.eligible,1,ifelse(nhis$age>=50&nhis$age<=80,2,3))

nhis$age.cat <- factor(nhis$age>65, lab=c("65-",">65"))
nhis$female.cat <- factor(nhis$female, lab=c("Male","Female"))
nhis$smoke.stat.cat <- factor(nhis$current, lab=c("Former","Current"))
nhis$packyears.cat <- factor(cut(nhis$packyears, breaks=c(0,20,40,60,max(nhis$packyears,na.rm=T)+1),right=F), exclude=NULL)

range <- range(nhis$risk)

breaks <- pop.quantile(nhis$risk[nhis$eligible.group0==1], nhis$wt[nhis$eligible.group0==1], seq(0,1,by=by))
breaks[1] <- range[1]-1
breaks[length(breaks)] <- range[2]+2
nhis$riskgroup1 <- cut(nhis$risk, breaks=breaks)

master  <- svydesign(id=~psu, strata=~strata, weights=~wt, data=nhis, nest=TRUE)
overall5 <- popEO(master)


design <- subset(master, eligible.group==1)
pop.age15 <- popriskEO("age.cat", design )
pop.female15 <- popriskEO("female.cat", design)
pop.race15 <- popriskEO("race", design )
pop.smoke.stat15 <- popriskEO("smoke.stat.cat", design )
pop.packyears15 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==2)
pop.age25 <- popriskEO("age.cat", design )
pop.female25 <- popriskEO("female.cat", design)
pop.race25 <- popriskEO("race", design )
pop.smoke.stat25 <- popriskEO("smoke.stat.cat", design )
pop.packyears25 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==3)
pop.age35 <- popriskEO("age.cat", design )
pop.female35 <- popriskEO("female.cat", design)
pop.race35 <- popriskEO("race", design )
pop.smoke.stat35 <- popriskEO("smoke.stat.cat", design )
pop.packyears35 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group!=3)
pop.age45 <- popriskEO("age.cat", design )
pop.female45 <- popriskEO("female.cat", design)
pop.race45 <- popriskEO("race", design )
pop.smoke.stat45 <- popriskEO("smoke.stat.cat", design )
pop.packyears45 <- popriskEO("packyears.cat", design )

design <- subset(master, eligible.group==3)
risk5 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group==1)
risk15 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group==2)
risk25 <- popriskEO("riskgroup1", design) 

design <- subset(master, eligible.group!=3)
risk35 <- popriskEO("riskgroup1", design) 

#############################################
# SUMMARY
#############################################
# For Overall
#means rows: 1:6 observed for the 6 groups, 7:12 expected for the 6 groups, 13:18 E/O for the 6 groups
#means columns: the 5 simulation datasets
means <- cbind( as.numeric(overall1[,1:3]),
						as.numeric(overall2[,1:3]),
						as.numeric(overall3[,1:3]),
						as.numeric(overall4[,1:3]),
						as.numeric(overall5[,1:3]))

#vars rows: var(E/O) for the 6 groups						
#vars columns: the 5 simulation datasets
vars  <- cbind(overall1[,4],overall2[,4],overall3[,4],overall4[,4],overall5[,4])

#average over the 5 simulations							
means.bar <- rowMeans(means)
vars.bar <- rowMeans(vars)

#v.impute is the sample variance in E/O
#V(E/O) = (E/O - E(E/O))^2/n
v.impute <- rowMeans((means[13:18,]-means.bar[13:18])^2/5)

#Calculate confidence intervals for mean E/O
#logv - delta method derivation of log(E(E/O))?
logv <- (vars.bar+(1+1/5)*v.impute)/means.bar[13:18]^2
lower <- exp(log(means.bar[13:18])-1.96*sqrt(logv))
upper <- exp(log(means.bar[13:18])+1.96*sqrt(logv))

#summary of observed, E/O, and confidence intervals
calibration <- rbind(
				Observed = means.bar[1:6],
			    EO = means.bar[13:18],
			    lower = lower,
			    upper = upper,
			    group = row.names(overall1)
			    )

# RISK GROUPS
#combines the observed for each quartile for each simulated other ineligible group
O <- cbind(
	sapply(risk1, function(x)x[2]),
	sapply(risk2, function(x)x[2]),
	sapply(risk3, function(x)x[2]),
	sapply(risk4, function(x)x[2]),
	sapply(risk5, function(x)x[2])
	)

#combines the expected/observed for each quartile for each simulated other ineligible group
EO <- cbind(
	sapply(risk1, function(x)x[3]),
	sapply(risk2, function(x)x[3]),
	sapply(risk3, function(x)x[3]),
	sapply(risk4, function(x)x[3]),
	sapply(risk5, function(x)x[3])
	)	

#combines the V(E/O) for each quartile for each simulated other ineligible group
V	<- cbind(
	sapply(risk1, function(x)x[4]),
	sapply(risk2, function(x)x[4]),
	sapply(risk3, function(x)x[4]),
	sapply(risk4, function(x)x[4]),
	sapply(risk5, function(x)x[4])
	)	

#means of E/O, V(E/O), and sample variance of E/O	
means.bar <- rowMeans(EO)
vars.bar <- rowMeans(V)
v.impute <- rowMeans((EO-means.bar)^2/5)

#confidence intervals
logv <- (vars.bar+(1+1/5)*v.impute)/means.bar^2
lower <- exp(log(means.bar)-1.96*sqrt(logv))
upper <- exp(log(means.bar)+1.96*sqrt(logv))

#mean observed, E/O, lower and upper CI of E/O for each quartile in the other ineligible group
overall.risk.calibration <- cbind(
	O = rowMeans(O), 
	EO = means.bar,
	lower = lower,
	upper = upper
)


#for each quintile of eligible group
# RISK GROUPS 
O <- cbind(
	sapply(risk11, function(x)x[2]),
	sapply(risk12, function(x)x[2]),
	sapply(risk13, function(x)x[2]),
	sapply(risk14, function(x)x[2]),
	sapply(risk15, function(x)x[2])
	)
	
EO <- cbind(
	sapply(risk11, function(x)x[3]),
	sapply(risk12, function(x)x[3]),
	sapply(risk13, function(x)x[3]),
	sapply(risk14, function(x)x[3]),
	sapply(risk15, function(x)x[3])
	)	

V	<- cbind(
	sapply(risk11, function(x)x[4]),
	sapply(risk12, function(x)x[4]),
	sapply(risk13, function(x)x[4]),
	sapply(risk14, function(x)x[4]),
	sapply(risk15, function(x)x[4])
	)	
	
means.bar <- rowMeans(EO)
vars.bar <- rowMeans(V)
v.impute <- rowMeans((EO-means.bar)^2/5)

logv <- (vars.bar+(1+1/5)*v.impute)/means.bar^2
lower <- exp(log(means.bar)-1.96*sqrt(logv))
upper <- exp(log(means.bar)+1.96*sqrt(logv))

overall.risk1.calibration <- cbind(
	O = rowMeans(O), 
	EO = means.bar,
	lower = lower,
	upper = upper
)


#for each quintile of ineligible 50-80 group
O <- cbind(
	sapply(risk21, function(x)x[2]),
	sapply(risk22, function(x)x[2]),
	sapply(risk23, function(x)x[2]),
	sapply(risk24, function(x)x[2]),
	sapply(risk25, function(x)x[2])
	)
	
EO <- cbind(
	sapply(risk21, function(x)x[3]),
	sapply(risk22, function(x)x[3]),
	sapply(risk23, function(x)x[3]),
	sapply(risk24, function(x)x[3]),
	sapply(risk25, function(x)x[3])
	)	

V	<- cbind(
	sapply(risk21, function(x)x[4]),
	sapply(risk22, function(x)x[4]),
	sapply(risk23, function(x)x[4]),
	sapply(risk24, function(x)x[4]),
	sapply(risk25, function(x)x[4])
	)	
	
means.bar <- rowMeans(EO)
vars.bar <- rowMeans(V)
v.impute <- rowMeans((EO-means.bar)^2/5)

logv <- (vars.bar+(1+1/5)*v.impute)/means.bar^2
lower <- exp(log(means.bar)-1.96*sqrt(logv))
upper <- exp(log(means.bar)+1.96*sqrt(logv))

overall.risk2.calibration <- cbind(
	O = rowMeans(O), 
	EO = means.bar,
	lower = lower,
	upper = upper
)


#for each quintile of 50-80 group
O <- cbind(
	sapply(risk31, function(x)x[2]),
	sapply(risk32, function(x)x[2]),
	sapply(risk33, function(x)x[2]),
	sapply(risk34, function(x)x[2]),
	sapply(risk35, function(x)x[2])
	)
	
EO <- cbind(
	sapply(risk31, function(x)x[3]),
	sapply(risk32, function(x)x[3]),
	sapply(risk33, function(x)x[3]),
	sapply(risk34, function(x)x[3]),
	sapply(risk35, function(x)x[3])
	)	

V	<- cbind(
	sapply(risk31, function(x)x[4]),
	sapply(risk32, function(x)x[4]),
	sapply(risk33, function(x)x[4]),
	sapply(risk34, function(x)x[4]),
	sapply(risk35, function(x)x[4])
	)	
	
means.bar <- rowMeans(EO)
vars.bar <- rowMeans(V)
v.impute <- rowMeans((EO-means.bar)^2/5)

logv <- (vars.bar+(1+1/5)*v.impute)/means.bar^2
lower <- exp(log(means.bar)-1.96*sqrt(logv))
upper <- exp(log(means.bar)+1.96*sqrt(logv))

overall.risk3.calibration <- cbind(
	O = rowMeans(O), 
	EO = means.bar,
	lower = lower,
	upper = upper
)


#rows 1-4: quartiles of other ineligible group, overall other ineligible group, and quintile of ineligible 50-80, eligible 50-80, and all 50-80 group
calibration.risk.groups <- rbind(overall.risk.calibration,
	overall.risk1.calibration,
	overall.risk2.calibration,
	overall.risk3.calibration
)


risklabels.eligible <- Labels(nhis$riskgroup1)  #eligible
risklabels.ineligible <- Labels(nhis$riskgroup1)  #ineligible 50-80
risklabels.all <- Labels(nhis$riskgroup1)  #All 50-80
risklabels.other.ineligible <- Labels(nhis$riskgroup1)  #Other ineligible

row.names(calibration.risk.groups) <- 
		c(paste(rep("Other Ineligible",each=5),paste("Q",1:5)),
		  paste(rep(c("Eligible","Ineligible, 50-80","All 50-80"),each=5), paste("Q",1:5)))


CalibrationSummary <- function(...){
	objects <- list(...)
	
	O <- cbind(
	sapply(objects[[1]], function(x)x[2]),
	sapply(objects[[2]], function(x)x[2]),
	sapply(objects[[3]], function(x)x[2]),
	sapply(objects[[4]], function(x)x[2]),
	sapply(objects[[5]], function(x)x[2])
	)
	
EO <- cbind(
	sapply(objects[[1]], function(x)x[3]),
	sapply(objects[[2]], function(x)x[3]),
	sapply(objects[[3]], function(x)x[3]),
	sapply(objects[[4]], function(x)x[3]),
	sapply(objects[[5]], function(x)x[3])
	)

V	<- cbind(
	sapply(objects[[1]], function(x)x[4]),
	sapply(objects[[2]], function(x)x[4]),
	sapply(objects[[3]], function(x)x[4]),
	sapply(objects[[4]], function(x)x[4]),
	sapply(objects[[5]], function(x)x[4])
	)
	
means.bar <- rowMeans(EO)
vars.bar <- rowMeans(V)
v.impute <- rowMeans((EO-means.bar)^2/5)

logv <- (vars.bar+(1+1/5)*v.impute)/means.bar^2
lower <- exp(log(means.bar)-1.96*sqrt(logv))
upper <- exp(log(means.bar)+1.96*sqrt(logv))

cbind(
	O = rowMeans(O), 
	EO = means.bar,
	lower = lower,
	upper = upper
)

}

age.summary <- CalibrationSummary(pop.age11,pop.age12,pop.age13, pop.age14,pop.age15)
female.summary <- CalibrationSummary(pop.female11,pop.female12,pop.female13, pop.female14,pop.female15)
race.summary <- CalibrationSummary(pop.race11,pop.race12,pop.race13, pop.race14,pop.race15)
smoke.stat.summary <- CalibrationSummary(pop.smoke.stat11,pop.smoke.stat12,pop.smoke.stat13, pop.smoke.stat14,pop.smoke.stat15)
packyears.summary <- CalibrationSummary(pop.packyears11,pop.packyears12,pop.packyears13, pop.packyears14,pop.packyears15)

summary.matrix.eligible <- rbind(
		age.summary,
		female.summary,
		race.summary,
		smoke.stat.summary,
		packyears.summary
	)

age.summary <- CalibrationSummary(pop.age21,pop.age22,pop.age23, pop.age24,pop.age25)
female.summary <- CalibrationSummary(pop.female21,pop.female22,pop.female23, pop.female24,pop.female25)
race.summary <- CalibrationSummary(pop.race21,pop.race22,pop.race23, pop.race24,pop.race25)
smoke.stat.summary <- CalibrationSummary(pop.smoke.stat21,pop.smoke.stat22,pop.smoke.stat23, pop.smoke.stat24,pop.smoke.stat25)
packyears.summary <- CalibrationSummary(pop.packyears21,pop.packyears22,pop.packyears23, pop.packyears24,pop.packyears25)

summary.matrix.ineligible5080 <- rbind(
		age.summary,
		female.summary,
		race.summary,
		smoke.stat.summary,
		packyears.summary
	)
	
age.summary <- CalibrationSummary(pop.age31,pop.age32,pop.age33, pop.age34,pop.age35)
female.summary <- CalibrationSummary(pop.female31,pop.female32,pop.female33, pop.female34,pop.female35)
race.summary <- CalibrationSummary(pop.race31,pop.race32,pop.race33, pop.race34,pop.race35)
smoke.stat.summary <- CalibrationSummary(pop.smoke.stat31,pop.smoke.stat32,pop.smoke.stat33, pop.smoke.stat34,pop.smoke.stat35)
packyears.summary <- CalibrationSummary(pop.packyears31,pop.packyears32,pop.packyears33, pop.packyears34,pop.packyears35)

summary.matrix.ineligible2 <- rbind(
		age.summary,
		female.summary,
		race.summary,
		smoke.stat.summary,
		packyears.summary
	)	
	
	
age.summary <- CalibrationSummary(pop.age41,pop.age42,pop.age43, pop.age44,pop.age45)
female.summary <- CalibrationSummary(pop.female41,pop.female42,pop.female43, pop.female44,pop.female45)
race.summary <- CalibrationSummary(pop.race41,pop.race42,pop.race43, pop.race44,pop.race45)
smoke.stat.summary <- CalibrationSummary(pop.smoke.stat41,pop.smoke.stat42,pop.smoke.stat43, pop.smoke.stat44,pop.smoke.stat45)
packyears.summary <- CalibrationSummary(pop.packyears41,pop.packyears42,pop.packyears43, pop.packyears44,pop.packyears45)

summary.matrix.all5080 <- rbind(
		age.summary,
		female.summary,
		race.summary,
		smoke.stat.summary,
		packyears.summary
	)

rownames(summary.matrix.eligible) <- c("Age <=65","Age >65",
                                      "Male","Female",
                                      "White, not hispanic","Black, not hispanic","Hispanic","Other",
                                      "Former","Current",
                                      "30-<40 pack-years","40-<60 pack-years","60+ pack-years")

rownames(summary.matrix.ineligible5080) <- c("Age <=65","Age >65",
                                      "Male","Female",
                                      "White, not hispanic","Black, not hispanic","Hispanic","Other",
                                      "Former","Current",
                                      "<20 pack-years","20-<40 pack-years","40-<60 pack-years","60+ pack-years")

rownames(summary.matrix.ineligible2) <- c("Age <=65","Age >65",
                                      "Male","Female",
                                      "White, not hispanic","Black, not hispanic","Hispanic","Other",
                                      "Former","Current",
                                      "<20 pack-years","20-<40 pack-years","40-<60 pack-years","60+ pack-years")	
	
rownames(summary.matrix.all5080) <- c("Age <=65","Age >65",
                                      "Male","Female",
                                      "White, not hispanic","Black, not hispanic","Hispanic","Other",
                                      "Former","Current",
                                      "<20 pack-years","20-<40 pack-years","40-<60 pack-years","60+ pack-years")	