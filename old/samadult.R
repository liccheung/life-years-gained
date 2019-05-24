# CREATE ANALYTIC DATA SETS FOR NHIS
nhis <- read.csv(file="~/master/project/nhis/data/samadult/nhis_1997.csv")

#NAME CHANGES
names(nhis)[names(nhis)=="PUBLICID"] <- "pid"
names(nhis)[names(nhis)=="AGE_P"] <- "age"
names(nhis)[names(nhis)=="BMI"] <- "bmi"

names(nhis)[names(nhis)=="WTFA_SA"]  <- "wt"
names(nhis)[names(nhis)=="PSU"]  <- "psu"
names(nhis)[names(nhis)=="STRATUM"]  <- "strata"


#DERIVED VARIABLES
nhis$year <- 1997
nhis$intyear <- nhis$year+nhis$INTV_QRT*1.5/12
nhis$female <- ifelse(nhis$SEX=="Female",1,0)
nhis$race <- ifelse(nhis$ORIGIN=="Yes","Hispanic",
				ifelse(nhis$RACE=="White","Non-Hispanic White",
					ifelse(nhis$RACE=="Black","Non-Hispanic Black",
						ifelse(nhis$RACE=="AIAN*","American Indian/Alaskan Native",
							ifelse(nhis$RACE=="API*","Asian",NA)))))

nhis$bmi[nhis$bmi==99.99] <- NA
							  
nhis$mining <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==10,1,0)
nhis$textile <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==31,1,0)
nhis$chemical <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==33,1,0)
nhis$metal <- ifelse(!is.na(nhis$INDSTRY1)&(nhis$INDSTRY1==41|nhis$INDSTRY1)==42,1,0)
nhis$handler <- ifelse(!is.na(nhis$OCCUP1)&(nhis$OCCUP1==33|nhis$OCCUP1==36|nhis$OCCUP1==40),1,0)

nhis$asb <- ifelse(nhis$mining==1|(nhis$chemical==1&nhis$handler==1),1,0)
nhis$dust <- ifelse((nhis$textile==1|nhis$metal==1)&nhis$handler==1,1,0)
nhis$emp <- ifelse(nhis$EPHEV=="Yes",1,0)
nhis$copd <- nhis$emp

nhis$current <- ifelse(nhis$SMKSTAT1=="Current"|nhis$SMKSTAT1=="Smoker, current status unknown",1,0)										  
nhis$former <- ifelse(nhis$SMKSTAT1=="Former",1,0)											  
nhis$never <- ifelse(nhis$SMKSTAT1=="Never",1,0)											  
nhis$unknown <-  ifelse(nhis$SMKSTAT1=="Unknown if ever smoked",1,0)	

nhis$edu <- ifelse(nhis$EDUC=="12th grade, no diploma"|nhis$EDUC=="Grades 1 - 11","8-11 years",
								ifelse(nhis$EDUC=="HIGH SCHOOL GRADUATE","12 years or completed high school",
									ifelse(nhis$EDUC=="AA degree: technical or vocational"|nhis$EDUC=="GED or equivalent",
									"Post-high school training other than college",
									ifelse(nhis$EDUC=="Never attended/ kindergarten only","Less than 8 years",
										ifelse(nhis$EDUC=="Some college, no degree","Some college",
											ifelse(nhis$EDUC== "AA degree: academic program"|nhis$EDUC=="Bachelor's degree (BA, AB, BS, BBA)",
												"College graduate",
												ifelse(nhis$EDUC=="Not Ascertained"|nhis$EDUC=="Refused"|nhis$EDUC=="Don't know","Missing","Post graduate")))))))
nhis$fam.cancer <- 0
nhis$fam.lung.cancer <- 0
nhis$fam.lung.trend <- 0
nhis$fam.smoke.cancer <- 0
nhis$fam.cancer.onset <- 0
nhis$no.asthma <- ifelse(nhis$AHAYFYR=="Yes",0,1) # 12 month
nhis$prior.cancer <- ifelse(nhis$CANEV=="Yes",1,0)
nhis$lung.cancer.before <- ifelse(nhis$CNKIND14=="Mentioned",1,0)

nhis$cpd <- nhis$CIGSDAY 
nhis$cpd[nhis$cpd==97|nhis$cpd==98|nhis$cpd==99] <- NA
	
nhis$smoke.age.start <- nhis$SMKREG
nhis$smoke.age.start[nhis$smoke.age.start==96|nhis$smoke.age.start==97|nhis$smoke.age.start==98|nhis$smoke.age.start==99] <- NA
	
nhis$qtyears <- nhis$SMKQTY 	
nhis$qtyears[nhis$qtyears==97|nhis$qtyears==98|nhis$qtyears==99] <- NA
nhis$qtyears[!is.na(nhis$qtyears)&nhis$qtyears==0] <- 1/2
# TREAT FORMER QUITTERS WITH UNKNOWN AS CURRENT
nhis$former[nhis$SMKQTD=="Yes"|nhis$SMKQTD=="No"|nhis$SMKQTD=="Don't know"] <- 0
nhis$current[nhis$SMKQTD=="Yes"|nhis$SMKQTD=="No"|nhis$SMKQTD=="Don't know"] <- 1
nhis$qtyears[nhis$current==1] <- 0

#IF FORMER SMOKER WITH UNKNOWN SMOKE AGE START BUT KNOWN QUIT YEARS
nhis$smoke.age.start[is.na(nhis$smoke.age.start)&nhis$former==1&!is.na(nhis$qtyears)] <- (nhis$age-nhis$qtyears)[is.na(nhis$smoke.age.start)&nhis$former==1&!is.na(nhis$qtyears)] 
nhis$smkyears <- nhis$age-nhis$smoke.age.start-nhis$qtyears 
nhis$smkyears <- ifelse(!is.na(nhis$smkyears)&nhis$smkyears<0,1/12,nhis$smkyears) 
nhis$birthyear <- 1997-nhis$age

# APPLY EXCLUSION CRITERIA; EVER SMOKERS, NO HISTORY OF LUNG CANCER
sum(nhis$lung.cancer.before)
sum(nhis$never)
sum(nhis$unknown)

nhis <- subset(nhis, unknown==0&never==0&lung.cancer.before==0)

keep <- c("pid",
	"psu",
	"strata",
    "intyear",
    "birthyear",
    "year",
	"age",
	"asb",
	"bmi",
	"copd",
	"cpd",
	"current",
	"dust",
	"edu",
	"emp",
	"fam.cancer",
	"fam.lung.cancer",
	"fam.lung.trend",
	"fam.smoke.cancer",
	"fam.cancer.onset",
	"female",
	"former",
	"no.asthma",
	"prior.cancer",
	"qtyears",
	"smkyears",
	"smoke.age.start",
	"race"
)

nhis1997 <- nhis[,keep]
save(nhis1997, file="~/master/project/nhis/data/samadult/nhis1997.RData")
rm(list=ls())



nhis <- read.csv(file="~/master/project/nhis/data/samadult/nhis_1998.csv")
#NAME CHANGES
names(nhis)[names(nhis)=="PUBLICID"] <- "pid"
names(nhis)[names(nhis)=="AGE_P"] <- "age"
names(nhis)[names(nhis)=="BMI"] <- "bmi"

names(nhis)[names(nhis)=="WTFA_SA"]  <- "wt"
names(nhis)[names(nhis)=="PSU"]  <- "psu"
names(nhis)[names(nhis)=="STRATUM"]  <- "strata"


#DERIVED VARIABLES
nhis$year <- 1998
nhis$intyear <- nhis$year+nhis$INTV_QRT*1.5/12
nhis$female <- ifelse(nhis$SEX=="Female",1,0)
nhis$race <- ifelse(nhis$ORIGIN=="Yes","Hispanic",
				ifelse(nhis$RACE=="White","Non-Hispanic White",
					ifelse(nhis$RACE=="Black","Non-Hispanic Black",
						ifelse(nhis$RACE=="AIAN*","American Indian/Alaskan Native",
							ifelse(nhis$RACE=="API*","Asian",NA)))))

nhis$bmi[nhis$bmi==99.99] <- NA
							  
nhis$mining <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==10,1,0)
nhis$textile <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==31,1,0)
nhis$chemical <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==33,1,0)
nhis$metal <- ifelse(!is.na(nhis$INDSTRY1)&(nhis$INDSTRY1==41|nhis$INDSTRY1)==42,1,0)
nhis$handler <- ifelse(!is.na(nhis$OCCUP1)&(nhis$OCCUP1==33|nhis$OCCUP1==36|nhis$OCCUP1==40),1,0)

nhis$asb <- ifelse(nhis$mining==1|(nhis$chemical==1&nhis$handler==1),1,0)
nhis$dust <- ifelse((nhis$textile==1|nhis$metal==1)&nhis$handler==1,1,0)
nhis$emp <- ifelse(nhis$EPHEV=="Yes",1,0)
nhis$copd <- nhis$emp

nhis$current <- ifelse(nhis$SMKSTAT1=="Current"|nhis$SMKSTAT1=="Smoker, current status unknown",1,0)										  
nhis$former <- ifelse(nhis$SMKSTAT1=="Former",1,0)											  
nhis$never <- ifelse(nhis$SMKSTAT1=="Never",1,0)											  
nhis$unknown <-  ifelse(nhis$SMKSTAT1=="Unknown if ever smoked",1,0)	

nhis$edu <- ifelse(nhis$EDUC=="12th grade, no diploma"|nhis$EDUC=="Grades 1 - 11","8-11 years",
								ifelse(nhis$EDUC=="HIGH SCHOOL GRADUATE","12 years or completed high school",
									ifelse(nhis$EDUC=="AA degree: technical or vocational"|nhis$EDUC=="GED or equivalent",
									"Post-high school training other than college",
									ifelse(nhis$EDUC=="Never attended/ kindergarten only","Less than 8 years",
										ifelse(nhis$EDUC=="Some college, no degree","Some college",
											ifelse(nhis$EDUC== "AA degree: academic program"|nhis$EDUC=="Bachelor's degree (BA, AB, BS, BBA)",
												"College graduate",
												ifelse(nhis$EDUC=="Not Ascertained"|nhis$EDUC=="Refused"|nhis$EDUC=="Don't know","Missing","Post graduate")))))))
nhis$fam.cancer <- 0
nhis$fam.lung.cancer <- 0
nhis$fam.lung.trend <- 0
nhis$fam.smoke.cancer <- 0
nhis$fam.cancer.onset <- 0
nhis$no.asthma <- ifelse(nhis$AHAYFYR=="Yes",0,1) # 12 month
nhis$prior.cancer <- ifelse(nhis$CANEV=="Yes",1,0)
nhis$lung.cancer.before <- ifelse(nhis$CNKIND14=="Mentioned",1,0)


nhis$cpd <- nhis$CIGSDAY 
nhis$cpd[nhis$cpd==97|nhis$cpd==98|nhis$cpd==99] <- NA
	
nhis$smoke.age.start <- nhis$SMKREG
nhis$smoke.age.start[nhis$smoke.age.start==96|nhis$smoke.age.start==97|nhis$smoke.age.start==98|nhis$smoke.age.start==99] <- NA
	
nhis$qtyears <- nhis$SMKQTY 	
nhis$qtyears[nhis$qtyears==97|nhis$qtyears==98|nhis$qtyears==99] <- NA
nhis$qtyears[!is.na(nhis$qtyears)&nhis$qtyears==0] <- 1/2
# TREAT FORMER QUITTERS WITH UNKNOWN AS CURRENT
nhis$former[nhis$SMKQTD=="Yes"|nhis$SMKQTD=="No"|nhis$SMKQTD=="Don't know"] <- 0
nhis$current[nhis$SMKQTD=="Yes"|nhis$SMKQTD=="No"|nhis$SMKQTD=="Don't know"] <- 1
nhis$qtyears[nhis$current==1] <- 0

#IF FORMER SMOKER WITH UNKNOWN SMOKE AGE START BUT KNOWN QUIT YEARS
nhis$smoke.age.start[is.na(nhis$smoke.age.start)&nhis$former==1&!is.na(nhis$qtyears)] <- (nhis$age-nhis$qtyears)[is.na(nhis$smoke.age.start)&nhis$former==1&!is.na(nhis$qtyears)] 
nhis$smkyears <- nhis$age-nhis$smoke.age.start-nhis$qtyears 
nhis$smkyears <- ifelse(!is.na(nhis$smkyears)&nhis$smkyears<0,1/12,nhis$smkyears) 
nhis$birthyear <- 1998-nhis$age

# APPLY EXCLUSION CRITERIA; EVER SMOKERS, NO HISTORY OF LUNG CANCER
sum(nhis$lung.cancer.before)
sum(nhis$never)
sum(nhis$unknown)

nhis <- subset(nhis, unknown==0&never==0&lung.cancer.before==0)

keep <- c("pid",
	"psu",
	"strata",
	"intyear",
    "birthyear",
    "year",
	"age",
	"asb",
	"bmi",
	"copd",
	"cpd",
	"current",
	"dust",
	"edu",
	"emp",
	"fam.cancer",
	"fam.lung.cancer",
	"fam.lung.trend",
	"fam.smoke.cancer",
	"fam.cancer.onset",
	"female",
	"former",
	"no.asthma",
	"prior.cancer",
	"qtyears",
	"smkyears",
	"smoke.age.start",
	"race"
)

nhis1998 <- nhis[,keep]
save(nhis1998, file="~/master/project/nhis/data/samadult/nhis1998.RData")
rm(list=ls())



nhis <- read.csv(file="~/master/project/nhis/data/samadult/nhis_1999.csv")
#NAME CHANGES
names(nhis)[names(nhis)=="PUBLICID"] <- "pid"
names(nhis)[names(nhis)=="AGE_P"] <- "age"
names(nhis)[names(nhis)=="BMI"] <- "bmi"
names(nhis)[names(nhis)=="RC_SUM_P"] <- "RACE"

names(nhis)[names(nhis)=="WTFA_SA"]  <- "wt"
names(nhis)[names(nhis)=="PSU"]  <- "psu"
names(nhis)[names(nhis)=="STRATUM"]  <- "strata"


#DERIVED VARIABLES
nhis$year <- 1999
nhis$intyear <- nhis$year+nhis$INTV_QRT*1.5/12
nhis$female <- ifelse(nhis$SEX=="Female",1,0)
nhis$race <- ifelse(nhis$ORIGIN=="Yes","Hispanic",
				ifelse(nhis$RACE=="White only","Non-Hispanic White",
					ifelse(nhis$RACE=="Black/African American only","Non-Hispanic Black",
						ifelse(nhis$RACE=="AIAN only*","American Indian/Alaskan Native",
							ifelse(nhis$RACE=="Asian only","Asian",NA)))))

nhis$bmi[nhis$bmi==99.99] <- NA
							  
nhis$mining <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==10,1,0)
nhis$textile <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==31,1,0)
nhis$chemical <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==33,1,0)
nhis$metal <- ifelse(!is.na(nhis$INDSTRY1)&(nhis$INDSTRY1==41|nhis$INDSTRY1)==42,1,0)
nhis$handler <- ifelse(!is.na(nhis$OCCUP1)&(nhis$OCCUP1==33|nhis$OCCUP1==36|nhis$OCCUP1==40),1,0)

nhis$asb <- ifelse(nhis$mining==1|(nhis$chemical==1&nhis$handler==1),1,0)
nhis$dust <- ifelse((nhis$textile==1|nhis$metal==1)&nhis$handler==1,1,0)
nhis$emp <- ifelse(nhis$EPHEV=="Yes",1,0)
nhis$copd <- nhis$emp

nhis$current <- ifelse(nhis$SMKSTAT1=="Current"|nhis$SMKSTAT1=="Smoker, current status unknown",1,0)										  
nhis$former <- ifelse(nhis$SMKSTAT1=="Former",1,0)											  
nhis$never <- ifelse(nhis$SMKSTAT1=="Never",1,0)											  
nhis$unknown <-  ifelse(nhis$SMKSTAT1=="Unknown if ever smoked",1,0)	

nhis$edu <- ifelse(nhis$EDUC=="12th grade, no diploma"|nhis$EDUC=="Grades 1 - 11","8-11 years",
								ifelse(nhis$EDUC=="HIGH SCHOOL GRADUATE","12 years or completed high school",
									ifelse(nhis$EDUC=="AA degree: technical or vocational"|nhis$EDUC=="GED or equivalent",
									"Post-high school training other than college",
									ifelse(nhis$EDUC=="Never attended/ kindergarten only","Less than 8 years",
										ifelse(nhis$EDUC=="Some college, no degree","Some college",
											ifelse(nhis$EDUC== "AA degree: academic program"|nhis$EDUC=="Bachelor's degree (BA, AB, BS, BBA)",
												"College graduate",
												ifelse(nhis$EDUC=="Not Ascertained"|nhis$EDUC=="Refused"|nhis$EDUC=="Don't know","Missing","Post graduate")))))))
nhis$fam.cancer <- 0
nhis$fam.lung.cancer <- 0
nhis$fam.lung.trend <- 0
nhis$fam.smoke.cancer <- 0
nhis$fam.cancer.onset <- 0
nhis$no.asthma <- ifelse(nhis$AHAYFYR=="Yes",0,1) # 12 month
nhis$prior.cancer <- ifelse(nhis$CANEV=="Yes",1,0)
nhis$lung.cancer.before <- ifelse(nhis$CNKIND14=="Mentioned",1,0)


nhis$cpd <- nhis$CIGSDAY 
nhis$cpd[nhis$cpd==97|nhis$cpd==98|nhis$cpd==99] <- NA
	
nhis$smoke.age.start <- nhis$SMKREG
nhis$smoke.age.start[nhis$smoke.age.start==96|nhis$smoke.age.start==97|nhis$smoke.age.start==98|nhis$smoke.age.start==99] <- NA
	
nhis$qtyears <- nhis$SMKQTY 	
nhis$qtyears[nhis$qtyears==97|nhis$qtyears==98|nhis$qtyears==99] <- NA
nhis$qtyears[!is.na(nhis$qtyears)&nhis$qtyears==0] <- 1/2
# TREAT FORMER QUITTERS WITH UNKNOWN AS CURRENT
nhis$former[nhis$SMKQTD=="Yes"|nhis$SMKQTD=="No"|nhis$SMKQTD=="Don't know"] <- 0
nhis$current[nhis$SMKQTD=="Yes"|nhis$SMKQTD=="No"|nhis$SMKQTD=="Don't know"] <- 1
nhis$qtyears[nhis$current==1] <- 0

#IF FORMER SMOKER WITH UNKNOWN SMOKE AGE START BUT KNOWN QUIT YEARS
nhis$smoke.age.start[is.na(nhis$smoke.age.start)&nhis$former==1&!is.na(nhis$qtyears)] <- (nhis$age-nhis$qtyears)[is.na(nhis$smoke.age.start)&nhis$former==1&!is.na(nhis$qtyears)] 
nhis$smkyears <- nhis$age-nhis$smoke.age.start-nhis$qtyears 
nhis$smkyears <- ifelse(!is.na(nhis$smkyears)&nhis$smkyears<0,1/12,nhis$smkyears) 
nhis$birthyear <- 1999-nhis$age

# APPLY EXCLUSION CRITERIA; EVER SMOKERS, NO HISTORY OF LUNG CANCER
sum(nhis$lung.cancer.before)
sum(nhis$never)
sum(nhis$unknown)

nhis <- subset(nhis, unknown==0&never==0&lung.cancer.before==0)

keep <- c("pid",
	"psu",
	"strata",
	"intyear",
    "birthyear",
    "year",
	"age",
	"asb",
	"bmi",
	"copd",
	"cpd",
	"current",
	"dust",
	"edu",
	"emp",
	"fam.cancer",
	"fam.lung.cancer",
	"fam.lung.trend",
	"fam.smoke.cancer",
	"fam.cancer.onset",
	"female",
	"former",
	"no.asthma",
	"prior.cancer",
	"qtyears",
	"smkyears",
	"smoke.age.start",
	"race"
)

nhis1999 <- nhis[,keep]
save(nhis1999, file="~/master/project/nhis/data/samadult/nhis1999.RData")
rm(list=ls())


nhis <- read.csv(file="~/master/project/nhis/data/samadult/nhis_2000.csv")
#NAME CHANGES
names(nhis)[names(nhis)=="RACERP_I"] <- "RACE"
names(nhis)[names(nhis)=="ORIGIN_I"] <- "ORIGIN"

names(nhis)[names(nhis)=="BMI"] <- "bmi"

names(nhis)[names(nhis)=="PUBLICID"] <- "pid"
names(nhis)[names(nhis)=="AGE_P"] <- "age"
names(nhis)[names(nhis)=="BMI"] <- "bmi"

names(nhis)[names(nhis)=="WTFA_SA"]  <- "wt"
names(nhis)[names(nhis)=="PSU"]  <- "psu"
names(nhis)[names(nhis)=="STRATUM"]  <- "strata"


#DERIVED VARIABLES
nhis$year <- 2000
nhis$intyear <- nhis$year+nhis$INTV_QRT*1.5/12
nhis$female <- ifelse(nhis$SEX=="Female",1,0)
nhis$race <- ifelse(nhis$ORIGIN=="Yes","Hispanic",
				ifelse(nhis$RACE=="White only","Non-Hispanic White",
					ifelse(nhis$RACE=="Black/African American only","Non-Hispanic Black",
						ifelse(nhis$RACE=="AIAN* only","American Indian/Alaskan Native",
							ifelse(nhis$RACE=="Asian only","Asian",NA)))))

nhis$bmi[nhis$bmi==99.99] <- NA
							  
nhis$mining <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==10,1,0)
nhis$textile <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==31,1,0)
nhis$chemical <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==33,1,0)
nhis$metal <- ifelse(!is.na(nhis$INDSTRY1)&(nhis$INDSTRY1==41|nhis$INDSTRY1)==42,1,0)
nhis$handler <- ifelse(!is.na(nhis$OCCUP1)&(nhis$OCCUP1==33|nhis$OCCUP1==36|nhis$OCCUP1==40),1,0)

nhis$asb <- ifelse(nhis$mining==1|(nhis$chemical==1&nhis$handler==1),1,0)
nhis$dust <- ifelse((nhis$textile==1|nhis$metal==1)&nhis$handler==1,1,0)
nhis$emp <- ifelse(nhis$EPHEV=="Yes",1,0)
nhis$copd <- nhis$emp

nhis$current <- ifelse(nhis$SMKSTAT1=="Current"|nhis$SMKSTAT1=="Smoker, current status unknown",1,0)										  
nhis$former <- ifelse(nhis$SMKSTAT1=="Former",1,0)											  
nhis$never <- ifelse(nhis$SMKSTAT1=="Never",1,0)											  
nhis$unknown <-  ifelse(nhis$SMKSTAT1=="Unknown if ever smoked",1,0)	

nhis$edu <- ifelse(nhis$EDUC=="12th grade, no diploma"|nhis$EDUC=="Grades 1 - 11","8-11 years",
								ifelse(nhis$EDUC=="HIGH SCHOOL GRADUATE","12 years or completed high school",
									ifelse(nhis$EDUC=="AA degree: technical or vocational"|nhis$EDUC=="GED or equivalent",
									"Post-high school training other than college",
									ifelse(nhis$EDUC=="Never attended/ kindergarten only","Less than 8 years",
										ifelse(nhis$EDUC=="Some college, no degree","Some college",
											ifelse(nhis$EDUC== "AA degree: academic program"|nhis$EDUC=="Bachelor's degree (BA, AB, BS, BBA)",
												"College graduate",
												ifelse(nhis$EDUC=="Not Ascertained"|nhis$EDUC=="Refused"|nhis$EDUC=="Don't know","Missing","Post graduate")))))))

nhis$fam.cancer <- ifelse(nhis$FHMCAN=="Yes"|nhis$FHFCAN=="Yes",1,0)
nhis$fam.lung.trend <- ifelse(nhis$FHFTYP14=="Mentioned",1,0)+ifelse(nhis$FHMTYP14=="Mentioned",1,0)
nhis$fam.lung.cancer <- ifelse(nhis$fam.lung.trend>=1,1,0)
nhis$fam.smoke.cancer <- nhis$fam.lung.cancer
nhis$fam.cancer.onset <- nhis$fam.lung.cancer
nhis$no.asthma <- ifelse(nhis$AHAYFYR=="Yes",0,1) # 12 month
nhis$prior.cancer <- ifelse(nhis$CANEV=="Yes",1,0)
nhis$lung.cancer.before <- ifelse(nhis$CNKIND14=="Mentioned",1,0)


nhis$cpd <- nhis$CIGSDAY 
nhis$cpd[nhis$cpd==97|nhis$cpd==98|nhis$cpd==99] <- NA
	
nhis$smoke.age.start <- nhis$SMKREG
nhis$smoke.age.start[nhis$smoke.age.start==96|nhis$smoke.age.start==97|nhis$smoke.age.start==98|nhis$smoke.age.start==99] <- NA
	
nhis$qtyears <- nhis$SMKQTY 	
nhis$qtyears[nhis$qtyears==97|nhis$qtyears==98|nhis$qtyears==99] <- NA
nhis$qtyears[!is.na(nhis$qtyears)&nhis$qtyears==0] <- 1/2
# TREAT FORMER QUITTERS WITH UNKNOWN AS CURRENT
nhis$former[nhis$SMKQTD=="Yes"|nhis$SMKQTD=="No"|nhis$SMKQTD=="Don't know"] <- 0
nhis$current[nhis$SMKQTD=="Yes"|nhis$SMKQTD=="No"|nhis$SMKQTD=="Don't know"] <- 1
nhis$qtyears[nhis$current==1] <- 0

#IF FORMER SMOKER WITH UNKNOWN SMOKE AGE START BUT KNOWN QUIT YEARS
nhis$smoke.age.start[is.na(nhis$smoke.age.start)&nhis$former==1&!is.na(nhis$qtyears)] <- (nhis$age-nhis$qtyears)[is.na(nhis$smoke.age.start)&nhis$former==1&!is.na(nhis$qtyears)] 
nhis$smkyears <- nhis$age-nhis$smoke.age.start-nhis$qtyears 
nhis$smkyears <- ifelse(!is.na(nhis$smkyears)&nhis$smkyears<0,1/12,nhis$smkyears) 
nhis$birthyear <- 2000-nhis$age

# APPLY EXCLUSION CRITERIA; EVER SMOKERS, NO HISTORY OF LUNG CANCER
sum(nhis$lung.cancer.before)
sum(nhis$never)
sum(nhis$unknown)

nhis <- subset(nhis, unknown==0&never==0&lung.cancer.before==0)

keep <- c("pid",
	"psu",
	"strata",
    "intyear",
    "birthyear",
    "year",
	"age",
	"asb",
	"bmi",
	"copd",
	"cpd",
	"current",
	"dust",
	"edu",
	"emp",
	"fam.cancer",
	"fam.lung.cancer",
	"fam.lung.trend",
	"fam.smoke.cancer",
	"fam.cancer.onset",
	"female",
	"former",
	"no.asthma",
	"prior.cancer",
	"qtyears",
	"smkyears",
	"smoke.age.start",
	"race"
)

nhis2000 <- nhis[,keep]
save(nhis2000, file="~/master/project/nhis/data/samadult/nhis2000.RData")
rm(list=ls())

nhis <- read.csv(file="~/master/project/nhis/data/samadult/nhis_2001.csv")
#NAME CHANGES
names(nhis)[names(nhis)=="RACERP_I"] <- "RACE"
names(nhis)[names(nhis)=="ORIGIN_I"] <- "ORIGIN"

names(nhis)[names(nhis)=="BMI"] <- "bmi"

names(nhis)[names(nhis)=="PUBLICID"] <- "pid"
names(nhis)[names(nhis)=="AGE_P"] <- "age"
names(nhis)[names(nhis)=="BMI"] <- "bmi"

names(nhis)[names(nhis)=="WTFA_SA"]  <- "wt"
names(nhis)[names(nhis)=="PSU"]  <- "psu"
names(nhis)[names(nhis)=="STRATUM"]  <- "strata"


#DERIVED VARIABLES
nhis$year <- 2001
nhis$intyear <- nhis$year+nhis$INTV_QRT*1.5/12
nhis$female <- ifelse(nhis$SEX=="Female",1,0)
nhis$race <- ifelse(nhis$ORIGIN=="Yes","Hispanic",
				ifelse(nhis$RACE=="White only","Non-Hispanic White",
					ifelse(nhis$RACE=="Black/African American only","Non-Hispanic Black",
						ifelse(nhis$RACE=="AIAN* only","American Indian/Alaskan Native",
							ifelse(nhis$RACE=="Asian only","Asian",NA)))))

nhis$bmi[nhis$bmi==99.99] <- NA
							  
nhis$mining <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==10,1,0)
nhis$textile <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==31,1,0)
nhis$chemical <- ifelse(!is.na(nhis$INDSTRY1)&nhis$INDSTRY1==33,1,0)
nhis$metal <- ifelse(!is.na(nhis$INDSTRY1)&(nhis$INDSTRY1==41|nhis$INDSTRY1)==42,1,0)
nhis$handler <- ifelse(!is.na(nhis$OCCUP1)&(nhis$OCCUP1==33|nhis$OCCUP1==36|nhis$OCCUP1==40),1,0)

nhis$asb <- ifelse(nhis$mining==1|(nhis$chemical==1&nhis$handler==1),1,0)
nhis$dust <- ifelse((nhis$textile==1|nhis$metal==1)&nhis$handler==1,1,0)
nhis$emp <- ifelse(nhis$EPHEV=="Yes",1,0)
nhis$copd <- nhis$emp

nhis$current <- ifelse(nhis$SMKSTAT1=="Current"|nhis$SMKSTAT1=="Smoker, current status unknown",1,0)										  
nhis$former <- ifelse(nhis$SMKSTAT1=="Former",1,0)											  
nhis$never <- ifelse(nhis$SMKSTAT1=="Never",1,0)											  
nhis$unknown <-  ifelse(nhis$SMKSTAT1=="Unknown if ever smoked",1,0)	

nhis$edu <- ifelse(nhis$EDUC=="12th grade, no diploma"|nhis$EDUC=="Grades 1 - 11","8-11 years",
								ifelse(nhis$EDUC=="HIGH SCHOOL GRADUATE","12 years or completed high school",
									ifelse(nhis$EDUC=="AA degree: technical or vocational"|nhis$EDUC=="GED or equivalent",
									"Post-high school training other than college",
									ifelse(nhis$EDUC=="Never attended/ kindergarten only","Less than 8 years",
										ifelse(nhis$EDUC=="Some college, no degree","Some college",
											ifelse(nhis$EDUC== "AA degree: academic program"|nhis$EDUC=="Bachelor's degree (BA, AB, BS, BBA)",
												"College graduate",
												ifelse(nhis$EDUC=="Not Ascertained"|nhis$EDUC=="Refused"|nhis$EDUC=="Don't know","Missing","Post graduate")))))))

nhis$fam.cancer <- 0
nhis$fam.lung.cancer <- 0
nhis$fam.lung.trend <- 0
nhis$fam.smoke.cancer <- 0
nhis$fam.cancer.onset <- 0
nhis$no.asthma <- ifelse(nhis$AHAYFYR=="Yes",0,1) # 12 month
nhis$prior.cancer <- ifelse(nhis$CANEV=="Yes",1,0)
nhis$lung.cancer.before <- ifelse(nhis$CNKIND14=="Mentioned",1,0)


nhis$cpd <- nhis$CIGSDAY 
nhis$cpd[nhis$cpd==97|nhis$cpd==98|nhis$cpd==99] <- NA
	
nhis$smoke.age.start <- nhis$SMKREG
nhis$smoke.age.start[nhis$smoke.age.start==96|nhis$smoke.age.start==97|nhis$smoke.age.start==98|nhis$smoke.age.start==99] <- NA
	
nhis$qtyears <- nhis$SMKQTY 	
nhis$qtyears[nhis$qtyears==97|nhis$qtyears==98|nhis$qtyears==99] <- NA
nhis$qtyears[!is.na(nhis$qtyears)&nhis$qtyears==0] <- 1/2
# TREAT FORMER QUITTERS WITH UNKNOWN AS CURRENT
nhis$former[nhis$SMKQTD=="Yes"|nhis$SMKQTD=="No"|nhis$SMKQTD=="Don't know"] <- 0
nhis$current[nhis$SMKQTD=="Yes"|nhis$SMKQTD=="No"|nhis$SMKQTD=="Don't know"] <- 1
nhis$qtyears[nhis$current==1] <- 0

#IF FORMER SMOKER WITH UNKNOWN SMOKE AGE START BUT KNOWN QUIT YEARS
nhis$smoke.age.start[is.na(nhis$smoke.age.start)&nhis$former==1&!is.na(nhis$qtyears)] <- (nhis$age-nhis$qtyears)[is.na(nhis$smoke.age.start)&nhis$former==1&!is.na(nhis$qtyears)] 
nhis$smkyears <- nhis$age-nhis$smoke.age.start-nhis$qtyears 
nhis$smkyears <- ifelse(!is.na(nhis$smkyears)&nhis$smkyears<0,1/12,nhis$smkyears) 
nhis$birthyear <- 2001-nhis$age

# APPLY EXCLUSION CRITERIA; EVER SMOKERS, NO HISTORY OF LUNG CANCER
sum(nhis$lung.cancer.before)
sum(nhis$never)
sum(nhis$unknown)

nhis <- subset(nhis, unknown==0&never==0&lung.cancer.before==0)

keep <- c("pid",
	"psu",
	"strata",
    "intyear",
    "birthyear",
    "year",
	"age",
	"asb",
	"bmi",
	"copd",
	"cpd",
	"current",
	"dust",
	"edu",
	"emp",
	"fam.cancer",
	"fam.lung.cancer",
	"fam.lung.trend",
	"fam.smoke.cancer",
	"fam.cancer.onset",
	"female",
	"former",
	"no.asthma",
	"prior.cancer",
	"qtyears",
	"smkyears",
	"smoke.age.start",
	"race"
)

nhis2001 <- nhis[,keep]
save(nhis2001, file="~/master/project/nhis/data/samadult/nhis2001.RData")
rm(list=ls())