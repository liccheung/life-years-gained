rm(list=ls())
setwd("~/Desktop/Lung cancer/lrisk/prog/lifeyearsgained")
load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis97_01.RData")
nhis <- nhis97_01

names(nhis)[names(nhis)=="AGE_P"] <- "age"
names(nhis)[names(nhis)=="WTFA_SA"]  <- "wt"
names(nhis)[names(nhis)=="stratum"]  <- "strata"
names(nhis)[names(nhis)=="NHIS_YR"] <- "year"

nhis$intyear <- nhis$year+.125*(nhis$INTV_QRT==1)+.375*(nhis$INTV_QRT==2)+
                .625*(nhis$INTV_QRT==3)+.875*(nhis$INTV_QRT==4)
nhis$female <- ifelse(nhis$sex==2,1,0)

nhis$edu6 <- ifelse(nhis$educ==0|nhis$educ==1|nhis$educ==2|nhis$educ==3|nhis$educ==4|nhis$educ==5|nhis$educ==6|nhis$educ==7|nhis$educ==8|nhis$educ==9|nhis$educ==10|nhis$educ==11|nhis$educ==12,1,
                    ifelse(nhis$educ==13,2,
                           ifelse(nhis$educ==16|nhis$educ==14,3,
                                  ifelse(nhis$educ==15,4,
                                         ifelse(nhis$educ==17|nhis$educ==18,5,
                                                ifelse(nhis$educ==98|nhis$educ==97|nhis$educ==99,NA,6))))))

nhis$bmi[nhis$bmi==99.99] <- NA

nhis$current <- ifelse(nhis$SMKSTAT1==1|nhis$SMKSTAT1==4,1,0)										  
nhis$former <- ifelse(nhis$SMKSTAT1==2,1,0)											  
nhis$never <- ifelse(nhis$SMKSTAT1==3,1,0)											  
nhis$unknown <-  ifelse(nhis$SMKSTAT1==9,1,0)	

nhis$fam.lung.trend <- ifelse(nhis$FHFTYP14==1,1,0)+ifelse(nhis$FHMTYP14==1,1,0)
nhis$fam.lung.trend[is.na(nhis$fam.lung.trend)] <- 0

nhis$prior.cancer <- ifelse(nhis$canev==1,1,0)
nhis$lung.cancer.before <- ifelse(!is.na(nhis$CNKIND14) & nhis$CNKIND14==1,1,0)
nhis$hypertension <- ifelse(nhis$hypev==1,1,0)
nhis$chd <- ifelse(nhis$chdev==1,1,0)
nhis$angina <- ifelse(nhis$angev==1,1,0)
nhis$heartattack <- ifelse(nhis$miev==1,1,0)	
nhis$heartdisease <- ifelse(nhis$hrtev==1,1,0)	
nhis$stroke <- ifelse(nhis$strev==1,1,0)	
nhis$emp <- ifelse(nhis$ephev==1,1,0)
nhis$diab <- ifelse(nhis$dibev==1,1,0)	 
nhis$bron <- ifelse(nhis$cbrchyr==1,1,0)	
nhis$kidney <- ifelse(nhis$kidwkyr==1,1,0)	
nhis$liver <- ifelse(nhis$livyr==1,1,0)	
nhis$speceq <- ifelse(nhis$speceq==1,1,0)

nhis$cpd <- nhis$cigsday
nhis$cpd[nhis$cpd==97|nhis$cpd==98|nhis$cpd==99] <- NA

nhis$smoke.age.start <- nhis$smkreg
nhis$smoke.age.start[nhis$smoke.age.start==96|nhis$smoke.age.start==97|nhis$smoke.age.start==98|nhis$smoke.age.start==99] <- NA

nhis$qtyears <- nhis$smkqty
nhis$qtyears[nhis$qtyears==97|nhis$qtyears==98|nhis$qtyears==99] <- NA
nhis$qtyears[!is.na(nhis$qtyears)&nhis$qtyears==0] <- 1/2

# TREAT QUITTERS IN PAST YEAR AS CURRENT
nhis$former[nhis$smkqtd==1|nhis$smkqtd==2|nhis$smkqtd==9] <- 0
nhis$current[nhis$smkqtd==1|nhis$smkqtd==2|nhis$smkqtd==9] <- 1
nhis$qtyears[nhis$current==1] <- 0

nhis$smkyears <- nhis$age-nhis$smoke.age.start-nhis$qtyears 
nhis$smkyears <- ifelse(!is.na(nhis$smkyears)&nhis$smkyears<0,1/12,nhis$smkyears)

# APPLY EXCLUSION CRITERIA; EVER SMOKERS, NO HISTORY OF LUNG CANCER
sum(nhis$lung.cancer.before)
sum(nhis$never)
sum(nhis$unknown)
sum(is.na(nhis$smoke.age.start))

nhis <- subset(nhis, unknown==0&never==0&lung.cancer.before==0&
                    is.na(smoke.age.start)==0)#&is.na(deathage)==0)


keep <- c("pid",
          "psu",
          "strata",
          "wt",
          "intyear",
          "year",
          "age",
          "bmi",
          "cpd",
          "current",
          "edu6",
          "fam.lung.trend",
          "hypertension",
          "chd",
          "angina",
          "heartattack",
          "heartdisease",
          "stroke",
          "emp",
          "diab",
          "bron",
          "kidney",
          "liver",
          "speceq",
          "female",
          "former",
          "prior.cancer",
          "qtyears",
          "smkyears",
          "smoke.age.start")
          #"deathage",
          #"wt_mort",
          #"wt_mort5")

nhis <- nhis[,keep]

# MERGE AND SUMMARIZE SAMPLE ADULT FILES
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/nhis1997.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/nhis1998.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/nhis1999.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/nhis2000.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/nhis2001.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/mortality2011/mortality2011.RData")

oldnhis <- rbind(nhis1997, nhis1998, nhis1999, nhis2000, nhis2001)

mortality2011$year <- NULL
nhis <- merge(mortality2011, nhis, by="pid") # MERGE THOSE WITH KNOWN MORTALITY OUTCOMES
oldnhis <- subset(oldnhis,select=c("pid","race"))
nhis <- merge(nhis, oldnhis, by="pid")

nhis$race <- ifelse(is.na(nhis$race),"Missing",nhis$race)
nhis$race <- factor(nhis$race, 
                    lev=c("Non-Hispanic White", 
                          "Non-Hispanic Black", 
                          "Hispanic", 
                          "Asian", 
                          "American Indian/Alaskan Native",
                          "Missing"),exclude=NULL,order=TRUE)

nhis$deathage <- nhis$deathyear-nhis$intyear+nhis$age

save(nhis,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/cleaned.nhis97_01.RData")