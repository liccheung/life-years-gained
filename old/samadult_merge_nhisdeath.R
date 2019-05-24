# MERGE AND SUMMARIZE SAMPLE ADULT FILES
rm(list=ls())
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/nhis1997.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/nhis1998.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/nhis1999.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/nhis2000.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis1997_2001/nhis2001.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/stephanie/nhis/data/mortality/mortality.RData")

nhis <- rbind(nhis1997, nhis1998, nhis1999, nhis2000, nhis2001)
mortality$year <- NULL
nhis <- merge(mortality, nhis, by="pid") # MERGE THOSE WITH KNOWN MORTALITY OUTCOMES
nhis$deathyears <- nhis$deathyear-nhis$intyear
nhis <- nhis[nhis$deathyears>0,] # ALIVE AT INTERVIEW
nrow(nhis)
nhis <- nhis[!is.na(nhis$smoke.age.start),]

nhis$five.year <- ifelse(nhis$deathyears<5&nhis$died==1,1,0)

nhis$edu <- factor(nhis$edu, 
							levels = c("Less than 8 years", 
											"8-11 years", 
											"12 years or completed high school", 
											"Post-high school training other than college", 
											"Some college", 
											"College graduate", "Post graduate", "Missing"),
							labels = c(1:7,"Missing"),order=TRUE
)

nhis$edu6 <- factor(ifelse(nhis$edu<=2,1,nhis$edu), lev=c(1,3:8),
								lab= c("<High School",
											"High School", 
											"Post-High School Training", 
											"Some College", 
											"College", "Post-Graduate","Missing"), order=TRUE)

nhis$race <- ifelse(is.na(nhis$race),"Missing",nhis$race)
nhis$race <- factor(nhis$race, 
								lev=c("Non-Hispanic White", 
									"Non-Hispanic Black", 
									"Hispanic", 
									"Asian", 
									"American Indian/Alaskan Native",
									"Missing"),exclude=NULL,order=TRUE)

save(nhis, file="P:/bb/lrisk/other/lifeyearsgained/nhis/nhis.RData")
