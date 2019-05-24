rm(list=ls(all=TRUE))
set.seed(11808901)
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/cleaned.nhis97_01.RData")
nhis97_01 <- nhis
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis2002.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis2003.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis2004.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis2005.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis2006.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis2007.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis2008.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis2009.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2013_15/nhis2013.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2013_15/nhis2014.RData")
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2013_15/nhis2015.RData")

nhis97_01$edu6 <- as.numeric(nhis97_01$edu6)
nhis97_01$edu6[nhis97_01$edu6==7] <- NA
nhis2002 <- nhis2002[,which(names(nhis2002) %in% names(nhis2009)==TRUE)]
nhis2003 <- nhis2003[,which(names(nhis2003) %in% names(nhis2009)==TRUE)]
nhis2004 <- nhis2004[,which(names(nhis2004) %in% names(nhis2009)==TRUE)]
nhis2005 <- nhis2005[,which(names(nhis2005) %in% names(nhis2009)==TRUE)]

nhis02_09 <- rbind(nhis2002,nhis2003,nhis2004,nhis2005,nhis2006,nhis2007,nhis2008,nhis2009)
nhis97_01 <- nhis97_01[which(names(nhis97_01) %in% names(nhis02_09)==TRUE)]
nhis02_09 <- nhis02_09[which(names(nhis02_09) %in% names(nhis97_01)==TRUE)]

train <- rbind(subset(nhis97_01,year %in% c(1997,1999,2001)),
               subset(nhis02_09,year %in% c(2003,2005,2007,2009)))
valid <- rbind(subset(nhis97_01,year %in% c(1998,2000)),
               subset(nhis02_09,year %in% c(2002,2004,2006,2008)))
proj <- rbind(nhis2013,nhis2014,nhis2015)

train <- subset(train,40<=train$age&train$age<=84&(train$current==1|train$former==1)&train$lung.cancer.before==0 & 
                  !(is.na(train$smoke.age.start)|train$smoke.age.start<0))

ncov1 <- c(table(40<=train$age&train$age<50)[2],
  table(50<=train$age&train$age<60)[2],
  table(60<=train$age&train$age<70)[2],
  table(70<=train$age&train$age<=79)[2],
  table(80<=train$age&train$age<=84)[2],
  table(train$female==1),
  table(train$race)[c(5,4,3)],
  sum(table(train$race)[c(1,2)]),
  table(train$bmi<=18.5)[2],
  table(18.5<train$bmi&train$bmi<=20)[2],
  table(20<train$bmi&train$bmi<=25)[2],
  table(25<train$bmi&train$bmi<=30)[2],
  table(30<train$bmi&train$bmi<=35)[2],
  table(35<train$bmi&train$bmi)[2],
  table(train$current==1)[2],
  table(train$current==1 & train$cpd<20)[2],
  table(train$current==1 & 20<=train$cpd&train$cpd<30)[2],
  table(train$current==1 & 30<=train$cpd&train$cpd<40)[2],  
  table(train$current==1 & 40<=train$cpd)[2],  
  table(train$current==1 & train$smkyears<20)[2],
  table(train$current==1 & 20<=train$smkyears&train$smkyears<30)[2],
  table(train$current==1 & 30<=train$smkyears&train$smkyears<40)[2],  
  table(train$current==1 & 40<=train$smkyears)[2],  
  table(train$former==1)[2],
  table(train$former==1 & train$cpd<20)[2],
  table(train$former==1 & 20<=train$cpd&train$cpd<30)[2],
  table(train$former==1 & 30<=train$cpd&train$cpd<40)[2],  
  table(train$former==1 & 40<=train$cpd)[2],
  table(train$former==1 & train$smkyears<20)[2],
  table(train$former==1 & 20<=train$smkyears&train$smkyears<30)[2],
  table(train$former==1 & 30<=train$smkyears&train$smkyears<40)[2],  
  table(train$former==1 & 40<=train$smkyears)[2],  
  table(train$former==1 & train$qtyears<=5)[2],
  table(train$former==1 & 5<train$qtyears&train$qtyears<=10)[2],
  table(train$former==1 & 10<train$qtyears&train$qtyears<=15)[2],
  table(train$former==1 & 15<train$qtyears&train$qtyears<=20)[2],  
  table(train$former==1 & 20<train$qtyears&train$qtyears<=30)[2],  
  table(train$former==1 & 30<train$qtyears)[2],  
  table(train$angina),
  table(train$bron),
  table(train$hypertension),
  table(train$heartdisease),
  table(train$chd),
  table(train$stroke),
  table(train$heartattack),
  table(train$prior.cancer),
  table(train$kidney),
  table(train$diab),
  table(train$liver),
  table(train$emp),
  table(train$speceq),
  table(train$fam.lung.trend))

valid <- subset(valid,40<=valid$age&valid$age<=84&(valid$current==1|valid$former==1)&valid$lung.cancer.before==0 & 
                  !(is.na(valid$smoke.age.start)|valid$smoke.age.start<0))

ncov2 <- c(table(40<=valid$age&valid$age<50)[2],
          table(50<=valid$age&valid$age<60)[2],
          table(60<=valid$age&valid$age<70)[2],
          table(70<=valid$age&valid$age<=79)[2],
          table(80<=valid$age&valid$age<=84)[2],
          table(valid$female==1),
          table(valid$race)[c(5,4,3)],
          sum(table(valid$race)[c(1,2)]),
          table(valid$bmi<=18.5)[2],
          table(18.5<valid$bmi&valid$bmi<=20)[2],
          table(20<valid$bmi&valid$bmi<=25)[2],
          table(25<valid$bmi&valid$bmi<=30)[2],
          table(30<valid$bmi&valid$bmi<=35)[2],
          table(35<valid$bmi&valid$bmi)[2],
          table(valid$current==1)[2],
          table(valid$current==1 & valid$cpd<20)[2],
          table(valid$current==1 & 20<=valid$cpd&valid$cpd<30)[2],
          table(valid$current==1 & 30<=valid$cpd&valid$cpd<40)[2],  
          table(valid$current==1 & 40<=valid$cpd)[2],  
          table(valid$current==1 & valid$smkyears<20)[2],
          table(valid$current==1 & 20<=valid$smkyears&valid$smkyears<30)[2],
          table(valid$current==1 & 30<=valid$smkyears&valid$smkyears<40)[2],  
          table(valid$current==1 & 40<=valid$smkyears)[2],  
          table(valid$former==1)[2],
          table(valid$former==1 & valid$cpd<20)[2],
          table(valid$former==1 & 20<=valid$cpd&valid$cpd<30)[2],
          table(valid$former==1 & 30<=valid$cpd&valid$cpd<40)[2],  
          table(valid$former==1 & 40<=valid$cpd)[2],
          table(valid$former==1 & valid$smkyears<20)[2],
          table(valid$former==1 & 20<=valid$smkyears&valid$smkyears<30)[2],
          table(valid$former==1 & 30<=valid$smkyears&valid$smkyears<40)[2],  
          table(valid$former==1 & 40<=valid$smkyears)[2],  
          table(valid$former==1 & valid$qtyears<=5)[2],
          table(valid$former==1 & 5<valid$qtyears&valid$qtyears<=10)[2],
          table(valid$former==1 & 10<valid$qtyears&valid$qtyears<=15)[2],
          table(valid$former==1 & 15<valid$qtyears&valid$qtyears<=20)[2],  
          table(valid$former==1 & 20<valid$qtyears&valid$qtyears<=30)[2],  
          table(valid$former==1 & 30<valid$qtyears)[2],  
          table(valid$angina),
          table(valid$bron),
          table(valid$hypertension),
          table(valid$heartdisease),
          table(valid$chd),
          table(valid$stroke),
          table(valid$heartattack),
          table(valid$prior.cancer),
          table(valid$kidney),
          table(valid$diab),
          table(valid$liver),          
          table(valid$emp),
          table(valid$speceq),
          table(valid$fam.lung.trend))

proj <- subset(proj,40<=proj$age&proj$age<=84&(proj$current==1|proj$former==1)&proj$lung.cancer.before==0 & 
                  !(is.na(proj$smoke.age.start)|proj$smoke.age.start<0))

ncov3 <- c(table(40<=proj$age&proj$age<50)[2],
          table(50<=proj$age&proj$age<60)[2],
          table(60<=proj$age&proj$age<70)[2],
          table(70<=proj$age&proj$age<=79)[2],
          table(80<=proj$age&proj$age<=84)[2],
          table(proj$female==1),
          table(proj$race)[c(5,4,3)],
          sum(table(proj$race)[c(1,2)]),
          table(proj$bmi<=18.5)[2],
          table(18.5<proj$bmi&proj$bmi<=20)[2],
          table(20<proj$bmi&proj$bmi<=25)[2],
          table(25<proj$bmi&proj$bmi<=30)[2],
          table(30<proj$bmi&proj$bmi<=35)[2],
          table(35<proj$bmi&proj$bmi)[2],
          table(proj$current==1)[2],
          table(proj$current==1 & proj$cpd<20)[2],
          table(proj$current==1 & 20<=proj$cpd&proj$cpd<30)[2],
          table(proj$current==1 & 30<=proj$cpd&proj$cpd<40)[2],  
          table(proj$current==1 & 40<=proj$cpd)[2],  
          table(proj$current==1 & proj$smkyears<20)[2],
          table(proj$current==1 & 20<=proj$smkyears&proj$smkyears<30)[2],
          table(proj$current==1 & 30<=proj$smkyears&proj$smkyears<40)[2],  
          table(proj$current==1 & 40<=proj$smkyears)[2],  
          table(proj$former==1)[2],
          table(proj$former==1 & proj$cpd<20)[2],
          table(proj$former==1 & 20<=proj$cpd&proj$cpd<30)[2],
          table(proj$former==1 & 30<=proj$cpd&proj$cpd<40)[2],  
          table(proj$former==1 & 40<=proj$cpd)[2],
          table(proj$former==1 & proj$smkyears<20)[2],
          table(proj$former==1 & 20<=proj$smkyears&proj$smkyears<30)[2],
          table(proj$former==1 & 30<=proj$smkyears&proj$smkyears<40)[2],  
          table(proj$former==1 & 40<=proj$smkyears)[2],  
          table(proj$former==1 & proj$qtyears<=5)[2],
          table(proj$former==1 & 5<proj$qtyears&proj$qtyears<=10)[2],
          table(proj$former==1 & 10<proj$qtyears&proj$qtyears<=15)[2],
          table(proj$former==1 & 15<proj$qtyears&proj$qtyears<=20)[2],  
          table(proj$former==1 & 20<proj$qtyears&proj$qtyears<=30)[2],  
          table(proj$former==1 & 30<proj$qtyears)[2],
          table(proj$angina),
          table(proj$bron),
          table(proj$hypertension),
          table(proj$heartdisease),
          table(proj$chd),
          table(proj$stroke),
          table(proj$heartattack),
          table(proj$prior.cancer),
          table(proj$kidney),
          table(proj$diab),
          table(proj$liver),
          table(proj$emp),
          table(proj$speceq),
          table(proj$fam.lung.trend))

res <- cbind(ncov1,100*ncov1/nrow(train),ncov2,100*ncov2/nrow(valid),ncov3,100*ncov3/nrow(proj))
write.table(res,"~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/cov_freq.csv",sep=",")


#Missings
table(is.na(train$race))
table(is.na(train$race))/nrow(train)
table(is.na(train$bmi))
table(is.na(train$bmi))/nrow(train)
table(is.na(train$qtyears))
table(is.na(train$qtyears))/nrow(train)
table(train$current & is.na(train$cpd))
table(train$current & is.na(train$cpd))/nrow(train)
table(train$former & is.na(train$cpd))
table(train$former & is.na(train$cpd))/nrow(train)
table(train$current & is.na(train$smkyears))
table(train$current & is.na(train$smkyears))/nrow(train)
table(train$former & is.na(train$smkyears))
table(train$former & is.na(train$smkyears))/nrow(train)


table(is.na(valid$race))
table(is.na(valid$race))/nrow(valid)
table(is.na(valid$bmi))
table(is.na(valid$bmi))/nrow(valid)
table(is.na(valid$qtyears))
table(is.na(valid$qtyears))/nrow(valid)
table(valid$current & is.na(valid$cpd))
table(valid$current & is.na(valid$cpd))/nrow(valid)
table(valid$former & is.na(valid$cpd))
table(valid$former & is.na(valid$cpd))/nrow(valid)
table(valid$current & is.na(valid$smkyears))
table(valid$current & is.na(valid$smkyears))/nrow(valid)
table(valid$former & is.na(valid$smkyears))
table(valid$former & is.na(valid$smkyears))/nrow(valid)


table(is.na(proj$race))
table(is.na(proj$race))/nrow(proj)
table(is.na(proj$bmi))
table(is.na(proj$bmi))/nrow(proj)
table(is.na(proj$qtyears))
table(is.na(proj$qtyears))/nrow(proj)
table(proj$current & is.na(proj$cpd))
table(proj$current & is.na(proj$cpd))/nrow(proj)
table(proj$former & is.na(proj$cpd))
table(proj$former & is.na(proj$cpd))/nrow(proj)
table(proj$current & is.na(proj$smkyears))
table(proj$current & is.na(proj$smkyears))/nrow(proj)
table(proj$former & is.na(proj$smkyears))
table(proj$former & is.na(proj$smkyears))/nrow(proj)


table(train$edu6,useNA="always")
table(train$edu6,useNA="always")/nrow(train)
table(valid$edu6,useNA="always")
table(valid$edu6,useNA="always")/nrow(valid)
table(proj$EDUC1,useNA="always")
table(proj$edu6,useNA="always")/nrow(proj)