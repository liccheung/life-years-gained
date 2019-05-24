# Preparation of 1997-2001 NHIS mortality-linked outcomes
# Datasets include all NHIS respondents, not only those completing the sample adult file

mortality1997 <- read.csv(file="~/master/project/nhis/data/mortality/nhis_mortality_1997.csv")
mortality1998 <- read.csv(file="~/master/project/nhis/data/mortality/nhis_mortality_1998.csv")
mortality1999 <- read.csv(file="~/master/project/nhis/data/mortality/nhis_mortality_1999.csv")
mortality2000 <- read.csv(file="~/master/project/nhis/data/mortality/nhis_mortality_2000.csv")
mortality2001 <- read.csv(file="~/master/project/nhis/data/mortality/nhis_mortality_2001.csv")

mortality1997$year <- 1997
mortality1998$year <- 1998
mortality1999$year <- 1999
mortality2000$year <- 2000
mortality2001$year <- 2001

mortality1997$DODYEAR <- as.numeric(as.character(mortality1997$DODYEAR))
mortality1998$DODYEAR <- as.numeric(as.character(mortality1998$DODYEAR))
mortality1999$DODYEAR <- as.numeric(as.character(mortality1999$DODYEAR))
mortality2000$DODYEAR <- as.numeric(as.character(mortality2000$DODYEAR))
mortality2001$DODYEAR <- as.numeric(as.character(mortality2001$DODYEAR))

mortality <- rbind(
	mortality1997, mortality1998, mortality1999, mortality2000, mortality2001
)

# SUBSET TO ELIGIBLE ADULTS
mortality <- subset(mortality, ELIGSTAT=="Eligible")
mortality$died <- ifelse(mortality$MORTSTAT=="Assumed deceased",1,0)
mortality$deathyear <- ifelse(is.na(mortality$DODYEAR),2007,mortality$DODYEAR) # ALIVE AT END OF Dec 31, 2006
mortality$dodqtr <- c(2,0,1,3,4)[mortality$DODQTR]
mortality$deathyear <- mortality$deathyear+mortality$dodqtr*1.5/12
mortality$lung.cancer.death <- ifelse(!is.na(mortality$UCOD_113)&mortality$UCOD_113==27,1,0) # lung cancer (C33-C34) 
names(mortality)[names(mortality)=="SA_WGT_NEW"]  <- "wt_mort"
names(mortality)[names(mortality)=="PUBLICID_2"]  <- "pid"

keep <- c(
	"pid",
	"wt_mort",
	"year",
	"died",
	"lung.cancer.death",
	"deathyear"
)

mortality <- mortality[,keep]

save(mortality, file="~/master/project/nhis/data/mortality/mortality.RData")
