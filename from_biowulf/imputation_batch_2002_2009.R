rm(list=ls(all=TRUE)) 
set.seed(11903808)

load(file="/home/cheunglc/lyg/nhis2002_09/nhis2002.Rdata")
load(file="/home/cheunglc/lyg/nhis2002_09/nhis2003.Rdata")
load(file="/home/cheunglc/lyg/nhis2002_09/nhis2004.Rdata")
load(file="/home/cheunglc/lyg/nhis2002_09/nhis2005.Rdata")
load(file="/home/cheunglc/lyg/nhis2002_09/nhis2006.Rdata")
load(file="/home/cheunglc/lyg/nhis2002_09/nhis2007.Rdata")
load(file="/home/cheunglc/lyg/nhis2002_09/nhis2008.Rdata")
load(file="/home/cheunglc/lyg/nhis2002_09/nhis2009.Rdata")
nhis2002_05 <- rbind(nhis2002,nhis2003,nhis2004,nhis2005)
nhis2002_05$strata <- 1000+nhis2002_05$strata
nhis2002_05$adj.wt <- nhis2002_05$wt/8 # ADJUSTED FOR POOLED ANALYSIS
nhis2002_05$adj.wt_mort <- nhis2002_05$wt_mort/8
nhis2006_09 <- rbind(nhis2006,nhis2007,nhis2008,nhis2009)
nhis2006_09$adj.wt <- nhis2006_09$wt/8 # ADJUSTED FOR POOLED ANALYSIS
nhis2006_09$adj.wt_mort <- nhis2006_09$wt_mort/8

nhis <- rbind(nhis2002_05,nhis2006_09)
source(file="/home/cheunglc/lyg/nhis_imputation_function.R")
save(nhis, file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_1.RData")

nhis <- rbind(nhis2002_05,nhis2006_09)
source(file="/home/cheunglc/lyg/nhis_imputation_function.R")
save(nhis, file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_2.RData")

nhis <- rbind(nhis2002_05,nhis2006_09)
source(file="/home/cheunglc/lyg/nhis_imputation_function.R")
save(nhis, file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_3.RData")

nhis <- rbind(nhis2002_05,nhis2006_09)
source(file="/home/cheunglc/lyg/nhis_imputation_function.R")
save(nhis, file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_4.RData")

nhis <- rbind(nhis2002_05,nhis2006_09)
source(file="/home/cheunglc/lyg/nhis_imputation_function.R")
save(nhis, file="/home/cheunglc/lyg/nhis2002_09/nhis_imputed_5.RData")
