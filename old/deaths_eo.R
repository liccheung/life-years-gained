# Fit different LCDRAT models for final choice
rm(list=ls(all=TRUE))

# Set up
library(survival)
library(coxph.risk)
library(glmnet)
library(pROC)
source("P:/bb/lrisk/prog/lifeyearsgained/cuts.R")
source("P:/bb/lrisk/prog/lifeyearsgained/kovalchik_LCC2014.R")

setwd("P:/bb/lrisk/prog/lifeyearsgained/")
load("P:/bb/lrisk/prog/jama/best_AIC_models.RData")


##################################
### External Validation - NHIS ###
##################################

source("P:/bb/lrisk/prog/lifeyearsgained/projectdeathRiskNHIS.R")
source("P:/bb/lrisk/prog/lifeyearsgained/EOfromNHIS_expanded2_LCC.R")

setwd("P:/bb/lrisk/other/lifeyearsgained/")
risklabels.eligible
save(risklabels.eligible, risklabels.ineligible, risklabels.all, risklabels.other.ineligible, file="risklabels2.RData")
save(calibration, calibration.risk.groups, file="calibration_risk10_2.RData")
write.table(calibration, file="calibration_overall10_2.csv",sep=",")
write.table(calibration.risk.groups, file="calibration_risk_groups10_2.csv",sep=",")

source("P:/bb/lrisk/prog/pop_discrimination.R")
write.table(aucs, file="LCDRAT_pop_auc.csv", sep=",")

save(summary.matrix.eligible, summary.matrix.ineligible5080, summary.matrix.ineligible2, summary.matrix.all5080, file="calibration_subgroups_2.RData")

