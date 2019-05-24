library(coxph.risk)
setwd("P:/bb/lrisk/other/lifeyearsgained/nhis")

load("P:/bb/lrisk/other/lifeyearsgained/nhis/nhis_imputed1.RData")
nhis$risk1 <- risk.kovalchik(0, 5, nhis, LCDRAT, cox.death.original)
nhis$risk2 <- risk.kovalchik(0, 5, nhis, cox.death.original, LCDRAT)
nhis$risk <- nhis$risk1+nhis$risk2
nhis$wt <- nhis$wt_mort/5
save(nhis, file="nhis_risk1.RData")

load("P:/bb/lrisk/other/lifeyearsgained/nhis/nhis_imputed2.RData")
nhis$risk1 <- risk.kovalchik(0, 5, nhis, LCDRAT, cox.death.original)
nhis$risk2 <- risk.kovalchik(0, 5, nhis, cox.death.original, LCDRAT)
nhis$risk <- nhis$risk1+nhis$risk2
nhis$wt <- nhis$wt_mort/5
save(nhis, file="nhis_risk2.RData")

load("P:/bb/lrisk/other/lifeyearsgained/nhis/nhis_imputed3.RData")
nhis$risk1 <- risk.kovalchik(0, 5, nhis, LCDRAT, cox.death.original)
nhis$risk2 <- risk.kovalchik(0, 5, nhis, cox.death.original, LCDRAT)
nhis$risk <- nhis$risk1+nhis$risk2
nhis$wt <- nhis$wt_mort/5
save(nhis, file="nhis_risk3.RData")

load("P:/bb/lrisk/other/lifeyearsgained/nhis/nhis_imputed4.RData")
nhis$risk1 <- risk.kovalchik(0, 5, nhis, LCDRAT, cox.death.original)
nhis$risk2 <- risk.kovalchik(0, 5, nhis, cox.death.original, LCDRAT)
nhis$risk <- nhis$risk1+nhis$risk2
nhis$wt <- nhis$wt_mort/5
save(nhis, file="nhis_risk4.RData")

load("P:/bb/lrisk/other/lifeyearsgained/nhis/nhis_imputed5.RData")
nhis$risk1 <- risk.kovalchik(0, 5, nhis, LCDRAT, cox.death.original)
nhis$risk2 <- risk.kovalchik(0, 5, nhis, cox.death.original, LCDRAT)
nhis$risk <- nhis$risk1+nhis$risk2
nhis$wt <- nhis$wt_mort/5
save(nhis, file="nhis_risk5.RData")
