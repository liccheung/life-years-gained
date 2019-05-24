# IMPUTATION 
rm(list=ls(all=TRUE))
set.seed(11808901)
load(file="P:/bb/lrisk/other/nhis1997_2001/cpd.RData")

load(file="P:/bb/lrisk/other/lifeyearsgained/nhis/nhis.RData")
source(file="P:/bb/lrisk/prog/lifeyearsgained/nhis_imputation_function.death.R")
save(nhis, file="P:/bb/lrisk/other/lifeyearsgained/nhis/nhis_imputed1.RData")

set.seed(11808902)
load(file="P:/bb/lrisk/other/lifeyearsgained/nhis/nhis.RData")
source(file="P:/bb/lrisk/prog/lifeyearsgained/nhis_imputation_function.death.R")
save(nhis, file="P:/bb/lrisk/other/lifeyearsgained/nhis/nhis_imputed2.RData")

set.seed(11808903)
load(file="P:/bb/lrisk/other/lifeyearsgained/nhis/nhis.RData")
source(file="P:/bb/lrisk/prog/lifeyearsgained/nhis_imputation_function.death.R")
save(nhis, file="P:/bb/lrisk/other/lifeyearsgained/nhis/nhis_imputed3.RData")

set.seed(11808904)
load(file="P:/bb/lrisk/other/lifeyearsgained/nhis/nhis.RData")
source(file="P:/bb/lrisk/prog/lifeyearsgained/nhis_imputation_function.death.R")
save(nhis, file="P:/bb/lrisk/other/lifeyearsgained/nhis/nhis_imputed4.RData")

set.seed(11808905)
load(file="P:/bb/lrisk/other/lifeyearsgained/nhis/nhis.RData")
source(file="P:/bb/lrisk/prog/lifeyearsgained/nhis_imputation_function.death.R")
save(nhis, file="P:/bb/lrisk/other/lifeyearsgained/nhis/nhis_imputed5.RData")


