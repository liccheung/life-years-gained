rm(list=ls(all=TRUE)) 
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis.imputed.1.RData")
nhis97_01 <- nhis
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_2002_05_1.RData")
nhis02_05 <- nhis

nhis97_01 <- nhis97_01[which(names(nhis97_01) %in% names(nhis02_05)==TRUE)]
nhis02_05 <- nhis02_05[which(names(nhis02_05) %in% names(nhis97_01)==TRUE)]
nhis <- rbind(nhis97_01,nhis02_05)
nhis$adj.wt <- nhis$wt/9 # ADJUSTED FOR POOLED ANALYSIS
nhis$adj.wt_mort <- nhis$wt_mort/9
save(nhis, file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_1997_05_1.RData")



rm(list=ls(all=TRUE)) 
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis.imputed.2.RData")
nhis97_01 <- nhis
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_2002_05_2.RData")
nhis02_05 <- nhis

nhis97_01 <- nhis97_01[which(names(nhis97_01) %in% names(nhis02_05)==TRUE)]
nhis02_05 <- nhis02_05[which(names(nhis02_05) %in% names(nhis97_01)==TRUE)]
nhis <- rbind(nhis97_01,nhis02_05)
nhis$adj.wt <- nhis$wt/9 # ADJUSTED FOR POOLED ANALYSIS
nhis$adj.wt_mort <- nhis$wt_mort/9
save(nhis, file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_1997_05_2.RData")



rm(list=ls(all=TRUE)) 
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis.imputed.3.RData")
nhis97_01 <- nhis
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_2002_05_3.RData")
nhis02_05 <- nhis

nhis97_01 <- nhis97_01[which(names(nhis97_01) %in% names(nhis02_05)==TRUE)]
nhis02_05 <- nhis02_05[which(names(nhis02_05) %in% names(nhis97_01)==TRUE)]
nhis <- rbind(nhis97_01,nhis02_05)
nhis$adj.wt <- nhis$wt/9 # ADJUSTED FOR POOLED ANALYSIS
nhis$adj.wt_mort <- nhis$wt_mort/9
save(nhis, file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_1997_05_3.RData")



rm(list=ls(all=TRUE)) 
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis.imputed.4.RData")
nhis97_01 <- nhis
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_2002_05_4.RData")
nhis02_05 <- nhis

nhis97_01 <- nhis97_01[which(names(nhis97_01) %in% names(nhis02_05)==TRUE)]
nhis02_05 <- nhis02_05[which(names(nhis02_05) %in% names(nhis97_01)==TRUE)]
nhis <- rbind(nhis97_01,nhis02_05)
nhis$adj.wt <- nhis$wt/9 # ADJUSTED FOR POOLED ANALYSIS
nhis$adj.wt_mort <- nhis$wt_mort/9
save(nhis, file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_1997_05_4.RData")



rm(list=ls(all=TRUE)) 
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis.imputed.5.RData")
nhis97_01 <- nhis
load(file="~/Desktop/Lung cancer/lrisk/other/nhis2002_09/nhis_imputed_2002_05_5.RData")
nhis02_05 <- nhis

nhis97_01 <- nhis97_01[which(names(nhis97_01) %in% names(nhis02_05)==TRUE)]
nhis02_05 <- nhis02_05[which(names(nhis02_05) %in% names(nhis97_01)==TRUE)]
nhis <- rbind(nhis97_01,nhis02_05)
nhis$adj.wt <- nhis$wt/9 # ADJUSTED FOR POOLED ANALYSIS
nhis$adj.wt_mort <- nhis$wt_mort/9
save(nhis, file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_1997_05_5.RData")