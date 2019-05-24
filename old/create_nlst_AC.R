rm(list=ls(all=TRUE))
library("haven")
library("lcrisks")
nlst_new <- read_sas("~/Desktop/Lung cancer/lrisk/sasfile/nlst/nlst364_prsn_20171023.sas7bdat")
load(file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/anlst_mod.RData")

persons <- data.frame(age=anlst$age,
                      female=anlst$female,
                      smkyears=anlst$smkyears,
                      qtyears=anlst$qtyears,
                      cpd=anlst$cpd,
                      race=as.numeric(as.character(anlst$race)),
                      emp=anlst$emp,
                      fam.lung.trend=anlst$fam.lung.trend,
                      bmi=anlst$bmi,
                      edu6=anlst$edu6)

missval <- ifelse(!is.na(persons$age) & !is.na(persons$female) & !is.na(persons$smkyears) & !is.na(persons$qtyears) &
                    !is.na(persons$cpd) & !is.na(persons$race) & !is.na(persons$emp) &
                    !is.na(persons$fam.lung.trend) & !is.na(persons$bmi) & !is.na(persons$edu6),0,1)

anlst$lcrat[missval==0] <- lcrisk(persons[missval==0,],5)[,5]/1000

quant.cp <- quantile(anlst$lcrat,c(.2,.4,.6,.8))
anlst$lcrat.cat <- ifelse(anlst$lcrat<quant.cp[1],1,
                          ifelse(anlst$lcrat<quant.cp[2],2,
                                 ifelse(anlst$lcrat<quant.cp[3],3,
                                        ifelse(anlst$lcrat<quant.cp[4],4,5))))

nlst_new2 <- nlst_new[c("pid","rand_month","rand_year","de_stag")]
nlst <- merge(anlst,nlst_new2,by="pid")
nlst <- data.frame(
  pid = nlst$pid,
  screen_group = nlst$screen_group,
  age = nlst$age,
  age.stopped = nlst$age.stopped,
  asb=nlst$asb,
  bmi = nlst$bmi,
  copd = nlst$copd,
  cpd = nlst$cpd,
  current = nlst$current,
  dust = nlst$dust,
  emp = nlst$emp,
  fam.cancer = nlst$fam.cancer,
  fam_lung_trend = nlst$fam.lung.trend,
  female = nlst$female,
  former = nlst$former,
  event_lung_cancer = nlst$event.lung.cancer,
  lung_cancer_death = nlst$lung.cancer.death,
  no.asthma = nlst$no.asthma,
  death = nlst$death,
  death.years = nlst$death.years,
  packyears = nlst$packyears,
  pneu = nlst$pneu,
  diab = nlst$diab,
  hypertension = nlst$hypertension,
  stroke = nlst$stroke,
  prior_cancer = nlst$prior.cancer,
  qtyears = nlst$qtyears,
  smkyears = nlst$smkyears,
  years_followed = nlst$years.followed,
  lung.comorbid.count = nlst$lung.comorbid.count,
  race = nlst$race,
  stage = nlst$stage,
  truefalse_scrnres_ly0 = nlst$truefalse_scrnres_ly0,
  truefalse_scrnres_ly1 = nlst$truefalse_scrnres_ly1,
  truefalse_scrnres_ly2 = nlst$truefalse_scrnres_ly2,
  eligible = nlst$eligible,
  lc_years = nlst$mhqcohort_incdays/365.25,
  chd = nlst$chd,
  angina = nlst$angina,
  heartattack = nlst$heartattack,
  heartdisease = nlst$heartdisease,
  bron = nlst$bron,
  kidney = nlst$kidney,
  liver = nlst$liver,
  speceq = nlst$speceq,
  lcrat5 = nlst$lcrat,
  lcrat5_cat = nlst$lcrat.cat,
  lcd5 = nlst$lcd5,
  ely = nlst$ely,
  ely_ct = nlst$ely_ct,
  lyg = nlst$lyg,
  rand_month = nlst$rand_month,
  rand_year = nlst$rand_year,
  de_stag = nlst$de_stag
)
save(nlst,file="~/Desktop/Lung cancer/lrisk/sasfile/nlst/nlst_lcc.RData")