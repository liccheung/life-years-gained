#marginal selection of USPSTF, Risk-based, and Life-years gained
rm(list=ls(all=TRUE))
library(lcrisks)
library(survey)
library(data.table)

inf_eo1 <- function(eo1,eo2,eo3,eo4,eo5){
  ave_eo <- (eo1+eo2+eo3+eo4+eo5)/5  #average estimate
  vars.bar <- (eo1[2,]+eo2[2,]+eo3[2,]+eo4[2,]+eo5[2,])/5  #within imputation variance
  v.impute <- ((eo1[1,]-ave_eo[1,])^2+(eo2[1,]-ave_eo[1,])^2+(eo3[1,]-ave_eo[1,])^2+(eo4[1,]-ave_eo[1,])^2+(eo5[1,]-ave_eo[1,])^2)/4  #between-imputation variance
  var.ave_eo <- (vars.bar+(1+1/5)*v.impute)  #variance for average estimate
  #confidence intervals for average estimate
  lower <- ave_eo[1,]-1.96*sqrt(var.ave_eo)
  upper <- ave_eo[1,]+1.96*sqrt(var.ave_eo)
  ave_eo <- rbind(ave_eo[1,],var.ave_eo,lower,upper)
  rownames(ave_eo) <- c("Ave","Var","LCL","UCL")
  return(ave_eo)
}


inf_eo2 <- function(eo1,eo2,eo3,eo4,eo5){
  ave_eo <- (eo1+eo2+eo3+eo4+eo5)/5  #average proportion
  vars.bar <- (eo1[2,]+eo2[2,]+eo3[2,]+eo4[2,]+eo5[2,])/5  #within imputation variance
  v.impute <- ((eo1[1,]-ave_eo[1,])^2+(eo2[1,]-ave_eo[1,])^2+(eo3[1,]-ave_eo[1,])^2+(eo4[1,]-ave_eo[1,])^2+(eo5[1,]-ave_eo[1,])^2)/4  #between-imputation variance
  var.ave_eo <- (vars.bar+(1+1/5)*v.impute)  #variance for average proportion
  logv <- var.ave_eo/ave_eo[1,]^2  #variance of log of average proportion
  #confidence intervals for average proportion
  lower <- exp(log(ave_eo[1,])-1.96*sqrt(logv))
  upper <- exp(log(ave_eo[1,])+1.96*sqrt(logv))
  ave_eo <- rbind(ave_eo[1,],var.ave_eo,lower,upper)
  rownames(ave_eo) <- c("Ave","Var","LCL","UCL")
  return(ave_eo)
}


inf_diff <- function(diff1,diff2,diff3,diff4,diff5){
  ave_diff <- (diff1+diff2+diff3+diff4+diff5)/5  #average diff
  vars.bar <-  ave_diff[2,]  #within imputation variance
  v.impute <- ((diff1[1,]-ave_diff[1,])^2+(diff2[1,]-ave_diff[1,])^2+(diff3[1,]-ave_diff[1,])^2+(diff4[1,]-ave_diff[1,])^2+(diff5[1,]-ave_diff[1,])^2)/4  #between-imputation variance
  var.ave_diff <- (vars.bar+(1+1/5)*v.impute)  #variance for average diff
  #confidence intervals for average E/O
  lower <- ave_diff[1,]-1.96*sqrt(var.ave_diff)
  upper <- ave_diff[1,]+1.96*sqrt(var.ave_diff)
  ave_diff <- rbind(ave_diff[1,],lower,upper,2*(1-pnorm(ave_diff[1,]/sqrt(var.ave_diff))))
  rownames(ave_diff) <- c("mean","lower","upper","2-sided pvalue")
  return(ave_diff)
}

inf_diff2 <- function(diff1,diff2,diff3,diff4,diff5){
  ave_diff <- (diff1+diff2+diff3+diff4+diff5)/5  #average diff
  vars.bar <-  ave_diff[2]  #within imputation variance
  v.impute <- ((diff1[1]-ave_diff[1])^2+(diff2[1]-ave_diff[1])^2+(diff3[1]-ave_diff[1])^2+(diff4[1]-ave_diff[1])^2+(diff5[1]-ave_diff[1])^2)/4  #between-imputation variance
  var.ave_diff <- (vars.bar+(1+1/5)*v.impute)  #variance for average diff
  logv <- var.ave_diff/ave_diff[1]^2  #variance of log of average proportion
  #confidence intervals for average proportion
  lower <- exp(log(ave_diff[1])-1.96*sqrt(logv))
  upper <- exp(log(ave_diff[1])+1.96*sqrt(logv))
  ave_diff <- rbind(ave_diff[1],lower,upper,2*(pnorm(log(ave_diff[1])/sqrt(logv)))) 
  rownames(ave_diff) <- c("mean","lower","upper","2-sided pvalue")
  return(ave_diff)
}

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod1.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>84)
nhis <- subset(nhis,analypop==1&age>=40&age<=84)  #impute variables for analysis population

persons <- data.frame(age=nhis$age,
                      female=nhis$female,
                      smkyears=nhis$smkyears,
                      qtyears=nhis$qtyears,
                      cpd=nhis$cpd,
                      race=as.numeric(as.character(nhis$race)),
                      emp=nhis$emp,
                      fam.lung.trend=nhis$fam.lung.trend,
                      bmi=nhis$bmi,
                      edu6=nhis$edu6)
predict <- lcrisk(persons,5)
predict2 <- data.frame(matrix(0,nrow(notanaly),7))             
colnames(predict2) <- colnames(predict)
nhis <- rbind(nhis,notanaly)
nhis$predict <- rbind(predict,predict2)
nhis$lcd5 <- nhis$predict[,3]/1000
nhis$LCRAT <- nhis$predict[,5]/1000
nhis$cxLCRAT <- 1.124*nhis$LCRAT
nhis$expected_falsepos <- 0
nhis$expected_falsepos[nhis$cxLCRAT>0] <- as.numeric(predict(polytmod,type="response",newdata=nhis[nhis$cxLCRAT>0,]) %*% c(0,1,2,3))


nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0
nhis$lcds5[is.na(nhis$lcds5)] <- 0
nhis$uspstf.eligible[is.na(nhis$uspstf.eligible)] <- 0
  
master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$lcd5
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$lcd5>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

#select USPSTF size population based on LYG model (40-80)
risk <- nhis$lyg
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lyg.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lyg.eligible <- ifelse(nhis$lyg>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

nhis$lc.uspstf <- nhis$uspstf.eligible*nhis$cxLCRAT
nhis$lcds.uspstf <- nhis$uspstf.eligible*nhis$lcds5
nhis$lyg.uspstf <- nhis$uspstf.eligible*nhis$lyg
nhis$fp.uspstf <- nhis$uspstf.eligible*nhis$expected_falsepos

nhis$lc.lcdrat <- nhis$lcdrat.eligible*nhis$cxLCRAT
nhis$lcds.lcdrat <- nhis$lcdrat.eligible*nhis$lcds5
nhis$lyg.lcdrat <- nhis$lcdrat.eligible*nhis$lyg
nhis$fp.lcdrat <- nhis$lcdrat.eligible*nhis$expected_falsepos

nhis$lc.lyg <- nhis$lyg.eligible*nhis$cxLCRAT
nhis$lcds.lyg <- nhis$lyg.eligible*nhis$lcds5
nhis$lyg.lyg <- nhis$lyg.eligible*nhis$lyg
nhis$fp.lyg <- nhis$lyg.eligible*nhis$expected_falsepos

nhis$diff1 <- nhis$lyg.lyg - nhis$lyg.lcdrat
nhis$diff2 <- nhis$fp.lcdrat - nhis$fp.lyg
nhis$diff3 <- nhis$lcds.lcdrat - nhis$lcds.lyg

nhis$age50_59 <- ifelse(nhis$age>=50 & nhis$age<=59,1,0)
nhis$age60_69 <- ifelse(nhis$age>=60 & nhis$age<=69,1,0)
nhis$age70plus <- ifelse(nhis$age>=70,1,0)
nhis$age50_59.lcdrat <- nhis$lcdrat.eligible*nhis$age50_59
nhis$age60_69.lcdrat <- nhis$lcdrat.eligible*nhis$age60_69
nhis$age70plus.lcdrat <- nhis$lcdrat.eligible*nhis$age70plus
nhis$age50_59.lyg <- nhis$lyg.eligible*nhis$age50_59
nhis$age60_69.lyg <- nhis$lyg.eligible*nhis$age60_69
nhis$age70plus.lyg <- nhis$lyg.eligible*nhis$age70plus
nhis$diff4 <- nhis$age50_59.lcdrat - nhis$age50_59.lyg
nhis$diff5 <- nhis$age60_69.lcdrat - nhis$age60_69.lyg
nhis$diff6 <- nhis$age70plus.lcdrat - nhis$age70plus.lyg

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)

uspstf.1 <- svytotal(~uspstf.eligible,master)[[1]]
se.uspstf.1 <- SE(svytotal(~uspstf.eligible,master))
uspstf.2 <- svyratio(~lc.uspstf,~cxLCRAT,master)[[1]]
se.uspstf.2 <- SE(svyratio(~lc.uspstf,~cxLCRAT,master))
uspstf.3 <- svyratio(~lcds.uspstf,~lcds5,master)[[1]]
se.uspstf.3 <- SE(svyratio(~lcds.uspstf,~lcds5,master))
uspstf.4 <- svyratio(~lyg.uspstf,~lyg,master)[[1]]
se.uspstf.4 <- SE(svyratio(~lyg.uspstf,~lyg,master))
uspstf.5 <- svytotal(~lc.uspstf,master)[[1]]
se.uspstf.5 <- SE(svytotal(~lc.uspstf,master))
uspstf.6 <- svytotal(~lcds.uspstf,master)[[1]]
se.uspstf.6 <- SE(svytotal(~lcds.uspstf,master))
uspstf.7 <- svytotal(~lyg.uspstf,master)[[1]]
se.uspstf.7 <- SE(svytotal(~lyg.uspstf,master))
uspstf.8 <- svyratio(~uspstf.eligible,~lcds.uspstf,uspstf)[[1]]
se.uspstf.8 <- SE(svyratio(~uspstf.eligible,~lcds.uspstf,uspstf))
uspstf.9 <- svyratio(~uspstf.eligible,~lyg.uspstf,uspstf)[[1]]
se.uspstf.9 <- SE(svyratio(~uspstf.eligible,~lyg.uspstf,uspstf))
uspstf.10 <- svyratio(~lyg.uspstf,~lc.uspstf,uspstf)[[1]]
se.uspstf.10 <- SE(svyratio(~lyg.uspstf,~lc.uspstf,uspstf))
uspstf.11 <- svyratio(~lyg.uspstf,~lcds.uspstf,uspstf)[[1]]
se.uspstf.11 <- SE(svyratio(~lyg.uspstf,~lcds.uspstf,uspstf))
uspstf.12 <- svyratio(~fp.uspstf,~lcds.uspstf,uspstf)[[1]]
se.uspstf.12 <- SE(svyratio(~fp.uspstf,~lcds.uspstf,uspstf))
uspstf.13 <- svyratio(~fp.uspstf,~lyg.uspstf,uspstf)[[1]]
se.uspstf.13 <- SE(svyratio(~fp.uspstf,~lyg.uspstf,uspstf))
uspstf.14 <- svytotal(~fp.uspstf,uspstf)[1]
se.uspstf.14 <- SE(svytotal(~fp.uspstf,uspstf))

lcdrat.1 <- svytotal(~lcdrat.eligible,master)[[1]]
se.lcdrat.1 <- SE(svytotal(~lcdrat.eligible,master))
lcdrat.2 <- svyratio(~lc.lcdrat,~cxLCRAT,master)[[1]]
se.lcdrat.2 <- SE(svyratio(~lc.lcdrat,~cxLCRAT,master))
lcdrat.3 <- svyratio(~lcds.lcdrat,~lcds5,master)[[1]]
se.lcdrat.3 <- SE(svyratio(~lcds.lcdrat,~lcds5,master))
lcdrat.4 <- svyratio(~lyg.lcdrat,~lyg,master)[[1]]
se.lcdrat.4 <- SE(svyratio(~lyg.lcdrat,~lyg,master))
lcdrat.5 <- svytotal(~lc.lcdrat,master)[[1]]
se.lcdrat.5 <- SE(svytotal(~lc.lcdrat,master))
lcdrat.6 <- svytotal(~lcds.lcdrat,master)[[1]]
se.lcdrat.6 <- SE(svytotal(~lcds.lcdrat,master))
lcdrat.7 <- svytotal(~lyg.lcdrat,master)[[1]]
se.lcdrat.7 <- SE(svytotal(~lyg.lcdrat,master))
lcdrat.8 <- svyratio(~lcdrat.eligible,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.8 <- SE(svyratio(~lcdrat.eligible,~lcds.lcdrat,lcdrat))
lcdrat.9 <- svyratio(~lcdrat.eligible,~lyg.lcdrat,lcdrat)[[1]]
se.lcdrat.9 <- SE(svyratio(~lcdrat.eligible,~lyg.lcdrat,lcdrat))
lcdrat.10 <- svyratio(~lyg.lcdrat,~lc.lcdrat,lcdrat)[[1]]
se.lcdrat.10 <- SE(svyratio(~lyg.lcdrat,~lc.lcdrat,lcdrat))
lcdrat.11 <- svyratio(~lyg.lcdrat,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.11 <- SE(svyratio(~lyg.lcdrat,~lcds.lcdrat,lcdrat))
lcdrat.12 <- svyratio(~fp.lcdrat,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.12 <- SE(svyratio(~fp.lcdrat,~lcds.lcdrat,lcdrat))
lcdrat.13 <- svyratio(~fp.lcdrat,~lyg.lcdrat,lcdrat)[[1]]
se.lcdrat.13 <- SE(svyratio(~fp.lcdrat,~lyg.lcdrat,lcdrat))
lcdrat.14 <- svytotal(~fp.lcdrat,lcdrat)[1]
se.lcdrat.14 <- SE(svytotal(~fp.lcdrat,lcdrat))
lcdrat.15 <- svyratio()

lyg.1 <- svytotal(~lyg.eligible,master)[[1]]
se.lyg.1 <- SE(svytotal(~lyg.eligible,master))
lyg.2 <- svyratio(~lc.lyg,~cxLCRAT,master)[[1]]
se.lyg.2 <- SE(svyratio(~lc.lyg,~cxLCRAT,master))
lyg.3 <- svyratio(~lcds.lyg,~lcds5,master)[[1]]
se.lyg.3 <- SE(svyratio(~lcds.lyg,~lcds5,master))
lyg.4 <- svyratio(~lyg.lyg,~lyg,master)[[1]]
se.lyg.4 <- SE(svyratio(~lyg.lyg,~lyg,master))
lyg.5 <- svytotal(~lc.lyg,master)[[1]]
se.lyg.5 <- SE(svytotal(~lc.lyg,master))
lyg.6 <- svytotal(~lcds.lyg,master)[[1]]
se.lyg.6 <- SE(svytotal(~lcds.lyg,master))
lyg.7 <- svytotal(~lyg.lyg,master)[[1]]
se.lyg.7 <- SE(svytotal(~lyg.lyg,master))
lyg.8 <- svyratio(~lyg.eligible,~lcds.lyg,lyg)[[1]]
se.lyg.8 <- SE(svyratio(~lyg.eligible,~lcds.lyg,lyg))
lyg.9 <- svyratio(~lyg.eligible,~lyg.lyg,lyg)[[1]]
se.lyg.9 <- SE(svyratio(~lyg.eligible,~lyg.lyg,lyg))
lyg.10 <- svyratio(~lyg.lyg,~lc.lyg,lyg)[[1]]
se.lyg.10 <- SE(svyratio(~lyg.lyg,~lc.lyg,lyg))
lyg.11 <- svyratio(~lyg.lyg,~lcds.lyg,lyg)[[1]]
se.lyg.11 <- SE(svyratio(~lyg.lyg,~lcds.lyg,lyg))
lyg.12 <- svyratio(~fp.lyg,~lcds.lyg,lyg)[[1]]
se.lyg.12 <- SE(svyratio(~fp.lyg,~lcds.lyg,lyg))
lyg.13 <- svyratio(~fp.lyg,~lyg.lyg,lyg)[[1]]
se.lyg.13 <- SE(svyratio(~fp.lyg,~lyg.lyg,lyg))
lyg.14 <- svytotal(~fp.lyg,lyg)[1]
se.lyg.14 <- SE(svytotal(~fp.lyg,lyg))

diff.1 <- svytotal(~diff1,master)[1]
se.diff.1 <- SE(svytotal(~diff1,master))
diff.2 <- svytotal(~diff2,master)[1]
se.diff.2 <- SE(svytotal(~diff2,master))
diff.3 <- svyratio(~diff3,~lcds5,master)[[1]]
se.diff.3 <- SE(svyratio(~diff3,~lcds5,master))

diff.4 <- svyratio(~diff4,~age50_59,master)[[1]]
se.diff.4 <- SE(svyratio(~diff4,~age50_59,master))
diff.5 <- svyratio(~diff5,~age60_69,master)[[1]]
se.diff.5 <- SE(svyratio(~diff5,~age60_69,master))
diff.6 <- svyratio(~diff6,~age70plus,master)[[1]]
se.diff.6 <- SE(svyratio(~diff6,~age70plus,master))

est.1 <- rbind(c(uspstf.1,uspstf.5,uspstf.6,uspstf.7,uspstf.8,uspstf.9,uspstf.10,uspstf.11,uspstf.12,uspstf.13,uspstf.14,
                 lcdrat.1,lcdrat.5,lcdrat.6,lcdrat.7,lcdrat.8,lcdrat.9,lcdrat.10,lcdrat.11,lcdrat.12,lcdrat.13,lcdrat.14,
                 lyg.1,lyg.5,lyg.6,lyg.7,lyg.8,lyg.9,lyg.10,lyg.11,lyg.12,lyg.13,lyg.14),
               c(se.uspstf.1,se.uspstf.5,se.uspstf.6,se.uspstf.7,se.uspstf.8,se.uspstf.9,se.uspstf.10,se.uspstf.11,se.uspstf.12,se.uspstf.13,se.uspstf.14,
                 se.lcdrat.1,se.lcdrat.5,se.lcdrat.6,se.lcdrat.7,se.lcdrat.8,se.lcdrat.9,se.lcdrat.10,se.lcdrat.11,se.lcdrat.12,se.lcdrat.13,se.lcdrat.14,
                 se.lyg.1,se.lyg.5,se.lyg.6,se.lyg.7,se.lyg.8,se.lyg.9,se.lyg.10,se.lyg.11,se.lyg.12,se.lyg.13,se.lyg.14)^2)
rat.1 <- rbind(c(uspstf.2,uspstf.3,uspstf.4,
                 lcdrat.2,lcdrat.3,lcdrat.4,
                 lyg.2,lyg.3,lyg.4),
               c(se.uspstf.2,se.uspstf.3,se.uspstf.4,
                 se.lcdrat.2,se.lcdrat.3,se.lcdrat.4,
                 se.lyg.2,se.lyg.3,se.lyg.4)^2)
diff1 <- rbind(c(diff.1,diff.2,diff.3,diff.4,diff.5,diff.6),
                c(se.diff.1,se.diff.2,se.diff.3,se.diff.4,se.diff.5,se.diff.6)^2)

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod2.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>84)
nhis <- subset(nhis,analypop==1&age>=40&age<=84)  #impute variables for analysis population

persons <- data.frame(age=nhis$age,
                      female=nhis$female,
                      smkyears=nhis$smkyears,
                      qtyears=nhis$qtyears,
                      cpd=nhis$cpd,
                      race=as.numeric(as.character(nhis$race)),
                      emp=nhis$emp,
                      fam.lung.trend=nhis$fam.lung.trend,
                      bmi=nhis$bmi,
                      edu6=nhis$edu6)
predict <- lcrisk(persons,5)
predict2 <- data.frame(matrix(0,nrow(notanaly),7))             
colnames(predict2) <- colnames(predict)
nhis <- rbind(nhis,notanaly)
nhis$predict <- rbind(predict,predict2)
nhis$lcd5 <- nhis$predict[,3]/1000
nhis$LCRAT <- nhis$predict[,5]/1000
nhis$cxLCRAT <- 1.124*nhis$LCRAT
nhis$expected_falsepos <- 0
nhis$expected_falsepos[nhis$cxLCRAT>0] <- as.numeric(predict(polytmod,type="response",newdata=nhis[nhis$cxLCRAT>0,]) %*% c(0,1,2,3))


nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0
nhis$lcds5[is.na(nhis$lcds5)] <- 0
nhis$uspstf.eligible[is.na(nhis$uspstf.eligible)] <- 0

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$lcd5
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$lcd5>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

#select USPSTF size population based on LYG model (40-80)
risk <- nhis$lyg
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lyg.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lyg.eligible <- ifelse(nhis$lyg>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

nhis$lc.uspstf <- nhis$uspstf.eligible*nhis$cxLCRAT
nhis$lcds.uspstf <- nhis$uspstf.eligible*nhis$lcds5
nhis$lyg.uspstf <- nhis$uspstf.eligible*nhis$lyg
nhis$fp.uspstf <- nhis$uspstf.eligible*nhis$expected_falsepos

nhis$lc.lcdrat <- nhis$lcdrat.eligible*nhis$cxLCRAT
nhis$lcds.lcdrat <- nhis$lcdrat.eligible*nhis$lcds5
nhis$lyg.lcdrat <- nhis$lcdrat.eligible*nhis$lyg
nhis$fp.lcdrat <- nhis$lcdrat.eligible*nhis$expected_falsepos

nhis$lc.lyg <- nhis$lyg.eligible*nhis$cxLCRAT
nhis$lcds.lyg <- nhis$lyg.eligible*nhis$lcds5
nhis$lyg.lyg <- nhis$lyg.eligible*nhis$lyg
nhis$fp.lyg <- nhis$lyg.eligible*nhis$expected_falsepos

nhis$diff1 <- nhis$lyg.lyg - nhis$lyg.lcdrat
nhis$diff2 <- nhis$fp.lcdrat - nhis$fp.lyg
nhis$diff3 <- nhis$lcds.lcdrat - nhis$lcds.lyg

nhis$age50_59 <- ifelse(nhis$age>=50 & nhis$age<=59,1,0)
nhis$age60_69 <- ifelse(nhis$age>=60 & nhis$age<=69,1,0)
nhis$age70plus <- ifelse(nhis$age>=70,1,0)
nhis$age50_59.lcdrat <- nhis$lcdrat.eligible*nhis$age50_59
nhis$age60_69.lcdrat <- nhis$lcdrat.eligible*nhis$age60_69
nhis$age70plus.lcdrat <- nhis$lcdrat.eligible*nhis$age70plus
nhis$age50_59.lyg <- nhis$lyg.eligible*nhis$age50_59
nhis$age60_69.lyg <- nhis$lyg.eligible*nhis$age60_69
nhis$age70plus.lyg <- nhis$lyg.eligible*nhis$age70plus
nhis$diff4 <- nhis$age50_59.lcdrat - nhis$age50_59.lyg
nhis$diff5 <- nhis$age60_69.lcdrat - nhis$age60_69.lyg
nhis$diff6 <- nhis$age70plus.lcdrat - nhis$age70plus.lyg

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)

uspstf.1 <- svytotal(~uspstf.eligible,master)[[1]]
se.uspstf.1 <- SE(svytotal(~uspstf.eligible,master))
uspstf.2 <- svyratio(~lc.uspstf,~cxLCRAT,master)[[1]]
se.uspstf.2 <- SE(svyratio(~lc.uspstf,~cxLCRAT,master))
uspstf.3 <- svyratio(~lcds.uspstf,~lcds5,master)[[1]]
se.uspstf.3 <- SE(svyratio(~lcds.uspstf,~lcds5,master))
uspstf.4 <- svyratio(~lyg.uspstf,~lyg,master)[[1]]
se.uspstf.4 <- SE(svyratio(~lyg.uspstf,~lyg,master))
uspstf.5 <- svytotal(~lc.uspstf,master)[[1]]
se.uspstf.5 <- SE(svytotal(~lc.uspstf,master))
uspstf.6 <- svytotal(~lcds.uspstf,master)[[1]]
se.uspstf.6 <- SE(svytotal(~lcds.uspstf,master))
uspstf.7 <- svytotal(~lyg.uspstf,master)[[1]]
se.uspstf.7 <- SE(svytotal(~lyg.uspstf,master))
uspstf.8 <- svyratio(~uspstf.eligible,~lcds.uspstf,uspstf)[[1]]
se.uspstf.8 <- SE(svyratio(~uspstf.eligible,~lcds.uspstf,uspstf))
uspstf.9 <- svyratio(~uspstf.eligible,~lyg.uspstf,uspstf)[[1]]
se.uspstf.9 <- SE(svyratio(~uspstf.eligible,~lyg.uspstf,uspstf))
uspstf.10 <- svyratio(~lyg.uspstf,~lc.uspstf,uspstf)[[1]]
se.uspstf.10 <- SE(svyratio(~lyg.uspstf,~lc.uspstf,uspstf))
uspstf.11 <- svyratio(~lyg.uspstf,~lcds.uspstf,uspstf)[[1]]
se.uspstf.11 <- SE(svyratio(~lyg.uspstf,~lcds.uspstf,uspstf))
uspstf.12 <- svyratio(~fp.uspstf,~lcds.uspstf,uspstf)[[1]]
se.uspstf.12 <- SE(svyratio(~fp.uspstf,~lcds.uspstf,uspstf))
uspstf.13 <- svyratio(~fp.uspstf,~lyg.uspstf,uspstf)[[1]]
se.uspstf.13 <- SE(svyratio(~fp.uspstf,~lyg.uspstf,uspstf))
uspstf.14 <- svytotal(~fp.uspstf,uspstf)[1]
se.uspstf.14 <- SE(svytotal(~fp.uspstf,uspstf))

lcdrat.1 <- svytotal(~lcdrat.eligible,master)[[1]]
se.lcdrat.1 <- SE(svytotal(~lcdrat.eligible,master))
lcdrat.2 <- svyratio(~lc.lcdrat,~cxLCRAT,master)[[1]]
se.lcdrat.2 <- SE(svyratio(~lc.lcdrat,~cxLCRAT,master))
lcdrat.3 <- svyratio(~lcds.lcdrat,~lcds5,master)[[1]]
se.lcdrat.3 <- SE(svyratio(~lcds.lcdrat,~lcds5,master))
lcdrat.4 <- svyratio(~lyg.lcdrat,~lyg,master)[[1]]
se.lcdrat.4 <- SE(svyratio(~lyg.lcdrat,~lyg,master))
lcdrat.5 <- svytotal(~lc.lcdrat,master)[[1]]
se.lcdrat.5 <- SE(svytotal(~lc.lcdrat,master))
lcdrat.6 <- svytotal(~lcds.lcdrat,master)[[1]]
se.lcdrat.6 <- SE(svytotal(~lcds.lcdrat,master))
lcdrat.7 <- svytotal(~lyg.lcdrat,master)[[1]]
se.lcdrat.7 <- SE(svytotal(~lyg.lcdrat,master))
lcdrat.8 <- svyratio(~lcdrat.eligible,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.8 <- SE(svyratio(~lcdrat.eligible,~lcds.lcdrat,lcdrat))
lcdrat.9 <- svyratio(~lcdrat.eligible,~lyg.lcdrat,lcdrat)[[1]]
se.lcdrat.9 <- SE(svyratio(~lcdrat.eligible,~lyg.lcdrat,lcdrat))
lcdrat.10 <- svyratio(~lyg.lcdrat,~lc.lcdrat,lcdrat)[[1]]
se.lcdrat.10 <- SE(svyratio(~lyg.lcdrat,~lc.lcdrat,lcdrat))
lcdrat.11 <- svyratio(~lyg.lcdrat,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.11 <- SE(svyratio(~lyg.lcdrat,~lcds.lcdrat,lcdrat))
lcdrat.12 <- svyratio(~fp.lcdrat,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.12 <- SE(svyratio(~fp.lcdrat,~lcds.lcdrat,lcdrat))
lcdrat.13 <- svyratio(~fp.lcdrat,~lyg.lcdrat,lcdrat)[[1]]
se.lcdrat.13 <- SE(svyratio(~fp.lcdrat,~lyg.lcdrat,lcdrat))
lcdrat.14 <- svytotal(~fp.lcdrat,lcdrat)[1]
se.lcdrat.14 <- SE(svytotal(~fp.lcdrat,lcdrat))

lyg.1 <- svytotal(~lyg.eligible,master)[[1]]
se.lyg.1 <- SE(svytotal(~lyg.eligible,master))
lyg.2 <- svyratio(~lc.lyg,~cxLCRAT,master)[[1]]
se.lyg.2 <- SE(svyratio(~lc.lyg,~cxLCRAT,master))
lyg.3 <- svyratio(~lcds.lyg,~lcds5,master)[[1]]
se.lyg.3 <- SE(svyratio(~lcds.lyg,~lcds5,master))
lyg.4 <- svyratio(~lyg.lyg,~lyg,master)[[1]]
se.lyg.4 <- SE(svyratio(~lyg.lyg,~lyg,master))
lyg.5 <- svytotal(~lc.lyg,master)[[1]]
se.lyg.5 <- SE(svytotal(~lc.lyg,master))
lyg.6 <- svytotal(~lcds.lyg,master)[[1]]
se.lyg.6 <- SE(svytotal(~lcds.lyg,master))
lyg.7 <- svytotal(~lyg.lyg,master)[[1]]
se.lyg.7 <- SE(svytotal(~lyg.lyg,master))
lyg.8 <- svyratio(~lyg.eligible,~lcds.lyg,lyg)[[1]]
se.lyg.8 <- SE(svyratio(~lyg.eligible,~lcds.lyg,lyg))
lyg.9 <- svyratio(~lyg.eligible,~lyg.lyg,lyg)[[1]]
se.lyg.9 <- SE(svyratio(~lyg.eligible,~lyg.lyg,lyg))
lyg.10 <- svyratio(~lyg.lyg,~lc.lyg,lyg)[[1]]
se.lyg.10 <- SE(svyratio(~lyg.lyg,~lc.lyg,lyg))
lyg.11 <- svyratio(~lyg.lyg,~lcds.lyg,lyg)[[1]]
se.lyg.11 <- SE(svyratio(~lyg.lyg,~lcds.lyg,lyg))
lyg.12 <- svyratio(~fp.lyg,~lcds.lyg,lyg)[[1]]
se.lyg.12 <- SE(svyratio(~fp.lyg,~lcds.lyg,lyg))
lyg.13 <- svyratio(~fp.lyg,~lyg.lyg,lyg)[[1]]
se.lyg.13 <- SE(svyratio(~fp.lyg,~lyg.lyg,lyg))
lyg.14 <- svytotal(~fp.lyg,lyg)[1]
se.lyg.14 <- SE(svytotal(~fp.lyg,lyg))

diff.1 <- svytotal(~diff1,master)[1]
se.diff.1 <- SE(svytotal(~diff1,master))
diff.2 <- svytotal(~diff2,master)[1]
se.diff.2 <- SE(svytotal(~diff2,master))
diff.3 <- svyratio(~diff3,~lcds5,master)[[1]]
se.diff.3 <- SE(svyratio(~diff3,~lcds5,master))

diff.4 <- svyratio(~diff4,~age50_59,master)[[1]]
se.diff.4 <- SE(svyratio(~diff4,~age50_59,master))
diff.5 <- svyratio(~diff5,~age60_69,master)[[1]]
se.diff.5 <- SE(svyratio(~diff5,~age60_69,master))
diff.6 <- svyratio(~diff6,~age70plus,master)[[1]]
se.diff.6 <- SE(svyratio(~diff6,~age70plus,master))

est.2 <- rbind(c(uspstf.1,uspstf.5,uspstf.6,uspstf.7,uspstf.8,uspstf.9,uspstf.10,uspstf.11,uspstf.12,uspstf.13,uspstf.14,
                 lcdrat.1,lcdrat.5,lcdrat.6,lcdrat.7,lcdrat.8,lcdrat.9,lcdrat.10,lcdrat.11,lcdrat.12,lcdrat.13,lcdrat.14,
                 lyg.1,lyg.5,lyg.6,lyg.7,lyg.8,lyg.9,lyg.10,lyg.11,lyg.12,lyg.13,lyg.14),
               c(se.uspstf.1,se.uspstf.5,se.uspstf.6,se.uspstf.7,se.uspstf.8,se.uspstf.9,se.uspstf.10,se.uspstf.11,se.uspstf.12,se.uspstf.13,se.uspstf.14,
                 se.lcdrat.1,se.lcdrat.5,se.lcdrat.6,se.lcdrat.7,se.lcdrat.8,se.lcdrat.9,se.lcdrat.10,se.lcdrat.11,se.lcdrat.12,se.lcdrat.13,se.lcdrat.14,
                 se.lyg.1,se.lyg.5,se.lyg.6,se.lyg.7,se.lyg.8,se.lyg.9,se.lyg.10,se.lyg.11,se.lyg.12,se.lyg.13,se.lyg.14)^2)
rat.2 <- rbind(c(uspstf.2,uspstf.3,uspstf.4,
                 lcdrat.2,lcdrat.3,lcdrat.4,
                 lyg.2,lyg.3,lyg.4),
               c(se.uspstf.2,se.uspstf.3,se.uspstf.4,
                 se.lcdrat.2,se.lcdrat.3,se.lcdrat.4,
                 se.lyg.2,se.lyg.3,se.lyg.4)^2)
diff2 <- rbind(c(diff.1,diff.2,diff.3,diff.4,diff.5,diff.6),
                c(se.diff.1,se.diff.2,se.diff.3,se.diff.4,se.diff.5,se.diff.6)^2)


load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod3.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>84)
nhis <- subset(nhis,analypop==1&age>=40&age<=84)  #impute variables for analysis population

persons <- data.frame(age=nhis$age,
                      female=nhis$female,
                      smkyears=nhis$smkyears,
                      qtyears=nhis$qtyears,
                      cpd=nhis$cpd,
                      race=as.numeric(as.character(nhis$race)),
                      emp=nhis$emp,
                      fam.lung.trend=nhis$fam.lung.trend,
                      bmi=nhis$bmi,
                      edu6=nhis$edu6)
predict <- lcrisk(persons,5)
predict2 <- data.frame(matrix(0,nrow(notanaly),7))             
colnames(predict2) <- colnames(predict)
nhis <- rbind(nhis,notanaly)
nhis$predict <- rbind(predict,predict2)
nhis$lcd5 <- nhis$predict[,3]/1000
nhis$LCRAT <- nhis$predict[,5]/1000
nhis$cxLCRAT <- 1.124*nhis$LCRAT
nhis$expected_falsepos <- 0
nhis$expected_falsepos[nhis$cxLCRAT>0] <- as.numeric(predict(polytmod,type="response",newdata=nhis[nhis$cxLCRAT>0,]) %*% c(0,1,2,3))


nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0
nhis$lcds5[is.na(nhis$lcds5)] <- 0
nhis$uspstf.eligible[is.na(nhis$uspstf.eligible)] <- 0

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$lcd5
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$lcd5>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

#select USPSTF size population based on LYG model (40-80)
risk <- nhis$lyg
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lyg.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lyg.eligible <- ifelse(nhis$lyg>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

nhis$lc.uspstf <- nhis$uspstf.eligible*nhis$cxLCRAT
nhis$lcds.uspstf <- nhis$uspstf.eligible*nhis$lcds5
nhis$lyg.uspstf <- nhis$uspstf.eligible*nhis$lyg
nhis$fp.uspstf <- nhis$uspstf.eligible*nhis$expected_falsepos

nhis$lc.lcdrat <- nhis$lcdrat.eligible*nhis$cxLCRAT
nhis$lcds.lcdrat <- nhis$lcdrat.eligible*nhis$lcds5
nhis$lyg.lcdrat <- nhis$lcdrat.eligible*nhis$lyg
nhis$fp.lcdrat <- nhis$lcdrat.eligible*nhis$expected_falsepos

nhis$lc.lyg <- nhis$lyg.eligible*nhis$cxLCRAT
nhis$lcds.lyg <- nhis$lyg.eligible*nhis$lcds5
nhis$lyg.lyg <- nhis$lyg.eligible*nhis$lyg
nhis$fp.lyg <- nhis$lyg.eligible*nhis$expected_falsepos

nhis$diff1 <- nhis$lyg.lyg - nhis$lyg.lcdrat
nhis$diff2 <- nhis$fp.lcdrat - nhis$fp.lyg
nhis$diff3 <- nhis$lcds.lcdrat - nhis$lcds.lyg

nhis$age50_59 <- ifelse(nhis$age>=50 & nhis$age<=59,1,0)
nhis$age60_69 <- ifelse(nhis$age>=60 & nhis$age<=69,1,0)
nhis$age70plus <- ifelse(nhis$age>=70,1,0)
nhis$age50_59.lcdrat <- nhis$lcdrat.eligible*nhis$age50_59
nhis$age60_69.lcdrat <- nhis$lcdrat.eligible*nhis$age60_69
nhis$age70plus.lcdrat <- nhis$lcdrat.eligible*nhis$age70plus
nhis$age50_59.lyg <- nhis$lyg.eligible*nhis$age50_59
nhis$age60_69.lyg <- nhis$lyg.eligible*nhis$age60_69
nhis$age70plus.lyg <- nhis$lyg.eligible*nhis$age70plus
nhis$diff4 <- nhis$age50_59.lcdrat - nhis$age50_59.lyg
nhis$diff5 <- nhis$age60_69.lcdrat - nhis$age60_69.lyg
nhis$diff6 <- nhis$age70plus.lcdrat - nhis$age70plus.lyg

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)

uspstf.1 <- svytotal(~uspstf.eligible,master)[[1]]
se.uspstf.1 <- SE(svytotal(~uspstf.eligible,master))
uspstf.2 <- svyratio(~lc.uspstf,~cxLCRAT,master)[[1]]
se.uspstf.2 <- SE(svyratio(~lc.uspstf,~cxLCRAT,master))
uspstf.3 <- svyratio(~lcds.uspstf,~lcds5,master)[[1]]
se.uspstf.3 <- SE(svyratio(~lcds.uspstf,~lcds5,master))
uspstf.4 <- svyratio(~lyg.uspstf,~lyg,master)[[1]]
se.uspstf.4 <- SE(svyratio(~lyg.uspstf,~lyg,master))
uspstf.5 <- svytotal(~lc.uspstf,master)[[1]]
se.uspstf.5 <- SE(svytotal(~lc.uspstf,master))
uspstf.6 <- svytotal(~lcds.uspstf,master)[[1]]
se.uspstf.6 <- SE(svytotal(~lcds.uspstf,master))
uspstf.7 <- svytotal(~lyg.uspstf,master)[[1]]
se.uspstf.7 <- SE(svytotal(~lyg.uspstf,master))
uspstf.8 <- svyratio(~uspstf.eligible,~lcds.uspstf,uspstf)[[1]]
se.uspstf.8 <- SE(svyratio(~uspstf.eligible,~lcds.uspstf,uspstf))
uspstf.9 <- svyratio(~uspstf.eligible,~lyg.uspstf,uspstf)[[1]]
se.uspstf.9 <- SE(svyratio(~uspstf.eligible,~lyg.uspstf,uspstf))
uspstf.10 <- svyratio(~lyg.uspstf,~lc.uspstf,uspstf)[[1]]
se.uspstf.10 <- SE(svyratio(~lyg.uspstf,~lc.uspstf,uspstf))
uspstf.11 <- svyratio(~lyg.uspstf,~lcds.uspstf,uspstf)[[1]]
se.uspstf.11 <- SE(svyratio(~lyg.uspstf,~lcds.uspstf,uspstf))
uspstf.12 <- svyratio(~fp.uspstf,~lcds.uspstf,uspstf)[[1]]
se.uspstf.12 <- SE(svyratio(~fp.uspstf,~lcds.uspstf,uspstf))
uspstf.13 <- svyratio(~fp.uspstf,~lyg.uspstf,uspstf)[[1]]
se.uspstf.13 <- SE(svyratio(~fp.uspstf,~lyg.uspstf,uspstf))
uspstf.14 <- svytotal(~fp.uspstf,uspstf)[1]
se.uspstf.14 <- SE(svytotal(~fp.uspstf,uspstf))

lcdrat.1 <- svytotal(~lcdrat.eligible,master)[[1]]
se.lcdrat.1 <- SE(svytotal(~lcdrat.eligible,master))
lcdrat.2 <- svyratio(~lc.lcdrat,~cxLCRAT,master)[[1]]
se.lcdrat.2 <- SE(svyratio(~lc.lcdrat,~cxLCRAT,master))
lcdrat.3 <- svyratio(~lcds.lcdrat,~lcds5,master)[[1]]
se.lcdrat.3 <- SE(svyratio(~lcds.lcdrat,~lcds5,master))
lcdrat.4 <- svyratio(~lyg.lcdrat,~lyg,master)[[1]]
se.lcdrat.4 <- SE(svyratio(~lyg.lcdrat,~lyg,master))
lcdrat.5 <- svytotal(~lc.lcdrat,master)[[1]]
se.lcdrat.5 <- SE(svytotal(~lc.lcdrat,master))
lcdrat.6 <- svytotal(~lcds.lcdrat,master)[[1]]
se.lcdrat.6 <- SE(svytotal(~lcds.lcdrat,master))
lcdrat.7 <- svytotal(~lyg.lcdrat,master)[[1]]
se.lcdrat.7 <- SE(svytotal(~lyg.lcdrat,master))
lcdrat.8 <- svyratio(~lcdrat.eligible,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.8 <- SE(svyratio(~lcdrat.eligible,~lcds.lcdrat,lcdrat))
lcdrat.9 <- svyratio(~lcdrat.eligible,~lyg.lcdrat,lcdrat)[[1]]
se.lcdrat.9 <- SE(svyratio(~lcdrat.eligible,~lyg.lcdrat,lcdrat))
lcdrat.10 <- svyratio(~lyg.lcdrat,~lc.lcdrat,lcdrat)[[1]]
se.lcdrat.10 <- SE(svyratio(~lyg.lcdrat,~lc.lcdrat,lcdrat))
lcdrat.11 <- svyratio(~lyg.lcdrat,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.11 <- SE(svyratio(~lyg.lcdrat,~lcds.lcdrat,lcdrat))
lcdrat.12 <- svyratio(~fp.lcdrat,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.12 <- SE(svyratio(~fp.lcdrat,~lcds.lcdrat,lcdrat))
lcdrat.13 <- svyratio(~fp.lcdrat,~lyg.lcdrat,lcdrat)[[1]]
se.lcdrat.13 <- SE(svyratio(~fp.lcdrat,~lyg.lcdrat,lcdrat))
lcdrat.14 <- svytotal(~fp.lcdrat,lcdrat)[1]
se.lcdrat.14 <- SE(svytotal(~fp.lcdrat,lcdrat))

lyg.1 <- svytotal(~lyg.eligible,master)[[1]]
se.lyg.1 <- SE(svytotal(~lyg.eligible,master))
lyg.2 <- svyratio(~lc.lyg,~cxLCRAT,master)[[1]]
se.lyg.2 <- SE(svyratio(~lc.lyg,~cxLCRAT,master))
lyg.3 <- svyratio(~lcds.lyg,~lcds5,master)[[1]]
se.lyg.3 <- SE(svyratio(~lcds.lyg,~lcds5,master))
lyg.4 <- svyratio(~lyg.lyg,~lyg,master)[[1]]
se.lyg.4 <- SE(svyratio(~lyg.lyg,~lyg,master))
lyg.5 <- svytotal(~lc.lyg,master)[[1]]
se.lyg.5 <- SE(svytotal(~lc.lyg,master))
lyg.6 <- svytotal(~lcds.lyg,master)[[1]]
se.lyg.6 <- SE(svytotal(~lcds.lyg,master))
lyg.7 <- svytotal(~lyg.lyg,master)[[1]]
se.lyg.7 <- SE(svytotal(~lyg.lyg,master))
lyg.8 <- svyratio(~lyg.eligible,~lcds.lyg,lyg)[[1]]
se.lyg.8 <- SE(svyratio(~lyg.eligible,~lcds.lyg,lyg))
lyg.9 <- svyratio(~lyg.eligible,~lyg.lyg,lyg)[[1]]
se.lyg.9 <- SE(svyratio(~lyg.eligible,~lyg.lyg,lyg))
lyg.10 <- svyratio(~lyg.lyg,~lc.lyg,lyg)[[1]]
se.lyg.10 <- SE(svyratio(~lyg.lyg,~lc.lyg,lyg))
lyg.11 <- svyratio(~lyg.lyg,~lcds.lyg,lyg)[[1]]
se.lyg.11 <- SE(svyratio(~lyg.lyg,~lcds.lyg,lyg))
lyg.12 <- svyratio(~fp.lyg,~lcds.lyg,lyg)[[1]]
se.lyg.12 <- SE(svyratio(~fp.lyg,~lcds.lyg,lyg))
lyg.13 <- svyratio(~fp.lyg,~lyg.lyg,lyg)[[1]]
se.lyg.13 <- SE(svyratio(~fp.lyg,~lyg.lyg,lyg))
lyg.14 <- svytotal(~fp.lyg,lyg)[1]
se.lyg.14 <- SE(svytotal(~fp.lyg,lyg))

diff.1 <- svytotal(~diff1,master)[1]
se.diff.1 <- SE(svytotal(~diff1,master))
diff.2 <- svytotal(~diff2,master)[1]
se.diff.2 <- SE(svytotal(~diff2,master))
diff.3 <- svyratio(~diff3,~lcds5,master)[[1]]
se.diff.3 <- SE(svyratio(~diff3,~lcds5,master))

diff.4 <- svyratio(~diff4,~age50_59,master)[[1]]
se.diff.4 <- SE(svyratio(~diff4,~age50_59,master))
diff.5 <- svyratio(~diff5,~age60_69,master)[[1]]
se.diff.5 <- SE(svyratio(~diff5,~age60_69,master))
diff.6 <- svyratio(~diff6,~age70plus,master)[[1]]
se.diff.6 <- SE(svyratio(~diff6,~age70plus,master))

est.3 <- rbind(c(uspstf.1,uspstf.5,uspstf.6,uspstf.7,uspstf.8,uspstf.9,uspstf.10,uspstf.11,uspstf.12,uspstf.13,uspstf.14,
                 lcdrat.1,lcdrat.5,lcdrat.6,lcdrat.7,lcdrat.8,lcdrat.9,lcdrat.10,lcdrat.11,lcdrat.12,lcdrat.13,lcdrat.14,
                 lyg.1,lyg.5,lyg.6,lyg.7,lyg.8,lyg.9,lyg.10,lyg.11,lyg.12,lyg.13,lyg.14),
               c(se.uspstf.1,se.uspstf.5,se.uspstf.6,se.uspstf.7,se.uspstf.8,se.uspstf.9,se.uspstf.10,se.uspstf.11,se.uspstf.12,se.uspstf.13,se.uspstf.14,
                 se.lcdrat.1,se.lcdrat.5,se.lcdrat.6,se.lcdrat.7,se.lcdrat.8,se.lcdrat.9,se.lcdrat.10,se.lcdrat.11,se.lcdrat.12,se.lcdrat.13,se.lcdrat.14,
                 se.lyg.1,se.lyg.5,se.lyg.6,se.lyg.7,se.lyg.8,se.lyg.9,se.lyg.10,se.lyg.11,se.lyg.12,se.lyg.13,se.lyg.14)^2)
rat.3 <- rbind(c(uspstf.2,uspstf.3,uspstf.4,
                 lcdrat.2,lcdrat.3,lcdrat.4,
                 lyg.2,lyg.3,lyg.4),
               c(se.uspstf.2,se.uspstf.3,se.uspstf.4,
                 se.lcdrat.2,se.lcdrat.3,se.lcdrat.4,
                 se.lyg.2,se.lyg.3,se.lyg.4)^2)
diff3 <- rbind(c(diff.1,diff.2,diff.3,diff.4,diff.5,diff.6),
                c(se.diff.1,se.diff.2,se.diff.3,se.diff.4,se.diff.5,se.diff.6)^2)



load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod4.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>84)
nhis <- subset(nhis,analypop==1&age>=40&age<=84)  #impute variables for analysis population

persons <- data.frame(age=nhis$age,
                      female=nhis$female,
                      smkyears=nhis$smkyears,
                      qtyears=nhis$qtyears,
                      cpd=nhis$cpd,
                      race=as.numeric(as.character(nhis$race)),
                      emp=nhis$emp,
                      fam.lung.trend=nhis$fam.lung.trend,
                      bmi=nhis$bmi,
                      edu6=nhis$edu6)
predict <- lcrisk(persons,5)
predict2 <- data.frame(matrix(0,nrow(notanaly),7))             
colnames(predict2) <- colnames(predict)
nhis <- rbind(nhis,notanaly)
nhis$predict <- rbind(predict,predict2)
nhis$lcd5 <- nhis$predict[,3]/1000
nhis$LCRAT <- nhis$predict[,5]/1000
nhis$cxLCRAT <- 1.124*nhis$LCRAT
nhis$expected_falsepos <- 0
nhis$expected_falsepos[nhis$cxLCRAT>0] <- as.numeric(predict(polytmod,type="response",newdata=nhis[nhis$cxLCRAT>0,]) %*% c(0,1,2,3))


nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0
nhis$lcds5[is.na(nhis$lcds5)] <- 0
nhis$uspstf.eligible[is.na(nhis$uspstf.eligible)] <- 0

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$lcd5
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$lcd5>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

#select USPSTF size population based on LYG model (40-80)
risk <- nhis$lyg
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lyg.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lyg.eligible <- ifelse(nhis$lyg>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

nhis$lc.uspstf <- nhis$uspstf.eligible*nhis$cxLCRAT
nhis$lcds.uspstf <- nhis$uspstf.eligible*nhis$lcds5
nhis$lyg.uspstf <- nhis$uspstf.eligible*nhis$lyg
nhis$fp.uspstf <- nhis$uspstf.eligible*nhis$expected_falsepos

nhis$lc.lcdrat <- nhis$lcdrat.eligible*nhis$cxLCRAT
nhis$lcds.lcdrat <- nhis$lcdrat.eligible*nhis$lcds5
nhis$lyg.lcdrat <- nhis$lcdrat.eligible*nhis$lyg
nhis$fp.lcdrat <- nhis$lcdrat.eligible*nhis$expected_falsepos

nhis$lc.lyg <- nhis$lyg.eligible*nhis$cxLCRAT
nhis$lcds.lyg <- nhis$lyg.eligible*nhis$lcds5
nhis$lyg.lyg <- nhis$lyg.eligible*nhis$lyg
nhis$fp.lyg <- nhis$lyg.eligible*nhis$expected_falsepos

nhis$diff1 <- nhis$lyg.lyg - nhis$lyg.lcdrat
nhis$diff2 <- nhis$fp.lcdrat - nhis$fp.lyg
nhis$diff3 <- nhis$lcds.lcdrat - nhis$lcds.lyg

nhis$age50_59 <- ifelse(nhis$age>=50 & nhis$age<=59,1,0)
nhis$age60_69 <- ifelse(nhis$age>=60 & nhis$age<=69,1,0)
nhis$age70plus <- ifelse(nhis$age>=70,1,0)
nhis$age50_59.lcdrat <- nhis$lcdrat.eligible*nhis$age50_59
nhis$age60_69.lcdrat <- nhis$lcdrat.eligible*nhis$age60_69
nhis$age70plus.lcdrat <- nhis$lcdrat.eligible*nhis$age70plus
nhis$age50_59.lyg <- nhis$lyg.eligible*nhis$age50_59
nhis$age60_69.lyg <- nhis$lyg.eligible*nhis$age60_69
nhis$age70plus.lyg <- nhis$lyg.eligible*nhis$age70plus
nhis$diff4 <- nhis$age50_59.lcdrat - nhis$age50_59.lyg
nhis$diff5 <- nhis$age60_69.lcdrat - nhis$age60_69.lyg
nhis$diff6 <- nhis$age70plus.lcdrat - nhis$age70plus.lyg

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)

uspstf.1 <- svytotal(~uspstf.eligible,master)[[1]]
se.uspstf.1 <- SE(svytotal(~uspstf.eligible,master))
uspstf.2 <- svyratio(~lc.uspstf,~cxLCRAT,master)[[1]]
se.uspstf.2 <- SE(svyratio(~lc.uspstf,~cxLCRAT,master))
uspstf.3 <- svyratio(~lcds.uspstf,~lcds5,master)[[1]]
se.uspstf.3 <- SE(svyratio(~lcds.uspstf,~lcds5,master))
uspstf.4 <- svyratio(~lyg.uspstf,~lyg,master)[[1]]
se.uspstf.4 <- SE(svyratio(~lyg.uspstf,~lyg,master))
uspstf.5 <- svytotal(~lc.uspstf,master)[[1]]
se.uspstf.5 <- SE(svytotal(~lc.uspstf,master))
uspstf.6 <- svytotal(~lcds.uspstf,master)[[1]]
se.uspstf.6 <- SE(svytotal(~lcds.uspstf,master))
uspstf.7 <- svytotal(~lyg.uspstf,master)[[1]]
se.uspstf.7 <- SE(svytotal(~lyg.uspstf,master))
uspstf.8 <- svyratio(~uspstf.eligible,~lcds.uspstf,uspstf)[[1]]
se.uspstf.8 <- SE(svyratio(~uspstf.eligible,~lcds.uspstf,uspstf))
uspstf.9 <- svyratio(~uspstf.eligible,~lyg.uspstf,uspstf)[[1]]
se.uspstf.9 <- SE(svyratio(~uspstf.eligible,~lyg.uspstf,uspstf))
uspstf.10 <- svyratio(~lyg.uspstf,~lc.uspstf,uspstf)[[1]]
se.uspstf.10 <- SE(svyratio(~lyg.uspstf,~lc.uspstf,uspstf))
uspstf.11 <- svyratio(~lyg.uspstf,~lcds.uspstf,uspstf)[[1]]
se.uspstf.11 <- SE(svyratio(~lyg.uspstf,~lcds.uspstf,uspstf))
uspstf.12 <- svyratio(~fp.uspstf,~lcds.uspstf,uspstf)[[1]]
se.uspstf.12 <- SE(svyratio(~fp.uspstf,~lcds.uspstf,uspstf))
uspstf.13 <- svyratio(~fp.uspstf,~lyg.uspstf,uspstf)[[1]]
se.uspstf.13 <- SE(svyratio(~fp.uspstf,~lyg.uspstf,uspstf))
uspstf.14 <- svytotal(~fp.uspstf,uspstf)[1]
se.uspstf.14 <- SE(svytotal(~fp.uspstf,uspstf))

lcdrat.1 <- svytotal(~lcdrat.eligible,master)[[1]]
se.lcdrat.1 <- SE(svytotal(~lcdrat.eligible,master))
lcdrat.2 <- svyratio(~lc.lcdrat,~cxLCRAT,master)[[1]]
se.lcdrat.2 <- SE(svyratio(~lc.lcdrat,~cxLCRAT,master))
lcdrat.3 <- svyratio(~lcds.lcdrat,~lcds5,master)[[1]]
se.lcdrat.3 <- SE(svyratio(~lcds.lcdrat,~lcds5,master))
lcdrat.4 <- svyratio(~lyg.lcdrat,~lyg,master)[[1]]
se.lcdrat.4 <- SE(svyratio(~lyg.lcdrat,~lyg,master))
lcdrat.5 <- svytotal(~lc.lcdrat,master)[[1]]
se.lcdrat.5 <- SE(svytotal(~lc.lcdrat,master))
lcdrat.6 <- svytotal(~lcds.lcdrat,master)[[1]]
se.lcdrat.6 <- SE(svytotal(~lcds.lcdrat,master))
lcdrat.7 <- svytotal(~lyg.lcdrat,master)[[1]]
se.lcdrat.7 <- SE(svytotal(~lyg.lcdrat,master))
lcdrat.8 <- svyratio(~lcdrat.eligible,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.8 <- SE(svyratio(~lcdrat.eligible,~lcds.lcdrat,lcdrat))
lcdrat.9 <- svyratio(~lcdrat.eligible,~lyg.lcdrat,lcdrat)[[1]]
se.lcdrat.9 <- SE(svyratio(~lcdrat.eligible,~lyg.lcdrat,lcdrat))
lcdrat.10 <- svyratio(~lyg.lcdrat,~lc.lcdrat,lcdrat)[[1]]
se.lcdrat.10 <- SE(svyratio(~lyg.lcdrat,~lc.lcdrat,lcdrat))
lcdrat.11 <- svyratio(~lyg.lcdrat,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.11 <- SE(svyratio(~lyg.lcdrat,~lcds.lcdrat,lcdrat))
lcdrat.12 <- svyratio(~fp.lcdrat,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.12 <- SE(svyratio(~fp.lcdrat,~lcds.lcdrat,lcdrat))
lcdrat.13 <- svyratio(~fp.lcdrat,~lyg.lcdrat,lcdrat)[[1]]
se.lcdrat.13 <- SE(svyratio(~fp.lcdrat,~lyg.lcdrat,lcdrat))
lcdrat.14 <- svytotal(~fp.lcdrat,lcdrat)[1]
se.lcdrat.14 <- SE(svytotal(~fp.lcdrat,lcdrat))

lyg.1 <- svytotal(~lyg.eligible,master)[[1]]
se.lyg.1 <- SE(svytotal(~lyg.eligible,master))
lyg.2 <- svyratio(~lc.lyg,~cxLCRAT,master)[[1]]
se.lyg.2 <- SE(svyratio(~lc.lyg,~cxLCRAT,master))
lyg.3 <- svyratio(~lcds.lyg,~lcds5,master)[[1]]
se.lyg.3 <- SE(svyratio(~lcds.lyg,~lcds5,master))
lyg.4 <- svyratio(~lyg.lyg,~lyg,master)[[1]]
se.lyg.4 <- SE(svyratio(~lyg.lyg,~lyg,master))
lyg.5 <- svytotal(~lc.lyg,master)[[1]]
se.lyg.5 <- SE(svytotal(~lc.lyg,master))
lyg.6 <- svytotal(~lcds.lyg,master)[[1]]
se.lyg.6 <- SE(svytotal(~lcds.lyg,master))
lyg.7 <- svytotal(~lyg.lyg,master)[[1]]
se.lyg.7 <- SE(svytotal(~lyg.lyg,master))
lyg.8 <- svyratio(~lyg.eligible,~lcds.lyg,lyg)[[1]]
se.lyg.8 <- SE(svyratio(~lyg.eligible,~lcds.lyg,lyg))
lyg.9 <- svyratio(~lyg.eligible,~lyg.lyg,lyg)[[1]]
se.lyg.9 <- SE(svyratio(~lyg.eligible,~lyg.lyg,lyg))
lyg.10 <- svyratio(~lyg.lyg,~lc.lyg,lyg)[[1]]
se.lyg.10 <- SE(svyratio(~lyg.lyg,~lc.lyg,lyg))
lyg.11 <- svyratio(~lyg.lyg,~lcds.lyg,lyg)[[1]]
se.lyg.11 <- SE(svyratio(~lyg.lyg,~lcds.lyg,lyg))
lyg.12 <- svyratio(~fp.lyg,~lcds.lyg,lyg)[[1]]
se.lyg.12 <- SE(svyratio(~fp.lyg,~lcds.lyg,lyg))
lyg.13 <- svyratio(~fp.lyg,~lyg.lyg,lyg)[[1]]
se.lyg.13 <- SE(svyratio(~fp.lyg,~lyg.lyg,lyg))
lyg.14 <- svytotal(~fp.lyg,lyg)[1]
se.lyg.14 <- SE(svytotal(~fp.lyg,lyg))

diff.1 <- svytotal(~diff1,master)[1]
se.diff.1 <- SE(svytotal(~diff1,master))
diff.2 <- svytotal(~diff2,master)[1]
se.diff.2 <- SE(svytotal(~diff2,master))
diff.3 <- svyratio(~diff3,~lcds5,master)[[1]]
se.diff.3 <- SE(svyratio(~diff3,~lcds5,master))

diff.4 <- svyratio(~diff4,~age50_59,master)[[1]]
se.diff.4 <- SE(svyratio(~diff4,~age50_59,master))
diff.5 <- svyratio(~diff5,~age60_69,master)[[1]]
se.diff.5 <- SE(svyratio(~diff5,~age60_69,master))
diff.6 <- svyratio(~diff6,~age70plus,master)[[1]]
se.diff.6 <- SE(svyratio(~diff6,~age70plus,master))

est.4 <- rbind(c(uspstf.1,uspstf.5,uspstf.6,uspstf.7,uspstf.8,uspstf.9,uspstf.10,uspstf.11,uspstf.12,uspstf.13,uspstf.14,
                 lcdrat.1,lcdrat.5,lcdrat.6,lcdrat.7,lcdrat.8,lcdrat.9,lcdrat.10,lcdrat.11,lcdrat.12,lcdrat.13,lcdrat.14,
                 lyg.1,lyg.5,lyg.6,lyg.7,lyg.8,lyg.9,lyg.10,lyg.11,lyg.12,lyg.13,lyg.14),
               c(se.uspstf.1,se.uspstf.5,se.uspstf.6,se.uspstf.7,se.uspstf.8,se.uspstf.9,se.uspstf.10,se.uspstf.11,se.uspstf.12,se.uspstf.13,se.uspstf.14,
                 se.lcdrat.1,se.lcdrat.5,se.lcdrat.6,se.lcdrat.7,se.lcdrat.8,se.lcdrat.9,se.lcdrat.10,se.lcdrat.11,se.lcdrat.12,se.lcdrat.13,se.lcdrat.14,
                 se.lyg.1,se.lyg.5,se.lyg.6,se.lyg.7,se.lyg.8,se.lyg.9,se.lyg.10,se.lyg.11,se.lyg.12,se.lyg.13,se.lyg.14)^2)
rat.4 <- rbind(c(uspstf.2,uspstf.3,uspstf.4,
                 lcdrat.2,lcdrat.3,lcdrat.4,
                 lyg.2,lyg.3,lyg.4),
               c(se.uspstf.2,se.uspstf.3,se.uspstf.4,
                 se.lcdrat.2,se.lcdrat.3,se.lcdrat.4,
                 se.lyg.2,se.lyg.3,se.lyg.4)^2)
diff4 <- rbind(c(diff.1,diff.2,diff.3,diff.4,diff.5,diff.6),
                c(se.diff.1,se.diff.2,se.diff.3,se.diff.4,se.diff.5,se.diff.6)^2)


load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod5.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>84)
nhis <- subset(nhis,analypop==1&age>=40&age<=84)  #impute variables for analysis population

persons <- data.frame(age=nhis$age,
                      female=nhis$female,
                      smkyears=nhis$smkyears,
                      qtyears=nhis$qtyears,
                      cpd=nhis$cpd,
                      race=as.numeric(as.character(nhis$race)),
                      emp=nhis$emp,
                      fam.lung.trend=nhis$fam.lung.trend,
                      bmi=nhis$bmi,
                      edu6=nhis$edu6)
predict <- lcrisk(persons,5)
predict2 <- data.frame(matrix(0,nrow(notanaly),7))             
colnames(predict2) <- colnames(predict)
nhis <- rbind(nhis,notanaly)
nhis$predict <- rbind(predict,predict2)
nhis$lcd5 <- nhis$predict[,3]/1000
nhis$LCRAT <- nhis$predict[,5]/1000
nhis$cxLCRAT <- 1.124*nhis$LCRAT
nhis$expected_falsepos <- 0
nhis$expected_falsepos[nhis$cxLCRAT>0] <- as.numeric(predict(polytmod,type="response",newdata=nhis[nhis$cxLCRAT>0,]) %*% c(0,1,2,3))


nhis$lyg <- ifelse(nhis$analypop==1 & nhis$age>=40 & nhis$age<=84, nhis$lyg,0)
nhis$lyg[is.na(nhis$lyg)] <- 0
nhis$lcds5[is.na(nhis$lcds5)] <- 0
nhis$uspstf.eligible[is.na(nhis$uspstf.eligible)] <- 0

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)

uspstf.total <- svytable(~uspstf.eligible,master)[2]

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$lcd5
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$lcd5>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

#select USPSTF size population based on LYG model (40-80)
risk <- nhis$lyg
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
lyg.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lyg.eligible <- ifelse(nhis$lyg>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

nhis$lc.uspstf <- nhis$uspstf.eligible*nhis$cxLCRAT
nhis$lcds.uspstf <- nhis$uspstf.eligible*nhis$lcds5
nhis$lyg.uspstf <- nhis$uspstf.eligible*nhis$lyg
nhis$fp.uspstf <- nhis$uspstf.eligible*nhis$expected_falsepos

nhis$lc.lcdrat <- nhis$lcdrat.eligible*nhis$cxLCRAT
nhis$lcds.lcdrat <- nhis$lcdrat.eligible*nhis$lcds5
nhis$lyg.lcdrat <- nhis$lcdrat.eligible*nhis$lyg
nhis$fp.lcdrat <- nhis$lcdrat.eligible*nhis$expected_falsepos

nhis$lc.lyg <- nhis$lyg.eligible*nhis$cxLCRAT
nhis$lcds.lyg <- nhis$lyg.eligible*nhis$lcds5
nhis$lyg.lyg <- nhis$lyg.eligible*nhis$lyg
nhis$fp.lyg <- nhis$lyg.eligible*nhis$expected_falsepos

nhis$diff1 <- nhis$lyg.lyg - nhis$lyg.lcdrat
nhis$diff2 <- nhis$fp.lcdrat - nhis$fp.lyg
nhis$diff3 <- nhis$lcds.lcdrat - nhis$lcds.lyg

nhis$age50_59 <- ifelse(nhis$age>=50 & nhis$age<=59,1,0)
nhis$age60_69 <- ifelse(nhis$age>=60 & nhis$age<=69,1,0)
nhis$age70plus <- ifelse(nhis$age>=70,1,0)
nhis$age50_59.lcdrat <- nhis$lcdrat.eligible*nhis$age50_59
nhis$age60_69.lcdrat <- nhis$lcdrat.eligible*nhis$age60_69
nhis$age70plus.lcdrat <- nhis$lcdrat.eligible*nhis$age70plus
nhis$age50_59.lyg <- nhis$lyg.eligible*nhis$age50_59
nhis$age60_69.lyg <- nhis$lyg.eligible*nhis$age60_69
nhis$age70plus.lyg <- nhis$lyg.eligible*nhis$age70plus
nhis$diff4 <- nhis$age50_59.lcdrat - nhis$age50_59.lyg
nhis$diff5 <- nhis$age60_69.lcdrat - nhis$age60_69.lyg
nhis$diff6 <- nhis$age70plus.lcdrat - nhis$age70plus.lyg

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
uspstf <- subset(master,uspstf.eligible)
lcdrat <- subset(master,lcdrat.eligible)
lyg <- subset(master,lyg.eligible)

uspstf.1 <- svytotal(~uspstf.eligible,master)[[1]]
se.uspstf.1 <- SE(svytotal(~uspstf.eligible,master))
uspstf.2 <- svyratio(~lc.uspstf,~cxLCRAT,master)[[1]]
se.uspstf.2 <- SE(svyratio(~lc.uspstf,~cxLCRAT,master))
uspstf.3 <- svyratio(~lcds.uspstf,~lcds5,master)[[1]]
se.uspstf.3 <- SE(svyratio(~lcds.uspstf,~lcds5,master))
uspstf.4 <- svyratio(~lyg.uspstf,~lyg,master)[[1]]
se.uspstf.4 <- SE(svyratio(~lyg.uspstf,~lyg,master))
uspstf.5 <- svytotal(~lc.uspstf,master)[[1]]
se.uspstf.5 <- SE(svytotal(~lc.uspstf,master))
uspstf.6 <- svytotal(~lcds.uspstf,master)[[1]]
se.uspstf.6 <- SE(svytotal(~lcds.uspstf,master))
uspstf.7 <- svytotal(~lyg.uspstf,master)[[1]]
se.uspstf.7 <- SE(svytotal(~lyg.uspstf,master))
uspstf.8 <- svyratio(~uspstf.eligible,~lcds.uspstf,uspstf)[[1]]
se.uspstf.8 <- SE(svyratio(~uspstf.eligible,~lcds.uspstf,uspstf))
uspstf.9 <- svyratio(~uspstf.eligible,~lyg.uspstf,uspstf)[[1]]
se.uspstf.9 <- SE(svyratio(~uspstf.eligible,~lyg.uspstf,uspstf))
uspstf.10 <- svyratio(~lyg.uspstf,~lc.uspstf,uspstf)[[1]]
se.uspstf.10 <- SE(svyratio(~lyg.uspstf,~lc.uspstf,uspstf))
uspstf.11 <- svyratio(~lyg.uspstf,~lcds.uspstf,uspstf)[[1]]
se.uspstf.11 <- SE(svyratio(~lyg.uspstf,~lcds.uspstf,uspstf))
uspstf.12 <- svyratio(~fp.uspstf,~lcds.uspstf,uspstf)[[1]]
se.uspstf.12 <- SE(svyratio(~fp.uspstf,~lcds.uspstf,uspstf))
uspstf.13 <- svyratio(~fp.uspstf,~lyg.uspstf,uspstf)[[1]]
se.uspstf.13 <- SE(svyratio(~fp.uspstf,~lyg.uspstf,uspstf))
uspstf.14 <- svytotal(~fp.uspstf,uspstf)[1]
se.uspstf.14 <- SE(svytotal(~fp.uspstf,uspstf))

lcdrat.1 <- svytotal(~lcdrat.eligible,master)[[1]]
se.lcdrat.1 <- SE(svytotal(~lcdrat.eligible,master))
lcdrat.2 <- svyratio(~lc.lcdrat,~cxLCRAT,master)[[1]]
se.lcdrat.2 <- SE(svyratio(~lc.lcdrat,~cxLCRAT,master))
lcdrat.3 <- svyratio(~lcds.lcdrat,~lcds5,master)[[1]]
se.lcdrat.3 <- SE(svyratio(~lcds.lcdrat,~lcds5,master))
lcdrat.4 <- svyratio(~lyg.lcdrat,~lyg,master)[[1]]
se.lcdrat.4 <- SE(svyratio(~lyg.lcdrat,~lyg,master))
lcdrat.5 <- svytotal(~lc.lcdrat,master)[[1]]
se.lcdrat.5 <- SE(svytotal(~lc.lcdrat,master))
lcdrat.6 <- svytotal(~lcds.lcdrat,master)[[1]]
se.lcdrat.6 <- SE(svytotal(~lcds.lcdrat,master))
lcdrat.7 <- svytotal(~lyg.lcdrat,master)[[1]]
se.lcdrat.7 <- SE(svytotal(~lyg.lcdrat,master))
lcdrat.8 <- svyratio(~lcdrat.eligible,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.8 <- SE(svyratio(~lcdrat.eligible,~lcds.lcdrat,lcdrat))
lcdrat.9 <- svyratio(~lcdrat.eligible,~lyg.lcdrat,lcdrat)[[1]]
se.lcdrat.9 <- SE(svyratio(~lcdrat.eligible,~lyg.lcdrat,lcdrat))
lcdrat.10 <- svyratio(~lyg.lcdrat,~lc.lcdrat,lcdrat)[[1]]
se.lcdrat.10 <- SE(svyratio(~lyg.lcdrat,~lc.lcdrat,lcdrat))
lcdrat.11 <- svyratio(~lyg.lcdrat,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.11 <- SE(svyratio(~lyg.lcdrat,~lcds.lcdrat,lcdrat))
lcdrat.12 <- svyratio(~fp.lcdrat,~lcds.lcdrat,lcdrat)[[1]]
se.lcdrat.12 <- SE(svyratio(~fp.lcdrat,~lcds.lcdrat,lcdrat))
lcdrat.13 <- svyratio(~fp.lcdrat,~lyg.lcdrat,lcdrat)[[1]]
se.lcdrat.13 <- SE(svyratio(~fp.lcdrat,~lyg.lcdrat,lcdrat))
lcdrat.14 <- svytotal(~fp.lcdrat,lcdrat)[1]
se.lcdrat.14 <- SE(svytotal(~fp.lcdrat,lcdrat))

lyg.1 <- svytotal(~lyg.eligible,master)[[1]]
se.lyg.1 <- SE(svytotal(~lyg.eligible,master))
lyg.2 <- svyratio(~lc.lyg,~cxLCRAT,master)[[1]]
se.lyg.2 <- SE(svyratio(~lc.lyg,~cxLCRAT,master))
lyg.3 <- svyratio(~lcds.lyg,~lcds5,master)[[1]]
se.lyg.3 <- SE(svyratio(~lcds.lyg,~lcds5,master))
lyg.4 <- svyratio(~lyg.lyg,~lyg,master)[[1]]
se.lyg.4 <- SE(svyratio(~lyg.lyg,~lyg,master))
lyg.5 <- svytotal(~lc.lyg,master)[[1]]
se.lyg.5 <- SE(svytotal(~lc.lyg,master))
lyg.6 <- svytotal(~lcds.lyg,master)[[1]]
se.lyg.6 <- SE(svytotal(~lcds.lyg,master))
lyg.7 <- svytotal(~lyg.lyg,master)[[1]]
se.lyg.7 <- SE(svytotal(~lyg.lyg,master))
lyg.8 <- svyratio(~lyg.eligible,~lcds.lyg,lyg)[[1]]
se.lyg.8 <- SE(svyratio(~lyg.eligible,~lcds.lyg,lyg))
lyg.9 <- svyratio(~lyg.eligible,~lyg.lyg,lyg)[[1]]
se.lyg.9 <- SE(svyratio(~lyg.eligible,~lyg.lyg,lyg))
lyg.10 <- svyratio(~lyg.lyg,~lc.lyg,lyg)[[1]]
se.lyg.10 <- SE(svyratio(~lyg.lyg,~lc.lyg,lyg))
lyg.11 <- svyratio(~lyg.lyg,~lcds.lyg,lyg)[[1]]
se.lyg.11 <- SE(svyratio(~lyg.lyg,~lcds.lyg,lyg))
lyg.12 <- svyratio(~fp.lyg,~lcds.lyg,lyg)[[1]]
se.lyg.12 <- SE(svyratio(~fp.lyg,~lcds.lyg,lyg))
lyg.13 <- svyratio(~fp.lyg,~lyg.lyg,lyg)[[1]]
se.lyg.13 <- SE(svyratio(~fp.lyg,~lyg.lyg,lyg))
lyg.14 <- svytotal(~fp.lyg,lyg)[1]
se.lyg.14 <- SE(svytotal(~fp.lyg,lyg))

diff.1 <- svytotal(~diff1,master)[1]
se.diff.1 <- SE(svytotal(~diff1,master))
diff.2 <- svytotal(~diff2,master)[1]
se.diff.2 <- SE(svytotal(~diff2,master))
diff.3 <- svyratio(~diff3,~lcds5,master)[[1]]
se.diff.3 <- SE(svyratio(~diff3,~lcds5,master))

diff.4 <- svyratio(~diff4,~age50_59,master)[[1]]
se.diff.4 <- SE(svyratio(~diff4,~age50_59,master))
diff.5 <- svyratio(~diff5,~age60_69,master)[[1]]
se.diff.5 <- SE(svyratio(~diff5,~age60_69,master))
diff.6 <- svyratio(~diff6,~age70plus,master)[[1]]
se.diff.6 <- SE(svyratio(~diff6,~age70plus,master))

est.5 <- rbind(c(uspstf.1,uspstf.5,uspstf.6,uspstf.7,uspstf.8,uspstf.9,uspstf.10,uspstf.11,uspstf.12,uspstf.13,uspstf.14,
                 lcdrat.1,lcdrat.5,lcdrat.6,lcdrat.7,lcdrat.8,lcdrat.9,lcdrat.10,lcdrat.11,lcdrat.12,lcdrat.13,lcdrat.14,
                 lyg.1,lyg.5,lyg.6,lyg.7,lyg.8,lyg.9,lyg.10,lyg.11,lyg.12,lyg.13,lyg.14),
               c(se.uspstf.1,se.uspstf.5,se.uspstf.6,se.uspstf.7,se.uspstf.8,se.uspstf.9,se.uspstf.10,se.uspstf.11,se.uspstf.12,se.uspstf.13,se.uspstf.14,
                 se.lcdrat.1,se.lcdrat.5,se.lcdrat.6,se.lcdrat.7,se.lcdrat.8,se.lcdrat.9,se.lcdrat.10,se.lcdrat.11,se.lcdrat.12,se.lcdrat.13,se.lcdrat.14,
                 se.lyg.1,se.lyg.5,se.lyg.6,se.lyg.7,se.lyg.8,se.lyg.9,se.lyg.10,se.lyg.11,se.lyg.12,se.lyg.13,se.lyg.14)^2)
rat.5 <- rbind(c(uspstf.2,uspstf.3,uspstf.4,
                 lcdrat.2,lcdrat.3,lcdrat.4,
                 lyg.2,lyg.3,lyg.4),
               c(se.uspstf.2,se.uspstf.3,se.uspstf.4,
                 se.lcdrat.2,se.lcdrat.3,se.lcdrat.4,
                 se.lyg.2,se.lyg.3,se.lyg.4)^2)
diff5 <- rbind(c(diff.1,diff.2,diff.3,diff.4,diff.5,diff.6),
                c(se.diff.1,se.diff.2,se.diff.3,se.diff.4,se.diff.5,se.diff.6)^2)

inf_est <- inf_eo1(est.1,est.2,est.3,est.4,est.5)
inf_rat <- inf_eo2(rat.1,rat.2,rat.3,rat.4,rat.5)
inf_dif <- inf_diff(diff1,diff2,diff3,diff4,diff5)
inf_dif2 <- inf_diff2(diff1[,3],diff2[,3],diff3[,3],diff4[,3],diff5[,3])