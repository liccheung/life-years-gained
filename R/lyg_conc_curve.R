#marginal selection of USPSTF, Risk-based, and Life-years gained
rm(list=ls(all=TRUE))
library(lcrisks)
library(survey)
library(data.table)
load("~/Desktop/Lung cancer/lrisk/other/jama/polytomousmodel.RData")
load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod1.RData")
source(file="~/Desktop/Lung cancer/lrisk/prog/lifeyearsgained/lyg_conc_curve_fcn.R")
guidelines.1 <- guidelines
nhis.1 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod2.RData")
source(file="~/Desktop/Lung cancer/lrisk/prog/lifeyearsgained/lyg_conc_curve_fcn.R")
guidelines.2 <- guidelines
nhis.2 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod3.RData")
source(file="~/Desktop/Lung cancer/lrisk/prog/lifeyearsgained/lyg_conc_curve_fcn.R")
guidelines.3 <- guidelines
nhis.3 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod4.RData")
source(file="~/Desktop/Lung cancer/lrisk/prog/lifeyearsgained/lyg_conc_curve_fcn.R")
guidelines.4 <- guidelines
nhis.4 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod5.RData")
source(file="~/Desktop/Lung cancer/lrisk/prog/lifeyearsgained/lyg_conc_curve_fcn.R")
guidelines.5 <- guidelines
nhis.5 <- nhis[order(nhis$pid),]

guidelines <- (guidelines.1+guidelines.2+guidelines.3+guidelines.4+guidelines.5)/5
write.table(guidelines,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/conc_curve.csv", sep=",", row.names=TRUE)

LCRAT <- (nhis.1$LCRAT+nhis.2$LCRAT+nhis.3$LCRAT+nhis.4$LCRAT+nhis.5$LCRAT)/5
lyg <- (nhis.1$lyg+nhis.2$lyg+nhis.3$lyg+nhis.4$lyg+nhis.5$lyg)/5
order1 <- order(-lyg)
lyg <- lyg[order1]
adj.wt <- nhis.1$adj.wt[order1]
cumLYG <- cumsum(lyg*adj.wt)
cumpeople <- cumsum(adj.wt)
max(cumLYG[cumpeople<=guidelines[2,9]])
xmatch <- c(max(cumpeople[cumLYG<=.4*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.5*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.6*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.7*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.8*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.9*guidelines[1,11]]),
            max(cumpeople))
order2 <- order(-LCRAT[order1])
lyg2 <- lyg[order2]
adj.wt2 <- adj.wt[order2]
cumLYG2 <- cumsum(lyg2*adj.wt2)
cumpeople2 <- cumsum(adj.wt2)
max(cumLYG2[cumpeople2<=guidelines[2,9]])

guidelines[8,7]*guidelines[8,9]  #screening 16.6 mil detects this number of lung cancers 
(guidelines[8,7]*guidelines[8,9])/(guidelines[1,7]*guidelines[1,9])  

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/conc_curve.jpeg",width=8,height=8,units='in',res=900)
plot(cumpeople/1000000,100*cumLYG/guidelines[1,11],type="l",
     main="Outcomes From Different Life Gained-Based Thresholds \nin US Ever-Smokers Aged 40 to 80 Years",
     xaxt="n",xlab="Screened Ever-Smokers Aged 40-80y, millions",
     ylim=c(30,100),ylab="Gainable Life, %")
axis(side=1,at=c(0,5,10,15,20,25,30,35,40,45,50,55,60),labels=NULL)
points(guidelines[2,9]/1000000,100*guidelines[2,11]/guidelines[1,11],pch=16,cex=.75)
points(guidelines[3,9]/1000000,100*max(cumLYG2[cumpeople2<=guidelines[3,9]])/guidelines[1,11],cex=.5)
#need to fix these lines to be based on cumpeople at 40,50,60,70,80,90 for average
segments(xmatch/1000000,c(40,50,60,70,80,90,100),xmatch/1000000,rep(0,5),col="lightblue")
abline(h=c(30,40,50,60,70,80,90,100),col="grey")
dev.off()