#marginal selection of USPSTF, Risk-based, and Life-years gained
rm(list=ls(all=TRUE))
library(lcrisks)
library(survey)
library(data.table)
load("~/Desktop/Lung cancer/lrisk/other/jama/polytomousmodel.RData")
load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod1.RData")
#change to 58.5% of gainable life for same NNS as USPSTF
source(file="~/Desktop/Lung cancer/lrisk/prog/lifeyearsgained/lyg_conc_curve_fcn.R")
guidelines.1 <- guidelines
nhis.1 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod2.RData")
source(file="~/Desktop/Lung cancer/lrisk/prog/lifeyearsgained/lyg_conc_curve_fcn.R")
guidelines.2 <- guidelines
nhis.2 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod3.RData")
source(file="~/Desktop/Lung cancer/lrisk/prog/lifeyearsgained/lyg_conc_curve_fcn.R")
guidelines.3 <- guidelines
nhis.3 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod4.RData")
source(file="~/Desktop/Lung cancer/lrisk/prog/lifeyearsgained/lyg_conc_curve_fcn.R")
guidelines.4 <- guidelines
nhis.4 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/v6/nhis_imputed_mod5.RData")
source(file="~/Desktop/Lung cancer/lrisk/prog/lifeyearsgained/lyg_conc_curve_fcn.R")
guidelines.5 <- guidelines
nhis.5 <- nhis[order(nhis$pid),]

guidelines <- (guidelines.1+guidelines.2+guidelines.3+guidelines.4+guidelines.5)/5
write.table(guidelines,file="~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/conc_curve.csv", sep=",", row.names=TRUE)

#sensitivity LCs
guidelines[,6]*guidelines[,9]/(guidelines[1,6]*guidelines[1,9])
#sensitivity LCDs
guidelines[,4]*guidelines[,9]/(guidelines[1,4]*guidelines[1,9])
#sensitivity LYG
guidelines[,11]/guidelines[1,11]

#benefits LC detected
guidelines[,7]*guidelines[,9]
#benefits lives saved
guidelines[,10]
#benefits years gained
guidelines[,11]
#NNS lives saved
guidelines[,15]
#NNS 10 years gained
10*guidelines[,16]
#Years gained per lung cancer detected
guidelines[,11]/(guidelines[,7]*guidelines[,9])
#Years gained per prevented LCD
guidelines[,11]/guidelines[,10]
#false positives per LCD prevented
guidelines[,17]
#false positives per 10 LYG
10*guidelines[,18]


LCDRAT <- (nhis.1$LCDRAT+nhis.2$LCDRAT+nhis.3$LCDRAT+nhis.4$LCDRAT+nhis.5$LCDRAT)/5
lyg <- (nhis.1$lyg+nhis.2$lyg+nhis.3$lyg+nhis.4$lyg+nhis.5$lyg)/5
LCRAT <- (nhis.1$LCRAT+nhis.2$LCRAT+nhis.3$LCRAT+nhis.4$LCRAT+nhis.5$LCRAT)/5
order1 <- order(-lyg)
lyg <- lyg[order1]
adj.wt <- nhis.1$adj.wt[order1]
cumLYG <- cumsum(lyg*adj.wt)
LCDRAT <- LCDRAT[order1]
cumLCDS <- .204*cumsum(LCDRAT*adj.wt)
cumpeople <- cumsum(adj.wt)
max(cumLYG[cumpeople<=guidelines[2,9]])
age <- nhis$age[order1]
cumlt50 <- cumsum(adj.wt*I(age<50))
LCRAT <- LCRAT[order1]
cumLC <- .204*cumsum(LCRAT*adj.wt)
xmatch <- c(min(365.25*lyg[cumLYG<=.1*guidelines[1,11]]),
            min(365.25*lyg[cumLYG<=.2*guidelines[1,11]]),
            min(365.25*lyg[cumLYG<=.3*guidelines[1,11]]),
            min(365.25*lyg[cumLYG<=.4*guidelines[1,11]]),
            min(365.25*lyg[cumLYG<=.5*guidelines[1,11]]),
            min(365.25*lyg[cumLYG<=.6*guidelines[1,11]]),
            min(365.25*lyg[cumLYG<=.7*guidelines[1,11]]),
            min(365.25*lyg[cumLYG<=.8*guidelines[1,11]]),
            min(365.25*lyg[cumLYG<=.9*guidelines[1,11]]))

xmatch2 <- c(max(cumpeople[cumLYG<=.1*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.2*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.3*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.4*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.5*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.6*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.7*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.8*guidelines[1,11]]),
            max(cumpeople[cumLYG<=.9*guidelines[1,11]]))

ymatch2 <- c(max(cumLYG[cumLYG<=.1*guidelines[1,11]]),
             max(cumLYG[cumLYG<=.2*guidelines[1,11]]),
             max(cumLYG[cumLYG<=.3*guidelines[1,11]]),
             max(cumLYG[cumLYG<=.4*guidelines[1,11]]),
             max(cumLYG[cumLYG<=.5*guidelines[1,11]]),
             max(cumLYG[cumLYG<=.6*guidelines[1,11]]),
             max(cumLYG[cumLYG<=.7*guidelines[1,11]]),
             max(cumLYG[cumLYG<=.8*guidelines[1,11]]),
             max(cumLYG[cumLYG<=.9*guidelines[1,11]]))

ymatch3 <- c(max(cumLCDS[cumLYG<=.1*guidelines[1,11]]),
             max(cumLCDS[cumLYG<=.2*guidelines[1,11]]),
             max(cumLCDS[cumLYG<=.3*guidelines[1,11]]),
             max(cumLCDS[cumLYG<=.4*guidelines[1,11]]),
             max(cumLCDS[cumLYG<=.5*guidelines[1,11]]),
             max(cumLCDS[cumLYG<=.6*guidelines[1,11]]),
             max(cumLCDS[cumLYG<=.7*guidelines[1,11]]),
             max(cumLCDS[cumLYG<=.8*guidelines[1,11]]),
             max(cumLCDS[cumLYG<=.9*guidelines[1,11]]))

ymatch4 <- c(max(cumLC[cumLYG<=.1*guidelines[1,11]]),
             max(cumLC[cumLYG<=.2*guidelines[1,11]]),
             max(cumLC[cumLYG<=.3*guidelines[1,11]]),
             max(cumLC[cumLYG<=.4*guidelines[1,11]]),
             max(cumLC[cumLYG<=.5*guidelines[1,11]]),
             max(cumLC[cumLYG<=.6*guidelines[1,11]]),
             max(cumLC[cumLYG<=.7*guidelines[1,11]]),
             max(cumLC[cumLYG<=.8*guidelines[1,11]]),
             max(cumLC[cumLYG<=.9*guidelines[1,11]]))

ymatch4 <- c(max(cumlt50[cumLYG<=.1*guidelines[1,11]]),
             max(cumlt50[cumLYG<=.2*guidelines[1,11]]),
             max(cumlt50[cumLYG<=.3*guidelines[1,11]]),
             max(cumlt50[cumLYG<=.4*guidelines[1,11]]),
             max(cumlt50[cumLYG<=.5*guidelines[1,11]]),
             max(cumlt50[cumLYG<=.6*guidelines[1,11]]),
             max(cumlt50[cumLYG<=.7*guidelines[1,11]]),
             max(cumlt50[cumLYG<=.8*guidelines[1,11]]),
             max(cumlt50[cumLYG<=.9*guidelines[1,11]]))

lag(xmatch2)-xmatch

order2 <- order(-LCDRAT[order1])
lyg2 <- lyg[order2]
adj.wt2 <- adj.wt[order2]
cumLYG2 <- cumsum(lyg2*adj.wt2)
cumpeople2 <- cumsum(adj.wt2)
max(cumLYG2[cumpeople2<=guidelines[2,9]])

guidelines[8,7]*guidelines[8,9]  #screening 16.6 mil detects this number of lung cancers 
(guidelines[8,7]*guidelines[8,9])/(guidelines[1,7]*guidelines[1,9])  

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/conc_curve_v6.jpeg",width=8,height=8,units='in',res=900)
plot(365.25*lyg,100*cumLYG/guidelines[1,11],type="l",
     xlim=c(1.5,30),xaxt="n",xlab="Days of life-gained thresholds",
     ylim=c(5,100),yaxt="n",ylab="Sensitivity, % of all life-gained")
axis(side=1,at=round(xmatch[-1],digits=1))
axis(side=2,at=c(10,20,30,40,50,60,70,80,90,100))
#need to fix these lines to be based on cumpeople at 40,50,60,70,80,90 for average
segments(xmatch,c(10,20,30,40,50,60,70,80,90,100),xmatch,rep(0,5),col="lightblue")
abline(h=c(10,20,30,40,50,60,70,80,90,100),col="grey")
dev.off()