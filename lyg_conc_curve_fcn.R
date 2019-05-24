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


nhis$LCDRAT <- nhis$predict[,3]/1000
nhis$cxLCDRAT <- 0.796*nhis$LCDRAT
nhis$LCRAT <- nhis$predict[,5]/1000
nhis$cxLCRAT <- 1.124*nhis$LCRAT
nhis$expected_falsepos <- 0
nhis$expected_falsepos[nhis$cxLCRAT>0] <- as.numeric(predict(polytmod,type="response",newdata=nhis[nhis$cxLCRAT>0,]) %*% c(0,1,2,3))

nhis$lyg <- ifelse(nhis$age>=40 & nhis$age<=84, nhis$lyg,0)


screenstats <- function(data,group) {
  data <- data[group==1,]
  LCDRAT <- min(data$LCDRAT)
  LCRAT <- min(data$LCRAT)
  LYG <- min(data$lyg)
  cumpeople <- sum(data$adj.wt)                                       #number of ever-smokers screened
  cumLCDsaved <- sum ( (data$LCDRAT - data$cxLCDRAT) * data$adj.wt )  #lung cancer deaths prevented
  cumLYG <- sum(data$lyg * data$adj.wt)                               #life years gained
  cumLCextra <- sum ( (data$cxLCRAT - data$LCRAT) * data$adj.wt )     #additional lung cancers diagnosed
  avgfalsepos <- sum(data$expected_falsepos*data$adj.wt)              #total number of false positive
  cumge1falsepos <- sum(data$predict[,7]/1000 * data$adj.wt)          #number with at least one false positive
  NNSlcd <- cumpeople/cumLCDsaved                                     #NNS to prevent 1 death
  NNSlyg <- cumpeople/cumLYG                                          #NNS to gain 1 year of life
  efficiencylcd <- avgfalsepos/cumLCDsaved                            #number of false positives per LCD prevented
  efficiencylyg <- avgfalsepos/cumLYG                                 #number of false positives per year of life gained
  avgLCDRAT <- sum(data$LCDRAT * data$adj.wt)/cumpeople               #average lung cancer death risk without CT screening
  avgLCRAT <- sum(data$LCRAT * data$adj.wt)/cumpeople                 #average lung cancer risk without CT screening
  avgcxLCDRAT <- sum(data$cxLCDRAT * data$adj.wt)/cumpeople           #average lung cancer death risk with CT screening
  avgcxLCRAT <- sum(data$cxLCRAT * data$adj.wt)/cumpeople             #average lung cancer risk with CT screening
  avgLG <- sum(365.25*data$lyg * data$adj.wt)/cumpeople               #average days of life gained
  c(LCDRAT,LCRAT,LYG,avgLCDRAT,avgcxLCDRAT,avgLCRAT,avgcxLCRAT,avgLG,
    cumpeople,cumLCDsaved,cumLYG,cumLCextra,avgfalsepos,cumge1falsepos,
    NNSlcd,NNSlyg,efficiencylcd,efficiencylyg)
}
guidelines <- matrix(NA,nrow=4,ncol=18)
colnames(guidelines) <- c("LCDRAT","LCRAT","lyg","avgLCDRAT","avgcxLCDRAT","avgLCRAT","avgcxLCRAT","avgLG",
                          "cumpeople","cumLCDsaved","cumLYG","cumLCextra","avgfalsepos", "cumge1falsepos",
                          "NNSlcd","NNSlyg","efficiencylcd","efficiencylyg")

master  <- svydesign(id=~psu, strata=~strata, weights=~adj.wt, data=nhis, nest=TRUE)
master <- subset(master, analypop==1 & age>=40 & age <= 84)
uspstf.total <- svytable(~uspstf.eligible,master)[2]

nhis$smoker <- 1
nhis <- subset(nhis,analypop==1&age>=40&age<=84)
guidelines[1,] <- screenstats(nhis,nhis$smoker)    
guidelines[2,] <- screenstats(nhis,nhis$uspstf.eligible)    

#select USPSTF size population based on LCDRAT model (40-80, highest risk)
risk <- nhis$LCDRAT
x<-cbind(nhis$adj.wt,risk)
y<-x[order(-risk),] # Use order() to sort a matrix
y <- cbind(y,cumsum(y[,1]))
colnames(y) <- c("weight","risk","cumulative weight")

w<-y[,3] # use the cumsum of weights, find the nlst,uspstf,and medicare population sizes in these
dt = data.table(w,val=w)
setattr(dt, "sorted", "w")  # let data.table know that w is sorted
z<-dt[J(uspstf.total), roll = "nearest"] # Find the value nearest to uspstf population size in the data.table
LCDRAT.cutoff.1 <- y[uspstf.index<-match(z[[2]],y[,3]),]
nhis$lcdrat.eligible <- ifelse(nhis$LCDRAT>=y[uspstf.index<-match(z[[2]],y[,3]),][2],1,0)

guidelines[3,] <- screenstats(nhis,nhis$lcdrat.eligible)    


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

guidelines[4,] <- screenstats(nhis,nhis$lyg.eligible)    


##
# Sort NHIS by LCDRAT risk, calculate all screening statistics at each risk threshold
##
nhis<-nhis[order(-nhis$lyg),]
nhis$cumpeople <- cumsum(nhis$adj.wt)
nhis$LCDsaved <- (nhis$LCDRAT - nhis$cxLCDRAT) * nhis$adj.wt
nhis$cumLCDsaved <- cumsum(nhis$LCDsaved)
nhis$cumLYG <- cumsum(nhis$lyg * nhis$adj.wt)
nhis$LCextra <- (nhis$cxLCRAT - nhis$LCRAT) * nhis$adj.wt
nhis$cumLCextra <- cumsum(nhis$LCextra)
nhis$avgfalsepos <- cumsum(nhis$expected_falsepos*nhis$adj.wt)
nhis$cumge1falsepos <- cumsum(nhis$predict[,7]/1000 * nhis$adj.wt)
nhis$NNSlcd <- nhis$cumpeople/nhis$cumLCDsaved
nhis$NNSlyg <- nhis$cumpeople/nhis$cumLYG
nhis$efficiencylcd <- nhis$avgfalsepos/nhis$cumLCDsaved # false positives per LCD prevented
nhis$efficiencylyg <- nhis$avgfalsepos/nhis$cumLYG # false positives per LCD prevented
nhis$avgLCDRAT <- cumsum(nhis$adj.wt * nhis$LCDRAT)/nhis$cumpeople
nhis$avgcxLCDRAT <- cumsum(nhis$adj.wt * nhis$cxLCDRAT)/nhis$cumpeople
nhis$avgLCRAT <- cumsum(nhis$adj.wt * nhis$LCRAT)/nhis$cumpeople
nhis$avgcxLCRAT <- cumsum(nhis$adj.wt * nhis$cxLCRAT)/nhis$cumpeople
nhis$avgLG <- cumsum(365.25*nhis$lyg * nhis$adj.wt)/nhis$cumpeople 
cols <- match(colnames(guidelines),colnames(nhis))

# Include in the NHIS table the row for the risk threshold where we have 40,50,60,70,80,90% of the gainable life
end <- nrow(nhis)
pctpreventables <- c(0.4,0.5,0.6,0.7,0.8,0.9)
pctpreventableindices <- rep(0,6)
for (i in 1:6) {
  print(pctpreventable <- nhis$cumLYG[end]*pctpreventables[i])
  cumLYG <- nhis$cumLYG
  dt2 = data.table(cumLYG,val=cumLYG)
  setattr(dt2, "sorted", "cumLYG")
  z2<-dt2[J(pctpreventable), roll = "nearest"] 
  print(addguidelines <- as.vector(nhis[pctpreventableindices[i]<-match(z2[[2]],nhis$cumLYG),cols]))
  guidelines <- rbind(guidelines,addguidelines)
}
