rm(list=ls(all=TRUE))
library(lcrisks)
library(survey)
library(data.table)
library(ggplot2)

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod1.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>80)
nhis <- subset(nhis,analypop==1&age>=40&age<=80)  #impute variables for analysis population

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
nhis$lcrat <- nhis$predict[,5]/1000
nhis$lcdeath_benefit <- 5*nhis$predict[,4]/1000
nhis$lyg <- ifelse(nhis$age>=40 & nhis$age<=80, nhis$lyg,0)
nhis.1 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod2.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>80)
nhis <- subset(nhis,analypop==1&age>=40&age<=80)  #impute variables for analysis population

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
nhis$lcrat <- nhis$predict[,5]/1000
nhis$lcdeath_benefit <- 5*nhis$predict[,4]/1000
nhis$lyg <- ifelse(nhis$age>=40 & nhis$age<=80, nhis$lyg,0)
nhis.2 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod3.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>80)
nhis <- subset(nhis,analypop==1&age>=40&age<=80)  #impute variables for analysis population

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
nhis$lcrat <- nhis$predict[,5]/1000
nhis$lcdeath_benefit <- 5*nhis$predict[,4]/1000
nhis$lyg <- ifelse(nhis$age>=40 & nhis$age<=80, nhis$lyg,0)
nhis.3 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod4.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>80)
nhis <- subset(nhis,analypop==1&age>=40&age<=80)  #impute variables for analysis population

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
nhis$lcrat <- nhis$predict[,5]/1000
nhis$lcdeath_benefit <- 5*nhis$predict[,4]/1000
nhis$lyg <- ifelse(nhis$age>=40 & nhis$age<=80, nhis$lyg,0)
nhis.4 <- nhis[order(nhis$pid),]

load("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/nhis_imputed_mod5.RData")
nhis$comorbidities <- nhis$emp+nhis$hypertension+nhis$chd+nhis$angina+nhis$heartattack+
  nhis$heartdisease+nhis$stroke+nhis$diab+nhis$bron+nhis$kidney+
  nhis$liver+nhis$prior.cancer+nhis$speceq
notanaly <- subset(nhis,analypop==0|age<40|age>80)
nhis <- subset(nhis,analypop==1&age>=40&age<=80)  #impute variables for analysis population

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
nhis$lcrat <- nhis$predict[,5]/1000
nhis$lcdeath_benefit <- 5*nhis$predict[,4]/1000
nhis$lyg <- ifelse(nhis$age>=40 & nhis$age<=80, nhis$lyg,0)
nhis <- nhis[order(nhis$pid),]
nhis.5 <- nhis

nhis$lcrat <- (nhis.1$lcrat+nhis.2$lcrat+nhis.3$lcrat+nhis.4$lcrat+nhis.5$lcrat)/5
nhis$lcdeath_benefit <- (nhis.1$lcdeath_benefit+nhis.2$lcdeath_benefit+nhis.3$lcdeath_benefit+nhis.4$lcdeath_benefit+nhis.5$lcdeath_benefit)/5
nhis$lyg <- (nhis.1$lyg+nhis.2$lyg+nhis.3$lyg+nhis.4$lyg+nhis.5$lyg)/5

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas1.jpeg",width=8,height=8,units='in',res=900)
nhis <- subset(nhis, analypop==1 & age>=40 & age <= 80)
plot(cumsum(nhis$adj.wt[order(nhis$lcrat)])/sum(nhis$adj.wt),cumsum(nhis$adj.wt*nhis$lcrat[order(nhis$lcrat)])/sum(nhis$adj.wt*nhis$lcrat),
     type="l",lwd=2,col="red",ylim=c(0.035,0.9525),xlim=c(0.035,0.9525),
     main="Lorenz curve",cex.main=.95,
     ylab="Cumulative frequency",xlab="Cumulative frequency of smokers")
lines(cumsum(nhis$adj.wt[order(nhis$lyg)])/sum(nhis$adj.wt),cumsum(nhis$adj.wt*nhis$lyg[order(nhis$lyg)])/sum(nhis$adj.wt*nhis$lyg))
abline(a=0,b=1)
text(x=.4,y=.15,"Life gained based")
text(x=.7,y=.1,"Risk-based",col="red")
dev.off()

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas2.jpeg",width=11.5,height=8,units='in',res=900)
ggplot(nhis, aes(x=100*nhis$lcrat,y=nhis$lyg,size=adj.wt)) +
      geom_point(shape=21) +
      ggtitle("Lung cancer risk and Life gained") +
      labs(x = "5-year lung cancer risk", y = "Years of life gained")
dev.off()

nhis$comorbidities_cat <- ifelse(nhis$comorbidities>=5,5,nhis$comorbidities)
p1 <- ggplot(nhis, aes(x=comorbidities_cat,y=365.25*lyg)) +
  geom_boxplot(outlier.size=.5,aes(group=cut_width(comorbidities_cat,1),weight=adj.wt)) +
  geom_hline(yintercept=14.8) + 
  ggtitle("Life gained from screening among ever-smokers by comorbidities") +
  geom_text(x=42,y=17,label="Screening threshold=14.8 days") + 
  labs(x = "Number of comorbidities", y = "Days of life gained")

p2 <- ggplot(nhis, aes(x=comorbidities_cat,y=100*lcrat)) +
  geom_boxplot(outlier.size=.5,aes(group=cut_width(comorbidities_cat,1),weight=adj.wt)) +
  geom_hline(yintercept=2.04) + 
  ggtitle("Lung cancer risk among ever-smokers by comorbidities") +
  geom_text(x=42,y=2.5,label="Screening threshold=2.04%") + 
  labs(x = "Number of comorbidities", y = "% risk of lung cancer within 5 years")

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas3.jpeg",width=11.5,height=8,units='in',res=900)
multiplot(p1, p2, cols=1)
dev.off()

p1 <- ggplot(nhis, aes(x=age,y=365.25*lyg)) +
  geom_boxplot(outlier.size=.5,aes(group=cut_width(age,1),weight=adj.wt,size=adj.wt)) +
  geom_hline(yintercept=14.8) + 
  ggtitle("Life gained from screening among ever-smokers by age") +
  geom_text(x=42,y=17,label="Screening threshold=14.8 days") + 
  labs(x = "Age", y = "Days of life gained")

p2 <- ggplot(nhis, aes(x=age,y=100*lcrat)) +
  geom_boxplot(outlier.size=.5,aes(group=cut_width(age,1),weight=adj.wt,size=adj.wt)) +
  geom_hline(yintercept=2.04) + 
  ggtitle("Lung cancer risk among ever-smokers by age") +
  geom_text(x=42,y=2.5,label="Screening threshold=2.04%") + 
  labs(x = "Age", y = "% risk of lung cancer within 5 years")

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas4.jpeg",width=11.5,height=8,units='in',res=900)
multiplot(p1, p2, cols=1)
dev.off()

p1 <- ggplot(nhis, aes(x=100*lcrat,y=365.25*lyg)) +
  geom_boxplot(outlier.size=.5,aes(group=cut_width(100*lcrat,.5),weight=adj.wt)) +
  geom_hline(yintercept=14.8) + 
  ggtitle("Life gained from screening among ever-smokers by lung cancer risk") +
  geom_text(x=42,y=17,label="Screening threshold=14.8 days") + 
  labs(x = "% risk of lung cancer within 5-years", y = "Days of life gained")

jpeg("~/Desktop/Lung cancer/lrisk/other/lifeyearsgained/plot_ideas6.jpeg",width=11.5,height=8,units='in',res=900)
p1
dev.off()