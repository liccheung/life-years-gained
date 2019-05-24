load("best_AIC_models.RData")
load("polytomousmodel.RData")

risk.kovalchik <-
function (begin, end, newdata, coxph1, ...) 
{
  c.coxph.risk <- function (coxph, ...) 
  {
    models <- list(coxph)
    if (!missing(..1)) 
      models <- c(models, list(...))
    models
  }  

  coxph.relrisk.uncentered <- function (coxph.object, newdata) 
  {
    center <- coxph.object$means %*% coef(coxph.object)
    if (!missing(newdata)) 
      lp <- predict(coxph.object, newdata, type = "lp")
    else lp <- coxph.object$linear.predictor
    exp(lp + center)
  }  
  
  projection.relrisk <- function (object, data) 
  {
    if (is.numeric(object)) 
      return(object)
    else if (class(object) == "coxph") 
      if (missing(data)) 
        return(coxph.relrisk.uncentered(object))
    else return(coxph.relrisk.uncentered(object, data))
    else stop(cat("No method for class", class(object)))
  }  

    #calculate absolute risk given hazards/survival and relative risks
    risk.fixed.interval <- function(H, RR) {
        if (!is.list(H)) 
            0
        else {
            absrisk <- H[[1]]$surv^RR[1] * H[[1]]$haz * RR[1]
            for (i in 2:length(H)) {
                absrisk <- absrisk * (H[[i]]$surv^RR[i])
            }
            sum(absrisk)
        }
    }
    
    models <- c.coxph.risk(coxph1, ...)
    
    #check for missing
    which.kept <- complete.cases(newdata)
    if (!all(which.kept)) {
        warning("Missing cases excluded.")
        if (length(begin) > 1) {
            begin <- begin[which.kept]
            end <- end[which.kept]
        }
        newdata <- newdata[which.kept, ]
    }

    #calculate relative risk for each subject and store as list
    rr <- sapply(models, projection.relrisk, data = newdata)

    if (is.matrix(rr)) { rr.list <- lapply(1:nrow(rr), function(x) rr[x, ])
    } else { rr.list <- list(rr)}

    #estimate risk
    if (length(begin) == 1) {    
       AllVars <- unique(unlist(sapply(models, function(x) all.vars(x$formula))))
        in.interval <- function(x, begin, end) x >= begin & x <= end
        if ("lung.cancer.death" %in% AllVars){
          H <- list(models[[1]]$basehaz[in.interval(models[[1]]$basehaz$time, begin, end),], 
          models[[2]]$basehaz_LCDRAT[in.interval(models[[1]]$basehaz$time, begin, end),])
        } else if ("case" %in% AllVars){
          H <- list(models[[1]]$basehaz[in.interval(models[[1]]$basehaz$time, begin, end),], 
          models[[2]]$basehaz_LCRAT[in.interval(models[[1]]$basehaz$time, begin, end),])        
        }
        #calculate absolute risk given hazards/survival for average covariate values and relative risks
        risks <- mapply(risk.fixed.interval, RR = rr.list, MoreArgs = list(H = H))
    } else {
        risks <- mapply(risk, begin = begin, end = end, RR = rr.list, MoreArgs = list(models = models))
    }
    risks
}



#' Wrapper to the cut() function
#'
#' A wrapper to the cut() function, so that you can automatically break into quantiles as
#' the default behavior, otherwise if the breakpoints are included, then just break on those.
#' In all cases, include.lowest is set to True      
#' @export
cuts <- function(data,npieces,simple.labels=T,...) {

  if (length(npieces)==0 | any(is.na(npieces)) | any(!is.numeric(npieces)))
    stop("npieces must be a numeric scalar or a vector of breakpoints")
  
  if (length(npieces)==1)
    # Just break into quantiles.  Use quantile labelling: Q1,Q2,Q3,etc. if you want
    if (simple.labels)
      cutdata <- cut(data,breaks=quantile(data,seq(0,1,1/npieces),na.rm=T,...),
                     include.lowest=T,labels=paste("Q",1:npieces,sep=""),...)
    else
      cutdata <- cut(data,breaks=quantile(data,seq(0,1,1/npieces),na.rm=T,...),
                     include.lowest=T,...)
  else
    # Break into pieces as specified by the breakpoints
    if (simple.labels)
      cutdata <- cut(data,breaks=npieces,include.lowest=T,
                     labels=paste("Q",1:(length(npieces)-1),sep=""),...)
    else
      cutdata <- cut(data,breaks=npieces,include.lowest=T,...)
  
  return(cutdata)
}



#' Lung Cancer Death Risk Predictor
#'
#' In both the absence and presence of screening, the R package calculates individual risks 
#' of lung cancer and lung cancer death based on covariates: age, education, sex, race, 
#' smoking intensity/duration/quit-years, Body Mass Index, family history of lung-cancer, 
#' and self-reported emphysema.  In the presence of CT screening akin to the NLST 
#' (3 yearly screens, 5 years of follow-up), it uses the covariates to estimate risk of 
#' false-positive CT screen as well as the reduction in risk of lung cancer death and 
#' increase in risk of lung cancer screening.  
#'
#' @import coxph.risk
#'
#' @section Warning:
#'  VGAM is a required dependency of this package.
#'  VGAM may automatically be installed the first time this package is used. 
#'
#' @section Model Objects in Package:
#'  \itemize{
#'  \item LCDRAT - model for lung cancer death in absence of screening;
#'  \item LCRAT - model for lung cancer incidence in absence of screening;
#'  \item cox.death - model for deaths from causes other than lung cancer;
#'  \item polytmod - polytomous model for false positive CT lung screens.
#'  }
#' 
#' @param x A numeric matrix containing individuals' covariates for the model.  
#'  Covariates should be in the following column and format:
#'
#'  \itemize{
#'  \item column 1 - current age (numeric);
#'  \item column 2 - gender (1=Female, 0=Male);
#'  \item column 3 - years smoked (numeric);
#'  \item column 4 - years quit (numeric or NA);
#'  \item column 5 - cigarettes per day (numeric or NA);
#'  \item column 6 - race (0=Non-hispanic white,
#'                   1=Non-hispanic Black/African American, 
#'                   2=Hispanic, 
#'                   3=Other Ethnicity);
#'  \item column 7 - lung disease (1=COPD or Emphysema, 0=No COPD or Emphysema);
#'  \item column 8 - number of parents with lung cancer (0,1,2);
#'  \item column 9 - bmi;
#'  \item column 10 - highest education level (1=<12 grade, 
#'                                       2=HS graduate, 
#'                                       3=post hs, no college, 
#'                                       4=associate degree/some college, 
#'                                       5=bachelors degree, 
#'                                       6=graduate school);
#'  }
#' @param y Number of years to calculate risks for (numeric, max of 10).
#'
#' @return A numeric matrix containing individuals' predictions:
#'
#'  \itemize{
#'  \item column 1 - An indicator variable for whether the individual is eligible 
#'              for CT lung screening according to 
#'              US Preventive Services Task Force (USPSTF) recommendations.
#'  \item column 2 - Number of years predictions are for.
#'  \item column 3 - Among 1000 people in the US with this risk-factor profile, 
#'              this is the number who will die from lung cancer 
#'              if they do not attend screening.
#'  \item column 4 - In the NLST, those who underwent 3 rounds of annual CT 
#'              screening had their risk reduced by 20 percent.  
#'              Therefore, among those who would have died from lung cancer,
#'              this is the number who will not die from lung cancer death, 
#;              if they undergo 3 yearly CT lung screens as in the NLST.
#'  \item column 5 - Among 1000 people in the US with this risk-factor profile, 
#'              this is the number who will be diagnosed with lung cancer 
#'              if they do not attend screening (LCRAT).
#'  \item column 6 - In the NLST, those who underwent CT screening had 12.4 percent 
#'              more lung cancer diagnosed, all of which require treatment.  
#'              Therefore, among 1000 people with this risk-factor profile, 
#'              this is the number of extra lung cancer that would be 
#'              diagnosed, if they undergo 3 yearly CT lung screens as in the NLST.
#'  \item column 7 - Out of 1000 NLST participants with this risk profile, 
#'              this is the number who had at least one false-positive 
#'              CT screen out of 3 screens.
#'  }
#'
#' @author Li C. Cheung, \email{li.cheung@nih.gov}, Stephanie A. Kovalchik, Hormuzd A. Katki
#'
#' @references Katki HA, Kovalchik SA, Berg CD, Cheung LC, Chaturvedi AK.
#'             Development and validation of risk models to select ever-smokers
#'             for CT lung cancer screening. JAMA. 2016; 315(21):2300-2311.
#'             doi: 10.1001/jama.2016.6255.
#'
#' @export
#' @examples
#' age <- c(66,58,75,72,56)
#' bmi <- c(23,28,26,27,24)
#' cpd <- c(36,36,40,24,40)
#' emp <- c(0,1,1,0,1)
#' fam.lung.trend <- c(0,2,0,2,0)
#' female <- c(0,1,0,1,0)
#' smkyears <- c(43,37,45,42,29)
#' qtyears <- c(NA,NA,9,6,6)
#' race <- c(0,1,2,2,3)
#' edu6 <- c(3,5,4,5,5)
#' years <- 5
#' 
#' persons <- cbind(age,
#'                  female,
#'                  smkyears,
#'                  qtyears,
#'                  cpd,
#'                  race,
#'                  emp,
#'                  fam.lung.trend,
#'                  bmi,
#'                  edu6)
#' 
#' persons_predictions <- lcrisk(persons,years)
#' persons_predictions

lcrisk <- function(x,y) {
   if (!require("VGAM")) {
     install.packages("VGAM",repos="http://watson.nci.nih.gov/cran_mirror/")
   } 
   library("VGAM")

   x <- as.matrix(x)
   
   age <- x[,1]
   female <- x[,2]
   smkyears <- x[,3]
   qtyears <- x[,4]
   qtyears <- ifelse(is.na(qtyears),0,qtyears)
   cpd <- x[,5]
   race <- as.factor(x[,6])
   emp <- x[,7]
   fam.lung.trend <- x[,8]
   bmi <- x[,9]
   edu6 <- x[,10]
   pkyr.cat <- smkyears*cpd/20
   
   covar <- data.frame(age=age,
                       bmi=bmi,
                       cpd=cpd,
                       emp=emp,
                       fam.lung.trend=fam.lung.trend,
                       female=female,
                       qtyears=qtyears,
                       smkyears=smkyears,
                       race=race,
                       edu6=edu6,
                       pkyr.cat=pkyr.cat)
      
   uspstf_elig <- ifelse(age>=50 & age<=80 & pkyr.cat>=30 & qtyears <= 15,1,0)
   
   #calculate predicted risks of lung cancer death within y years based on LCDRAT with competing risk of death
   covar$LCDRAT <- risk.kovalchik(0, y, covar[,1:11], LCDRAT, cox.death)
   
   #calculate predicted risks of lung cancer death using chest xray
   covar$cxLCDRAT <- 0.796*covar$LCDRAT
   
   #calculate predicted risks of lung cancer incidence within y years based on LCRAT with competing risk of death
   covar$LCRAT <- risk.kovalchik(0, y, covar[,1:11], LCRAT, cox.death)
   covar$LCRAT <- pmax(covar$LCRAT,covar$LCDRAT)
   
   #calculate predicted risks of lung cancer incidence using chest xray
   covar$cxLCRAT <- 1.124*covar$LCRAT
   
   #calculate probability of Z, number of false positives, taking values of 0, 1, 2, or 3 
   #responses are reported as log(P(Z=1)/P(Z=0)), log(P(Z=2)/P(Z=0)), and log(P(Z=3)/P(Z=0))
   prob_numfalsepos <- predict(polytmod,type="response",newdata=covar)
   
   #solving for P(Z=0), P(Z=1), P(Z=2), and P(Z=3), we have:
   covar$prob_0falsepos <- prob_numfalsepos[,1]
   covar$prob_1falsepos <- prob_numfalsepos[,2]
   covar$prob_2falsepos <- prob_numfalsepos[,3]
   covar$prob_3falsepos <- prob_numfalsepos[,4]
   covar$expected_falsepos <- covar$prob_1falsepos + 2*covar$prob_2falsepos + 3*covar$prob_3falsepos
   
   ### US preventive service task force eligible - uspstf_elig
   
   ### number out of 1000 with this risk profile who will die from lung cancer death in y years without screening - 1000*covar$LCDRAT  
   
   ### Reduction in y years lung cancer death risk due to screening - 1000*(covar$LCDRAT - covar$cxLCDRAT)

   ### chance of lung cancer diagnosis without screening - 1000*covar$LCRAT
   
   ### increase in lung cancer diagnosis due to screening - 1000*(covar$cxLCRAT - covar$LCRAT)
   
   ### chance of false positive CT lung screen - 1000*(1-covar$prob_0falsepos)
   
   predicted <- data.frame(uspstf_elig,
                           rep(y,length(uspstf_elig)),
                           1000*covar$LCDRAT,
                           1000*(covar$LCDRAT - covar$cxLCDRAT),
                           1000*covar$LCRAT,
                           1000*(covar$cxLCRAT - covar$LCRAT),
                           1000*(1-covar$prob_0falsepos))
                           
   colnames(predicted) <- c("USPSTF eligible",
                            "Number of years predictions are for",
                            "Number of lung cancer deaths per 1000 (LCDRAT)",
                            "Screening reduced lung cancer deaths per 1000",
                            "Number with lung cancer diagnosed per 1000 (LCRAT)",
                            "Screening increase lung cancer diagnosis per 1000",
                            "False-positive CT lung screens per 1000")
   predicted                         
}
