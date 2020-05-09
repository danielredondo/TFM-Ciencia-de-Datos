nordpred.estimate <- function(cases,pyr,noperiod,startestage,linkfunc="power5") {
  Rplatform <- exists("R.Version")
  if ( dim(cases)[1]!=18 || dim(pyr)[1]!=18 ) {
    stop("\"cases\" and \"pyr\" must have data for 18 age groups")	
  }
  if ( dim(cases)[2]>dim(pyr)[2]) {
    stop("\"pyr\" must include information about all periods in \"cases\"")	
  }
  if (dim(pyr)[2]==dim(cases)[2]) {
    stop("\"pyr\" must include information on future rates (E.g. must be greater than \"cases\")")	
  }
  if ((dim(pyr)[2]-dim(cases)[2])>5) {
    stop("Package can not predict more than 5 periods (given by sizes of \"pyr\" and \"cases\")")	
  }
  if ((dim(cases)[2]-noperiod)<0) {
    stop("More periods specified in \"noperiod\" than columns in \"cases\"")   	
  }
  if (noperiod<3) {
    stop("\"noperiod\" must be 3 or larger to get enough for estimating")   	
  }  
  dnoperiods <- dim(cases)[2]
  dnoagegr   <- 18
  ageno    <- rep(1:dnoagegr,dnoperiods)
  periodno <- sort(rep(1:dnoperiods,dnoagegr))
  cohort   <- max(ageno)-ageno+periodno
  y        <- c(as.matrix(pyr[,1:dnoperiods]))
  apcdata <- data.frame(Cases= c(as.matrix(cases)),Age=ageno,Cohort=cohort,Period=periodno,y=y)
  apcdata <- apcdata[apcdata$Age>=startestage,]
  apcdata <- apcdata[apcdata$Period>(dnoperiods-noperiod),]
  if (Rplatform) {
    y <- apcdata$y
    power5link <- poisson()
    power5link$link <- "0.2 root link Poisson family"
    power5link$linkfun <- function(mu)  { (mu/y)^0.2 }
    power5link$linkinv <- function(eta) { pmax(.Machine$double.eps, y*eta^5) }
    power5link$mu.eta <- function(eta)  { pmax(.Machine$double.eps, 5*y*eta^4) }
  } else {
    y <<- apcdata$y
    power5link         <- poisson()
    power5link$link    <- function(mu)  { (mu/y)^0.2 }
    power5link$inverse <- function(eta) { y*eta^(1/0.2) }
    power5link$deriv   <- function(mu)  { 0.2*(1/y)*(mu/y)^(0.2 - 1) }
  }
  options(contrasts=c("contr.treatment","contr.poly"))
  if (linkfunc=="power5") {
    res.glm <- glm(Cases~as.factor(Age)+Period+as.factor(Period)+as.factor(Cohort) -1,family=power5link,data=apcdata)
  } else  if (linkfunc=="poisson") {
    res.glm <- glm(Cases~as.factor(Age)+Period+as.factor(Period)+as.factor(Cohort)+ offset(log(y)) -1,family=poisson(),data=apcdata)  	
  } else {
    stop("Unknown \"linkfunc\"")	
  }
  if (Rplatform) {
    pvalue <- 1-pchisq(res.glm$deviance,res.glm$df.residual)
  } else {
    pvalue <- 1-pchisq(res.glm$deviance,res.glm$df) 	
  }
  mod1 <- glm(Cases~as.factor(Age)+Period +as.factor(Cohort) + offset(log(y)) -1,family=poisson,data=apcdata)
  mod2 <- glm(Cases~as.factor(Age)+Period +I(Period^2)+as.factor(Cohort) + offset(log(y)) -1,family=poisson,data=apcdata)
  if (Rplatform) {
    pdiff<-anova(mod1,mod2,test="Chisq")$"P(>|Chi|)"[2]
    # Correction added 2015-09-19 for compabilty with newer R versions:
    if (is.null(pdiff))
      pdiff <- anova(mod1,mod2,test="Chisq")$"Pr(>Chi)"[2]  
  } else {
    pdiff<-anova(mod1,mod2,test="Chisq")[2,7] 	
  } 
  if (pdiff <0.05) {
    suggestionrecent <-T
  } else {
    suggestionrecent <-F
  }   
  res <- list(glm=res.glm,cases=cases,pyr=pyr,noperiod=noperiod,gofpvalue=pvalue,startestage=startestage,
              suggestionrecent=suggestionrecent,pvaluerecent=pdiff,linkfunc=linkfunc)
  class(res) <- "nordpred.estimate"
  attr(res,"Call") <- sys.call()
  return(res)
}

nordpred.prediction <- function(nordpred.estimate.object,startuseage,recent,cuttrend=c(0,.25,.5,.75,.75)) {
  if (class(nordpred.estimate.object)!="nordpred.estimate") {
    stop("Variable \"nordpred.estimate.object\" must be of type \"nordpred.estimate\"")	
  } 
  if (nordpred.estimate.object$startestage>startuseage) {
    stop("\"startuseage\" is set to high compared to \"startestage\" in \"nordpred.estimate.object\"") 	
  }
  if (length(cuttrend)<(dim(nordpred.estimate.object$pyr)[2]-dim(nordpred.estimate.object$cases)[2])) {
    err <- paste("\"cuttrend\" must always be at least the same length as")
    err <- paste(err,"the number of periods with population forecasts")
    stop(err)
  }
  cases     <-  nordpred.estimate.object$cases
  pyr	    <-  nordpred.estimate.object$pyr
  noperiod  <- nordpred.estimate.object$noperiod
  nototper  <- dim(pyr)[2]
  noobsper  <- dim(cases)[2]
  nonewpred <- nototper-noobsper
  if (length(cuttrend)>nonewpred) {
    cuttrend <- cuttrend[1:nonewpred] 
  }
  if (is.data.frame(pyr)) {
    years <- names(pyr)
  } else {
    if (is.null(dimnames(pyr))) {
      years <- paste("Periode",1:nototper)	
    } else {
      years <- dimnames(pyr)[[2]]
    }	
  }  
  datatable <- matrix(NA,18,nototper)
  datatable[,1:(nototper-nonewpred)] <- as.matrix(cases)
  datatable <- data.frame(datatable)
  row.names(datatable) <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85+")
  names(datatable)     <- years
  for (age in 1:(startuseage-1)) {
    obsinc <- cases[age,(noobsper-1):noobsper]/pyr[age,(noobsper-1):noobsper]
    if (sum(is.na(obsinc))) {
      obsinc[is.na(obsinc)] <- 0 	
    }
    datatable[age,(noobsper+1):nototper] <- ((obsinc[,1]+obsinc[,2])/2)*pyr[age,(noobsper+1):nototper]
  }
  for (age in startuseage:18) {
    startestage  <- nordpred.estimate.object$startestage
    coefficients <- nordpred.estimate.object$glm$coefficients
    coh       <- (18-startestage) - (age-startestage) + (noperiod+1:nonewpred)
    noages    <- 18-startestage+1
    driftmp   <- cumsum(1-cuttrend)
    cohfind   <- noages + (noperiod-1) + 1 + (coh-1)
    maxcoh    <- 18 - startuseage + noperiod 
    agepar    <- as.numeric(coefficients[age-startestage+1])
    driftfind <- pmatch("Period",attributes(coefficients)$names)
    driftpar  <- as.numeric(coefficients[driftfind])
    cohpar <- rep(NA,length(coh))
    for (i in 1:length(coh)) {
      if (coh[i] < maxcoh) {
        cohpar[i] <- as.numeric(coefficients[cohfind[i]])
      } else {
        # For young cohorts not estimated:
        cohpar[i] <- as.numeric(coefficients[length(coefficients)-(startuseage-startestage)])
        cohpar[i][is.na(cohpar[i])] <- 0
      }
    } 
    if (recent) {
      lpfind <- driftfind + noperiod-2
      lppar  <-as.numeric(coefficients[lpfind])
      driftrecent <- driftpar - lppar
    }
    if (nordpred.estimate.object$linkfunc=="power5") { 
      if (recent) {
        rate <- (agepar+driftpar*noobsper+driftrecent*driftmp+cohpar)^5      
      } else {
        rate <- (agepar+driftpar*(noobsper+driftmp)+cohpar)^5
      }
    } else { # Possion link:
      if (recent) {
        rate <- exp(agepar+driftpar*noobsper+driftrecent*driftmp+cohpar)      
      } else {
        rate <- exp(agepar+driftpar*(noobsper+driftmp)+cohpar)
      }
    }
    datatable[age,(noobsper+1):nototper] <- rate*pyr[age,(noobsper+1):nototper]
  }
  res <- list(predictions=datatable,pyr=pyr,nopred=nonewpred,noperiod=nordpred.estimate.object$noperiod,
              gofpvalue=nordpred.estimate.object$gofpvalue,recent=recent,pvaluerecent=nordpred.estimate.object$pvaluerecent,
              cuttrend=cuttrend,startuseage=startuseage,startestage=nordpred.estimate.object$startestage,
              glm=nordpred.estimate.object$glm)
  class(res) <- "nordpred"
  attr(res,"Call") <- sys.call()
  return(res)
}
