library(data.table)
library(survival)
library(smoothHR)

extractDOBFuntion <- function(inputFrame) {
  inputFrame$charID<-as.character(inputFrame$numericCHI)
  inputFrame$charID<-ifelse(nchar(inputFrame$charID)==9,paste("0",inputFrame$charID,sep=""),inputFrame$charID)
  
  yrvals <- substr(inputFrame$charID,5,6)
  yrvalsDF <- as.data.frame(yrvals);colnames(yrvalsDF) <- c("yr")
  
  yrvalsDF$yr.as.char <- as.character(yrvalsDF$yr)
  yrvalsDF$yr.as.num  <- as.numeric(yrvalsDF$yr.as.char)
  yrvalsDF$century    <- ifelse (((yrvalsDF$yr.as.num>16)&(yrvalsDF$yr.as.num<=99)),19,20)
  yrvalsDF$Y          <- paste(yrvalsDF$century,yrvalsDF$yr.as.char,sep="",collapse = NULL)
  yrvalsDF$Y.as.num   <- suppressWarnings(as.numeric(yrvalsDF$Y))
  
  datevalsdaymonth      <- substr(inputFrame$charID,1,4)
  datevals              <- paste(datevalsdaymonth,yrvalsDF$Y,sep='')
  inputFrame$dob        <- as.POSIXct(datevals,format="%d%m%Y")
  inputFrame$dob.as.num <- as.numeric(inputFrame$dob)
  
  return(inputFrame)
}

returnUnixDateTime<-function(date) {
  returnVal<-as.numeric(as.POSIXct(date, format="%Y-%m-%d", tz="GMT"))
  return(returnVal)
}

# not quite working yet
replaceCRNwithCHI<-function(inputCRN,chiCheckFrame) {
  #  inputCRN<-CRN_cpepDataset$PATnq[10]
  inputCRN<-as.character(inputCRN)
  print(inputCRN)
  chiNumber<-subset(chiCheckFrame,HPIc==inputCRN)$CHINUMn
  
  chiNumberReturn<-ifelse(length(chiNumber)>0,chiNumber,0)
  
  return(chiNumberReturn)
}

takeFirstValueFromEachPatient<-function(dateplustime1) {
  reportFrame<-data.frame(dateplustime1); colnames(reportFrame)<-c("dateplustime1")
  minDateplustime1<-min(dateplustime1)
  reportFrame$flag<-0
  reportFrame$flag<-ifelse(dateplustime1==minDateplustime1,1,0)
  
  return(reportFrame$flag)
}

simpleSurvivalPlot<-function(inputFrame,endDateUnix,startAnalysisDays,followUpDays,ylimMin) {
  
  # simpleSurvivalPlot(survivalSet,returnUnixDateTime('2016-12-16'),(365.25/12)*0, 365.25*11,0.8)
  # endDateUnix<-returnUnixDateTime('2016-12-16'); startAnalysisDays<-0; followUpDays<-365.25*10
  
  SurvivalData<-inputFrame
  SurvivalData$isDead<-ifelse(SurvivalData$DeathDateUnix>0,1,0)
  
  DaySeconds<-(60*60*24)
  shortCensorPeriodStartDay  <- DaySeconds*startAnalysisDays
  shortCensorPeriodEndDay    <- DaySeconds*followUpDays
  
  lastDOD<-endDateUnix
  SurvivalData$dateOfDischarge<-SurvivalData$dateplustime1
  SurvivalData$timeToDeath<-ifelse(SurvivalData$isDead==1,(SurvivalData$DeathDateUnix-SurvivalData$dateOfDischarge),0)
  #		SurvivalData$timeToDeath<-SurvivalData$timeToDeath/DaySeconds
  SurvivalData$timeToDeathInterval<-ifelse(SurvivalData$isDead==0,(lastDOD-SurvivalData$dateOfDischarge),SurvivalData$timeToDeath)
  SurvivalData$timeToDeathInterval[is.na(SurvivalData$timeToDeathInterval)]<-0; SurvivalData<-subset(SurvivalData,timeToDeathInterval>0)
  # SurvivalData$timeToDeathInterval<-SurvivalData$timeToDeathInterval/(60*60*24*365.25)
  
  SurvivalData$shortDeathEvent <- SurvivalData$isDead
  SurvivalData$shortDeathEvent <- ifelse(SurvivalData$isDead==1 & SurvivalData$timeToDeath>=(shortCensorPeriodStartDay) & SurvivalData$timeToDeath<(shortCensorPeriodEndDay),1,0)	
  
  xlimMax<-ifelse(shortCensorPeriodEndDay<DaySeconds*10000,shortCensorPeriodEndDay,round(max(SurvivalData$timeToDeathInterval)))
  
  #  SurvivalData$sexDigit<-ifelse(nchar(SurvivalData$charID==9),as.numeric(substr(SurvivalData$charID,8,8)),as.numeric(substr(SurvivalData$charID,9,9)))
  # SurvivalData$sexNumber<-ifelse(SurvivalData$sexDigit%%2==0,1,0)
  #  SurvivalData$sex<-factor(1*(SurvivalData$sexNumber <1),levels=0:1,labels=c("F","M"))
  
  
  mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ 1, data = SurvivalData)
  shortPlotTitle <- paste("Mortality, time ",round(shortCensorPeriodStartDay)/DaySeconds," to ",round(max(SurvivalData$timeToDeathInterval))/DaySeconds," days\n n= ",nrow(SurvivalData),", threshold: ",quantile(SurvivalData$hba1cIQRinRange)[3],sep="")
  plot(mfitAge50,mark.time=T,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,round(max(SurvivalData$timeToDeathInterval))),lwd=3,ylim=c(ylimMin,1))
  
  
  mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ (numericCPEP>quantile(numericCPEP)[3]), data = SurvivalData)
  shortPlotTitle <- paste("Mortality, time ",round(shortCensorPeriodStartDay)/DaySeconds," to ",round(max(SurvivalData$timeToDeathInterval))/DaySeconds," days\n n= ",nrow(SurvivalData),", threshold: ",quantile(SurvivalData$hba1cIQRinRange)[3],sep="")
  plot(mfitAge50,mark.time=T,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,xlimMax),lwd=3,ylim=c(ylimMin,1))
  mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ ageAtTimeOfTestYears+(numericCPEP>quantile(numericCPEP)[3]), data = SurvivalData)
  pVal <- summary(mfitAge50.coxph)$coef[,5]; HR <- round(exp(coef(mfitAge50.coxph)),2)
  legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
  summarySurvfit <- summary(mfitAge50); legendNames <- row.names(summarySurvfit$table)
  legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("bottomright",legendText,cex=0.6)
  
  print(mfitAge50.coxph)
  
  
  fit <- coxph(Surv(timeToDeathInterval, shortDeathEvent)~ageAtTimeOfTestYears+pspline(numericCPEP), data=SurvivalData, x=TRUE)
  hr1 <- smoothHR(data=SurvivalData, coxfit=fit)
  plot(hr1, predictor="numericCPEP", prob=0, conf.level=0.95)
  
  ## distance from lowest HR value
  SurvivalData$distanceFromValue <- sqrt((SurvivalData$numericCPEP-0.67)^2)
  
  mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ (distanceFromValue>quantile(distanceFromValue)[3]), data = SurvivalData)
  shortPlotTitle <- paste("Mortality, time ",round(shortCensorPeriodStartDay)/DaySeconds," to ",round(max(SurvivalData$timeToDeathInterval))/DaySeconds," days\n n= ",nrow(SurvivalData),", threshold: ",quantile(SurvivalData$hba1cIQRinRange)[3],sep="")
  plot(mfitAge50,mark.time=T,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,xlimMax),lwd=3,ylim=c(ylimMin,1))
  mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ ageAtTimeOfTestYears+(distanceFromValue>quantile(distanceFromValue)[3]), data = SurvivalData)
  pVal <- summary(mfitAge50.coxph)$coef[,5]; HR <- round(exp(coef(mfitAge50.coxph)),2)
  legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
  summarySurvfit <- summary(mfitAge50); legendNames <- row.names(summarySurvfit$table)
  legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("bottomright",legendText,cex=0.6)
  
  
  
}

simpleSurvivalPlot_timeToInsulin<-function(inputFrame,endDateUnix,ylimMin,startAnalysisDays,followUpDays,cpepThreshold) {
  
  # inputFrame<-subset(mergeCPEP_insulin_firstCPEP_firstIns,timeToFirstInsulinYears>=0)
  # endDateUnix<-returnUnixDateTime("2016-12-16")
  
  SurvivalData<-inputFrame
  
  DaySeconds<-(60*60*24)
  shortCensorPeriodStartDay  <- DaySeconds*startAnalysisDays
  shortCensorPeriodEndDay    <- DaySeconds*followUpDays
  
  SurvivalData$onInsulin<-ifelse(SurvivalData$nPrescriptionsPerID>0,1,0)
  
  lastDOD<-endDateUnix
  SurvivalData$dateOfDischarge<-SurvivalData$dateplustime1.x
  SurvivalData$timeToDeath<-ifelse(SurvivalData$onInsulin==1,(SurvivalData$dateplustime1.y-SurvivalData$dateplustime1.x),0)
  # SurvivalData$timeToDeath<-SurvivalData$timeToDeath/DaySeconds
  
  SurvivalData$timeToDeathInterval<-ifelse(SurvivalData$onInsulin==0,(lastDOD-SurvivalData$dateplustime1.x),SurvivalData$timeToDeath)
  SurvivalData$timeToDeathInterval[is.na(SurvivalData$timeToDeathInterval)]<-0; SurvivalData<-subset(SurvivalData,timeToDeathInterval>0)
  # SurvivalData$timeToDeathInterval<-SurvivalData$timeToDeathInterval/(60*60*24*365.25)
  
  SurvivalData$shortDeathEvent <- SurvivalData$onInsulin
  SurvivalData$shortDeathEvent <- ifelse(SurvivalData$onInsulin==1 & SurvivalData$timeToDeath>=(shortCensorPeriodStartDay) & SurvivalData$timeToDeath<(shortCensorPeriodEndDay),1,0)
  
  xlimMax<-ifelse(shortCensorPeriodEndDay<DaySeconds*10000,shortCensorPeriodEndDay,round(max(SurvivalData$timeToDeathInterval)))
  
  #  SurvivalData$sexDigit<-ifelse(nchar(SurvivalData$charID==9),as.numeric(substr(SurvivalData$charID,8,8)),as.numeric(substr(SurvivalData$charID,9,9)))
  # SurvivalData$sexNumber<-ifelse(SurvivalData$sexDigit%%2==0,1,0)
  #  SurvivalData$sex<-factor(1*(SurvivalData$sexNumber <1),levels=0:1,labels=c("F","M"))
  
  
  mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ (numericCPEP<cpepThreshold), data = SurvivalData)
  shortPlotTitle <- paste("insulin free survival - in those not on insulin at the time of the cpeptide test, time ",round(shortCensorPeriodStartDay)/DaySeconds," to ",xlimMax/DaySeconds," days\n n= ",nrow(SurvivalData),", threshold: ",cpepThreshold,"\ncovariables: age at time of test",sep="")
  plot(mfitAge50,mark.time=T,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,xlimMax),lwd=3,ylim=c(ylimMin,1))
  mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ ageAtTimeOfTestYears+(numericCPEP<cpepThreshold), data = SurvivalData)
  pVal <- summary(mfitAge50.coxph)$coef[,5]; HR <- round(exp(coef(mfitAge50.coxph)),2)
  legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
  summarySurvfit <- summary(mfitAge50); legendNames <- row.names(summarySurvfit$table)
  legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("bottomright",legendText,cex=0.6)
  
  print(mfitAge50.coxph)
  
  
  fit <- coxph(Surv(timeToDeathInterval, shortDeathEvent)~ageAtTimeOfTestYears+pspline(numericCPEP), data=SurvivalData, x=TRUE)
  hr1 <- smoothHR(data=SurvivalData, coxfit=fit)
  # plot(hr1, predictor="numericCPEP", prob=0, conf.level=0.95)
  
  print(summary(fit))
  
  
}

simpleSurvivalPlot_timeToInsulin_BMI<-function(inputFrame,endDateUnix,ylimMin,startAnalysisDays,followUpDays,cpepThreshold) {
  
  # inputFrame<-subset(mergeCPEP_insulin_firstCPEP_firstIns,timeToFirstInsulinYears>=0)
  # endDateUnix<-returnUnixDateTime("2016-12-16")
  
  SurvivalData<-inputFrame
  
  DaySeconds<-(60*60*24)
  shortCensorPeriodStartDay  <- DaySeconds*startAnalysisDays
  shortCensorPeriodEndDay    <- DaySeconds*followUpDays
  
  SurvivalData$onInsulin<-ifelse(SurvivalData$nPrescriptionsPerID>0,1,0)
  
  lastDOD<-endDateUnix
  SurvivalData$dateOfDischarge<-SurvivalData$dateplustime1.x
  SurvivalData$timeToDeath<-ifelse(SurvivalData$onInsulin==1,(SurvivalData$dateplustime1.y-SurvivalData$dateplustime1.x),0)
  # SurvivalData$timeToDeath<-SurvivalData$timeToDeath/DaySeconds
  
  SurvivalData$timeToDeathInterval<-ifelse(SurvivalData$onInsulin==0,(lastDOD-SurvivalData$dateplustime1.x),SurvivalData$timeToDeath)
  SurvivalData$timeToDeathInterval[is.na(SurvivalData$timeToDeathInterval)]<-0; SurvivalData<-subset(SurvivalData,timeToDeathInterval>0)
  # SurvivalData$timeToDeathInterval<-SurvivalData$timeToDeathInterval/(60*60*24*365.25)
  
  SurvivalData$shortDeathEvent <- SurvivalData$onInsulin
  SurvivalData$shortDeathEvent <- ifelse(SurvivalData$onInsulin==1 & SurvivalData$timeToDeath>=(shortCensorPeriodStartDay) & SurvivalData$timeToDeath<(shortCensorPeriodEndDay),1,0)
  
  xlimMax<-ifelse(shortCensorPeriodEndDay<DaySeconds*10000,shortCensorPeriodEndDay,round(max(SurvivalData$timeToDeathInterval)))
  
  #  SurvivalData$sexDigit<-ifelse(nchar(SurvivalData$charID==9),as.numeric(substr(SurvivalData$charID,8,8)),as.numeric(substr(SurvivalData$charID,9,9)))
  # SurvivalData$sexNumber<-ifelse(SurvivalData$sexDigit%%2==0,1,0)
  #  SurvivalData$sex<-factor(1*(SurvivalData$sexNumber <1),levels=0:1,labels=c("F","M"))
  
  
  mfitAge50<-survfit(Surv(timeToDeathInterval, shortDeathEvent) ~ (numericCPEP<cpepThreshold), data = SurvivalData)
  shortPlotTitle <- paste("insulin free survival, time ",round(shortCensorPeriodStartDay)/DaySeconds," to ",xlimMax/DaySeconds," days\n n= ",nrow(SurvivalData),", threshold: ",cpepThreshold,"\ncovariables: age at time of test",sep="")
  plot(mfitAge50,mark.time=T,lty=1:6,conf.int=F,col=c("black","red","blue","green","orange","purple"),main=shortPlotTitle,xlim=c(shortCensorPeriodStartDay,xlimMax),lwd=3,ylim=c(ylimMin,1))
  mfitAge50.coxph<-coxph(Surv(timeToDeathInterval, shortDeathEvent) ~ ageAtTimeOfTestYears+nearestBMI+(numericCPEP<cpepThreshold), data = SurvivalData)
  pVal <- summary(mfitAge50.coxph)$coef[,5]; HR <- round(exp(coef(mfitAge50.coxph)),2)
  legendText <- paste("p = ",pVal," | HR = ",HR,sep="")
  summarySurvfit <- summary(mfitAge50); legendNames <- row.names(summarySurvfit$table)
  legend("bottomleft",c(legendNames),lty=1:6,col=c("black","red","blue","green","orange","purple"),cex=0.8); legend("bottomright",legendText,cex=0.6)
  
  print(mfitAge50.coxph)
  
  
  fit <- coxph(Surv(timeToDeathInterval, shortDeathEvent)~ageAtTimeOfTestYears+nearestBMI+pspline(numericCPEP), data=SurvivalData, x=TRUE)
  hr1 <- smoothHR(data=SurvivalData, coxfit=fit)
  # plot(hr1, predictor="numericCPEP", prob=0, conf.level=0.95)
  
  print(summary(fit))
  
  
}


cpepDataset<-read.csv("~/R/GlCoSy/SD_workingSource/cpepGGC.csv",stringsAsFactors = TRUE)
chiCheck<-read.csv("~/R/GlCoSy/SD_workingSource/CHICheck.csv")
# merge with SD data
# diagnosisDataset<-read.csv("../GlCoSy/SDsource/diagnosisDateDeathDate.txt")
diagnosisDataset<-read.csv("~/R/GlCoSy/SDsource/demogALL.txt", quote = "", 
                           row.names = NULL, 
                           stringsAsFactors = FALSE)

#### ####
## alternative attempt to merge datasets
# use criteria other than age to pull from sci diabetes dataset
diagnosisDataset_truncated<-data.frame(diagnosisDataset$LinkId,diagnosisDataset$PatId,diagnosisDataset$Forename,diagnosisDataset$Surname,diagnosisDataset$BirthDate,diagnosisDataset$DiabetesMellitusType_Mapped,diagnosisDataset$DeathDate); colnames(diagnosisDataset_truncated)<-c("LinkId", "PatId", "Forename", "Surname", "BirthDate","DiabetesMellitusType_Mapped","DeathDate")
diagnosisDataset_truncated$BirthDataUnix<-as.numeric(as.POSIXct(diagnosisDataset_truncated$BirthDate, format="%Y-%m-%d", tz="GMT"))
diagnosisDataset_truncated$DeathDateUnix<-as.numeric(as.POSIXct(diagnosisDataset_truncated$DeathDate, format="%Y-%m-%d", tz="GMT"))
diagnosisDataset_truncated$DeathDateUnix[is.na(diagnosisDataset_truncated$DeathDateUnix)]<-0

extractDOByrs<-substr(cpepDataset$DOB.Age,7,8)
DOB_Years<-ifelse(extractDOByrs>15 & extractDOByrs<=99,paste("19",extractDOByrs,sep=""),paste("20",extractDOByrs,sep=""))
DOB_Full<-paste(DOB_Years,"-",substr(cpepDataset$DOB.Age,4,5),"-",substr(cpepDataset$DOB.Age,1,2),sep="")
cpepDataset$DOB_Full<-DOB_Full
cpepDataset$DOBunix<-as.numeric(as.POSIXct(cpepDataset$DOB_Full, format="%Y-%m-%d", tz="GMT"))

cpepMerge<-merge(cpepDataset,diagnosisDataset_truncated,by.x=c("Surname", "Forename", "DOBunix"),by.y=c("Surname", "Forename", "BirthDataUnix"))

# convert CPEP to numeric, and convert <0.1 to 0
cpepMerge$CPEPcharacter<-as.character(cpepMerge$CPEP)
cpepMerge$numericCPEP<-as.numeric(levels(cpepMerge$CPEP))[cpepMerge$CPEP]
cpepMerge$numericCPEP<-ifelse(cpepMerge$CPEPcharacter=="<0.10 ",0.1,cpepMerge$numericCPEP)
cpepMerge$numericCPEP<-ifelse(cpepMerge$CPEPcharacter=="<0.03 ",0.03,cpepMerge$numericCPEP)

cpepMerge<-cpepMerge[is.na(cpepMerge$numericCPEP)==FALSE,]

cpepMerge$numericGlu<-as.numeric(levels(cpepMerge$Glucose))[cpepMerge$Glucose]
cpepMerge$cpepGluRatio<-cpepMerge$numericCPEP/cpepMerge$numericGlu

cpepMerge$LinkId<-as.numeric(levels(cpepMerge$LinkId))[cpepMerge$LinkId]


## numberic CHIs
cpepMerge$numericCHI<-as.numeric(levels(cpepMerge$PatId))[cpepMerge$PatId]

## convert to data table
cpepMergeDT<-data.table(cpepMerge)

## strip test dates
cpepMergeDT$dateplustime1<-as.numeric(as.POSIXct(cpepMergeDT$DTC, format="%d.%m.%y %H:%M", tz="GMT"))

## add age
cpepMergeDT$ageAtTimeOfTest<-cpepMergeDT$dateplustime1-cpepMergeDT$DOBunix
cpepMergeDT$ageAtTimeOfTestYears<-(cpepMergeDT$ageAtTimeOfTest)/(60*60*24*365.25)

cpepMergeDT[, c("flagFirstVal") := takeFirstValueFromEachPatient(dateplustime1) , by=.(numericCHI)]

#### ####
## survivalSet<-cpepMergeDT[flagFirstVal==1]

# only analyse patients with defined diabetes diagnoses
survivalSet<-cpepMergeDT[flagFirstVal==1 & (DiabetesMellitusType_Mapped=="Type 1 Diabetes Mellitus" | DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" | DiabetesMellitusType_Mapped=="Latent Autoimmune Diabetes of Adulthood" | DiabetesMellitusType_Mapped=="Impaired Fasting Glucose" | DiabetesMellitusType_Mapped=="Maturity Onset Diabetes of Youth" | DiabetesMellitusType_Mapped=="Impaired Glucose Tolerance" | DiabetesMellitusType_Mapped=="History of Gestational Diabetes" | DiabetesMellitusType_Mapped=="Current Gestational Diabetes" | DiabetesMellitusType_Mapped=="Secondary - Pancreatic Pathology" | DiabetesMellitusType_Mapped=="Secondary - Drug Induced")]

# remove duplicate IDs
survivalSet <- survivalSet[order(survivalSet$numericCHI),]
survivalSet$chiDiff<-1
survivalSet$chiDiff[2:nrow(survivalSet)] <- diff(survivalSet$numericCHI)

survivalSet<-survivalSet[chiDiff>0]


## survivalSet<-cpepMergeDT[flagFirstVal==1 & (DiabetesMellitusType_Mapped=="Type 1 Diabetes Mellitus" | DiabetesMellitusType_Mapped=="Latent Autoimmune Diabetes of Adulthood")]

## survivalSet<-cpepMergeDT[flagFirstVal==1 & (DiabetesMellitusType_Mapped=="Type 1 Diabetes Mellitus" | DiabetesMellitusType_Mapped=="Latent Autoimmune Diabetes of Adulthood") & ageAtTimeOfTestYears<80]

logitRegressionFunction<-function(inputFrame,timePoint) {
  
  inputFrame <- subset(inputFrame, (DeathDateUnix>0 & (DeathDateUnix-dateplustime1)<=(timePoint*(60*60*24*365.25))) | ((returnUnixDateTime('2016-12-16') - dateplustime1) > (timePoint*(60*60*24*365.25))) )
  
  dead_at_timePoint<-ifelse(inputFrame$DeathDateUnix>0 & (inputFrame$DeathDateUnix-inputFrame$dateplustime1)<=(timePoint*(60*60*24*365.25)),1,0)
  
  testing_set<-data.frame(inputFrame$numericCPEP, inputFrame$ageAtTimeOfTestYears, dead_at_timePoint)
  # testing_set[-3] = scale(testing_set[-3])
  
  logitRegression<-glm(testing_set$dead_at_timePoint ~ .,
                       family = binomial,
                       data = testing_set)
  
  return(summary(logitRegression))
}

logitRegressionFunctionWithBMI<-function(inputFrame,timePoint) {
  
  inputFrame <- subset(inputFrame, (DeathDateUnix>0 & (DeathDateUnix-dateplustime1)<=(timePoint*(60*60*24*365.25))) | ((returnUnixDateTime('2016-12-16') - dateplustime1) > (timePoint*(60*60*24*365.25))) )
  
  dead_at_timePoint<-ifelse(inputFrame$DeathDateUnix>0 & (inputFrame$DeathDateUnix-inputFrame$dateplustime1)<=(timePoint*(60*60*24*365.25)),1,0)
  
  testing_set<-data.frame(inputFrame$numericCPEP, inputFrame$ageAtTimeOfTestYears, inputFrame$bmi_cpep_merge.nearestBMI, dead_at_timePoint)
  # testing_set[-3] = scale(testing_set[-3])
  
  logitRegression<-glm(testing_set$dead_at_timePoint ~ .,
                       family = binomial,
                       data = testing_set)
  
  return(summary(logitRegression))
}


# ## feature scaling makes no difference to outcome but makes coeff difficult to interpret
#  #   logitRegressionFunctionWithFeatureScaling<-function(inputFrame,timePoint) {
#       
#       dead_at_timePoint<-ifelse(inputFrame$DeathDateUnix>0 & (inputFrame$DeathDateUnix-inputFrame$dateplustime1)<=(timePoint*(60*60*24*365.25)),1,0)
#       
#       testing_set<-data.frame(inputFrame$numericCPEP, inputFrame$ageAtTimeOfTestYears, dead_at_timePoint)
#       testing_set[-3] = scale(testing_set[-3])
#       
#       logitRegression<-glm(testing_set$dead_at_timePoint ~ .,
#                            family = binomial,
#                            data = testing_set)
#       
#       return(summary(logitRegression))
#     }

# 
logitRegressionFunction(survivalSet,3)


## 
## 
plotfilename <- paste("../GlCoSy/plots/cpep_test_output.pdf",sep="")
pdf(plotfilename, width=16, height=9)

simpleSurvivalPlot(survivalSet,returnUnixDateTime('2016-12-16'),(365.25/12)*0, 365.25*3,3) 

dev.off()
#### ####

##
boxplot(cpepMergeDT[flagFirstVal==1]$numericCPEP ~ cut(cpepMergeDT[flagFirstVal==1]$ageAtTimeOfTest/(60*60*24*365.25),breaks=seq(0,100,5)),varwidth=T,ylim=c(0,2))

boxplot(cpepMergeDT[flagFirstVal==1 & (DiabetesMellitusType_Mapped=="Type 1 Diabetes Mellitus" | DiabetesMellitusType_Mapped=="Latent Autoimmune Diabetes of Adulthood")]$numericCPEP ~ cut(cpepMergeDT[flagFirstVal==1 & (DiabetesMellitusType_Mapped=="Type 1 Diabetes Mellitus" | DiabetesMellitusType_Mapped=="Latent Autoimmune Diabetes of Adulthood")]$ageAtTimeOfTest/(60*60*24*365.25),breaks=seq(0,100,5)),varwidth=T,main=paste("type 1 first cpep in record. n=",nrow(cpepMergeDT[flagFirstVal==1 & (DiabetesMellitusType_Mapped=="Type 1 Diabetes Mellitus" | DiabetesMellitusType_Mapped=="Latent Autoimmune Diabetes of Adulthood")]),sep=""),ylim=c(0,1))

boxplot(cpepMergeDT[flagFirstVal==1 & (DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" | DiabetesMellitusType_Mapped=="Impaired Fasting Glucose" | DiabetesMellitusType_Mapped=="Maturity Onset Diabetes of Youth" | DiabetesMellitusType_Mapped=="Impaired Glucose Tolerance" | DiabetesMellitusType_Mapped=="History of Gestational Diabetes" | DiabetesMellitusType_Mapped=="Current Gestational Diabetes")]$numericCPEP ~ 
          cut(cpepMergeDT[flagFirstVal==1 & (DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" | DiabetesMellitusType_Mapped=="Impaired Fasting Glucose" | DiabetesMellitusType_Mapped=="Maturity Onset Diabetes of Youth" | DiabetesMellitusType_Mapped=="Impaired Glucose Tolerance" | DiabetesMellitusType_Mapped=="History of Gestational Diabetes" | DiabetesMellitusType_Mapped=="Current Gestational Diabetes")]$ageAtTimeOfTest/(60*60*24*365.25),breaks=seq(0,100,10)),
        varwidth=T,
        main=paste("type 2 first cpep in record. n=",nrow(cpepMergeDT[flagFirstVal==1 & (DiabetesMellitusType_Mapped=="Type 2 Diabetes Mellitus" | DiabetesMellitusType_Mapped=="Impaired Fasting Glucose" | DiabetesMellitusType_Mapped=="Maturity Onset Diabetes of Youth" | DiabetesMellitusType_Mapped=="Impaired Glucose Tolerance" | DiabetesMellitusType_Mapped=="History of Gestational Diabetes" | DiabetesMellitusType_Mapped=="Current Gestational Diabetes")]),sep=""),
        ylim=c(0,2))

cpepMergeDT_firstCPEPval<-cpepMergeDT[flagFirstVal==1]

## write out IDs for drug match
# outputName<-paste("../cPeptideData/cpepIDs2.csv",sep="")
# write.table(mergeCPepSD[flagFirstVal==1]$LinkId,file=outputName,sep=",",append=F,col.names=T)

## import drug data
cpepDrugs<-read.csv("~/R/GlCoSy/SD_workingSource/cpep_drugs_byLinkId.txt",stringsAsFactors = TRUE)
cpepDrugs$LinkId<-as.numeric(levels(cpepDrugs$LinkId))[cpepDrugs$LinkId]
cpepDrugsDT<-data.table(cpepDrugs)

cpepDrugsDT_insulin<-cpepDrugsDT[BNFCode=="6.1.1.1" | BNFCode=="6.1.1.51"| substr(BNFCode,1,7)=="0601011" | BNFCode=="6.1.1.2" | substr(BNFCode,1,7)=="0601012"]
cpepDrugsDT_insulin$dateplustime1<-as.numeric(as.POSIXct(cpepDrugsDT_insulin$PrescriptionDateTime, format="%Y-%m-%d", tz="GMT"))
cpepDrugsDT_insulin<-cpepDrugsDT_insulin[order(cpepDrugsDT_insulin$LinkId,cpepDrugsDT_insulin$dateplustime1)]

## insulin prescription timings
cpepDrugsDT_insulin[, c("timeFirstInsulin") := min(dateplustime1) , by=.(LinkId)]
cpepDrugsDT_insulin[, c("timeLastInsulin") := max(dateplustime1) , by=.(LinkId)]
cpepDrugsDT_insulin[, nPrescriptionsPerID := .N , by=.(LinkId)]
cpepDrugsDT_insulin[, nPrescriptionsInSequencePerID := seq(1,.N,1) , by=.(LinkId)]

cpepDrugsDT_insulin_firstPrescription<-cpepDrugsDT_insulin[nPrescriptionsInSequencePerID==1]

cpepMergeDT_firstCPEPval_firstIns<-merge(survivalSet,cpepDrugsDT_insulin_firstPrescription,by.x="LinkId",by.y="LinkId",all.x=T)

# cpepMergeDT_firstCPEPval_firstIns<-merge(cpepMergeDT_firstCPEPval,cpepDrugsDT_insulin_firstPrescription,by.x="LinkId",by.y="LinkId",all.x=T)

## need to add death to composite endpoint of prescription of insulin or death

mergeCPEP_insulin_firstCPEP_firstIns<-cpepMergeDT_firstCPEPval_firstIns

cpepMergeDT_firstCPEPval_firstIns$nPrescriptionsPerID[is.na(cpepMergeDT_firstCPEPval_firstIns$nPrescriptionsPerID)]<-0
cpepMergeDT_firstCPEPval_firstIns$timeToFirstInsulin<-cpepMergeDT_firstCPEPval_firstIns$timeFirstInsulin-cpepMergeDT_firstCPEPval_firstIns$dateplustime1.x
cpepMergeDT_firstCPEPval_firstIns$timeToFirstInsulinYears<-cpepMergeDT_firstCPEPval_firstIns$timeToFirstInsulin/(60*60*24*365.25)
cpepMergeDT_firstCPEPval_firstIns$timeToFirstInsulinYears[is.na(cpepMergeDT_firstCPEPval_firstIns$timeToFirstInsulinYears)]<-0

## start plots
hist(cpepMergeDT_firstCPEPval_firstIns$numericCPEP,breaks=seq(0,8,0.1))
boxplot(subset(cpepMergeDT_firstCPEPval_firstIns,timeToFirstInsulinYears>0)$timeToFirstInsulinYears ~ cut(subset(cpepMergeDT_firstCPEPval_firstIns,timeToFirstInsulinYears>0)$numericCPEP,breaks=seq(0,8,0.1)),ylim=c(0,1))

## logistic regression ? on insulin at 6m, 1y, 2y
insulin_at_0.5y<-ifelse(cpepMergeDT_firstCPEPval_firstIns$nPrescriptionsPerID>0 & cpepMergeDT_firstCPEPval_firstIns$timeToFirstInsulinYears<=0.5,1,0)
insulin_at_1y<-ifelse(cpepMergeDT_firstCPEPval_firstIns$nPrescriptionsPerID>0 & cpepMergeDT_firstCPEPval_firstIns$timeToFirstInsulinYears<=1,1,0)
insulin_at_2y<-ifelse(cpepMergeDT_firstCPEPval_firstIns$nPrescriptionsPerID>0 & cpepMergeDT_firstCPEPval_firstIns$timeToFirstInsulinYears<=1,2,0)

logitRegression<-glm(insulin_at_0.5y ~ cpepMergeDT_firstCPEPval_firstIns$numericCPEP + cpepMergeDT_firstCPEPval_firstIns$ageAtTimeOfTestYears); summary(logitRegression)
logitRegression<-glm(insulin_at_1y ~ cpepMergeDT_firstCPEPval_firstIns$numericCPEP + cpepMergeDT_firstCPEPval_firstIns$ageAtTimeOfTestYears); summary(logitRegression)
logitRegression<-glm(insulin_at_2y ~ cpepMergeDT_firstCPEPval_firstIns$numericCPEP + cpepMergeDT_firstCPEPval_firstIns$ageAtTimeOfTestYears); summary(logitRegression)

# subset to insulin naive
cpepMergeDT_firstCPEPval_firstIns <- subset(cpepMergeDT_firstCPEPval_firstIns,timeToFirstInsulinYears>=0)
    # logisitc regression loop:
logitPlotFrame <- as.data.frame(matrix(0,nrow=0, ncol=2))
colnames(logitPlotFrame) <- c("time", "estimate")

    for (i in seq(0,24,0.01)) {
      insulin_at_i_y<-ifelse(cpepMergeDT_firstCPEPval_firstIns$nPrescriptionsPerID>0 & cpepMergeDT_firstCPEPval_firstIns$timeToFirstInsulinYears>0 & cpepMergeDT_firstCPEPval_firstIns$timeToFirstInsulinYears<=(i/12),1,0)
      logitRegression<-glm(insulin_at_i_y ~ cpepMergeDT_firstCPEPval_firstIns$numericCPEP + cpepMergeDT_firstCPEPval_firstIns$ageAtTimeOfTestYears); summary(logitRegression)
      
      reportFrame <- data.frame(i, logitRegression$coefficients[2])
      colnames(reportFrame) <- c("time", "estimate")
      
      logitPlotFrame <- rbind(logitPlotFrame, reportFrame)

    }
plot(logitPlotFrame$time, logitPlotFrame$estimate)
lines(logitPlotFrame$time, logitPlotFrame$estimate)





plotfilename <- paste("../GlCoSy/plots/survivalPlot_timeToInsulin.pdf",sep="")
pdf(plotfilename, width=16, height=9)

medianCPEP<-quantile(subset(cpepMergeDT_firstCPEPval_firstIns,timeToFirstInsulinYears>=0)$numericCPEP)[3]
# insulin free survival
simpleSurvivalPlot_timeToInsulin(subset(cpepMergeDT_firstCPEPval_firstIns,timeToFirstInsulinYears>=0), returnUnixDateTime('2016-12-16'), 0, (365.25/12)*0, 365.25*1,medianCPEP)


dev.off()

logitRegressionFunction<-function(inputFrame,timePoint) {
  insulin_at_timePoint<-ifelse(inputFrame$nPrescriptionsPerID>0 & inputFrame$timeToFirstInsulinYears<=timePoint,1,0)
  logitRegression<-glm(insulin_at_timePoint ~ inputFrame$numericCPEP + inputFrame$ageAtTimeOfTestYears)
  
  return(summary(logitRegression))
}

logitRegressionFunction(subset(cpepMergeDT_firstCPEPval_firstIns,timeToFirstInsulinYears>=0),1)
logitRegressionFunction(subset(cpepMergeDT_firstCPEPval_firstIns,timeToFirstInsulinYears>=0),1)

################# ################
## start to add covariables

## first CPEP dataset
# cpepMergeDT_firstCPEPval_firstIns
# mergeCPEP_insulin_firstCPEP_firstIns   

## load in BMI dataset
bmiSetDF<-read.csv("../GlCoSy/SD_workingSource/BMISetDTclean.csv")
bmiSetDT<-data.table(bmiSetDF)
bmiSetDT$dateplustime1_bmi<-bmiSetDT$dateplustime1; bmiSetDT$dateplustime1<-NULL


## find most recent BMI prior to CPEP measurement for each ID
bmi_cpep_merge<-merge(cpepMergeDT_firstCPEPval_firstIns,bmiSetDT,by.x="LinkId",by.y="LinkId")
bmi_cpep_merge$testDate_bmiDate_difference<-bmi_cpep_merge$dateplustime1.x - bmi_cpep_merge$dateplustime1_bmi
## all BMIs performed before the CPEP or up to 90 days post CPEP
bmi_cpep_merge<-bmi_cpep_merge[testDate_bmiDate_difference>0 | testDate_bmiDate_difference>(60*60*24*365.25*10*-1)]

findNearestMeasure<-function(all_metric_dateplustime1, metricNumeric, testDate_metricDate_difference,postTestTimeLimitDays) {
  # all_metric_dateplustime1<-bmi_cpep_merge[LinkId==2147487389]$dateplustime1_bmi; metricNumeric<-bmi_cpep_merge[LinkId==2147487389]$bmiNumeric; testDate_metricDate_difference<-bmi_cpep_merge[LinkId==2147487389]$testDate_bmiDate_difference; postTestTimeLimitDays<-90

  testFrame<-data.frame(all_metric_dateplustime1, metricNumeric, testDate_metricDate_difference)
  testFrame<-subset(testFrame,testDate_metricDate_difference>(60*60*24*postTestTimeLimitDays*-1))
  testFrame$flag<-ifelse(testFrame$testDate_metricDate_difference>0,1,0)
  
  if (sum(testFrame$flag>0)) {outputVals<-subset(testFrame,testDate_metricDate_difference==min(subset(testFrame,flag==1)$testDate_metricDate_difference))}
  if (sum(testFrame$flag==0)) {outputVals<-subset(testFrame,testDate_metricDate_difference==max(subset(testFrame,flag==0)$testDate_metricDate_difference))}
  
  returnVals<-list(outputVals$all_metric_dateplustime1, outputVals$metricNumeric)
  
  
  return(returnVals)
  
}

bmi_cpep_merge[, c("nearestBMI_dateplustime1", "nearestBMI") := findNearestMeasure(dateplustime1_bmi, bmiNumeric, testDate_bmiDate_difference,90) , by=.(LinkId)]


bmi_cpep_merge[, nBMImeasuresPerID := seq(1,.N,1) , by=.(LinkId)]

bmi_cpep_merge<-bmi_cpep_merge[nBMImeasuresPerID==1]

bmi_cpep_merge <- bmi_cpep_merge[nearestBMI>0]
# write.table(bmi_cpep_merge, file = "../GlCoSy/SD_workingSource/cpep_bmi_insulinNaive.csv",sep = ",")

## plot bmi cpep threshold
plotfilename <- paste("../GlCoSy/plots/CPEP_all_withBMIdata_CPEPthresh_0.76.pdf",sep="")
pdf(plotfilename, width=16, height=9)

medianCPEP <- quantile(subset(bmi_cpep_merge,timeToFirstInsulinYears>=0)$numericCPEP)[3]
simpleSurvivalPlot_timeToInsulin_BMI(subset(bmi_cpep_merge,timeToFirstInsulinYears>=0), returnUnixDateTime('2016-12-16'), 0, (365.25/12)*0, 365.25*1,medianCPEP)


dev.off()

#######################
# logit with BMI as coV

# mortality
bmi_cpep_merge_forSurvival <- data.frame(bmi_cpep_merge$LinkId, bmi_cpep_merge$nearestBMI)
survivalSet_BMI <- merge(survivalSet, bmi_cpep_merge_forSurvival, by.x = "LinkId", by.y = "bmi_cpep_merge.LinkId")

logitRegressionFunctionWithBMI(survivalSet_BMI,3)


# time to insulin
## logistic regression ? on insulin at 6m, 1y, 2y
insulin_at_0.5y<-ifelse(bmi_cpep_merge$nPrescriptionsPerID>0 & bmi_cpep_merge$timeToFirstInsulinYears<=0.5,1,0)
insulin_at_1y<-ifelse(bmi_cpep_merge$nPrescriptionsPerID>0 & bmi_cpep_merge$timeToFirstInsulinYears<=1,1,0)
insulin_at_2y<-ifelse(bmi_cpep_merge$nPrescriptionsPerID>0 & bmi_cpep_merge$timeToFirstInsulinYears<=1,2,0)

logitRegression<-glm(insulin_at_0.5y ~ bmi_cpep_merge$numericCPEP + bmi_cpep_merge$ageAtTimeOfTestYears + bmi_cpep_merge$nearestBMI); summary(logitRegression)
logitRegression<-glm(insulin_at_1y ~ bmi_cpep_merge$numericCPEP + bmi_cpep_merge$ageAtTimeOfTestYears + bmi_cpep_merge$nearestBMI); summary(logitRegression)
logitRegression<-glm(insulin_at_2y ~ bmi_cpep_merge$numericCPEP + bmi_cpep_merge$ageAtTimeOfTestYears + bmi_cpep_merge$nearestBMI); summary(logitRegression)

#######################
# BMI subgroup analysis

cpep_BMI<-bmi_cpep_merge[bmiNumeric<30]
medianCPEP<-quantile(subset(cpep_BMI,timeToFirstInsulinYears>=0)$numericCPEP)[3]
##
simpleSurvivalPlot_timeToInsulin(subset(cpep_BMI,timeToFirstInsulinYears>=0), returnUnixDateTime('2016-12-16'), 0, (365.25/12)*0, 365.25*2,medianCPEP)

simpleSurvivalPlot_timeToInsulin_BMI(subset(cpep_BMI,timeToFirstInsulinYears>=0), returnUnixDateTime('2016-12-16'), 0, (365.25/12)*0, 365.25*2,medianCPEP)

##########################
# ppv of c-peptide value for insulin prescription at n years

ppv_function <- function(inputFrame, n_years, cpeptide_bin_size, xlimVals) {
  
  binNumber <- round((max(inputFrame$numericCPEP) - 0) / cpeptide_bin_size,0)
  reportingFrame <- data.frame(matrix(0,ncol=7, nrow=binNumber)); colnames(reportingFrame)<-c("binN", "lowerLimitCPEP", "upperLimitCPEP", "n_IDs", "n_on_insulin", "ppv", "npv")
  reportingFrame$binN <- c(1:binNumber)
  reportingFrame$lowerLimitCPEP <- seq(0,max(inputFrame$numericCPEP),cpeptide_bin_size)
  reportingFrame$upperLimitCPEP[1:(nrow(reportingFrame)-1)] <- reportingFrame$lowerLimitCPEP[2:nrow(reportingFrame)]
  reportingFrame$upperLimitCPEP[nrow(reportingFrame)] <- reportingFrame$lowerLimitCPEP[nrow(reportingFrame)] + cpeptide_bin_size
  
  for (j in seq(1,binNumber,1)) {
    sub <- subset(inputFrame, numericCPEP > reportingFrame$lowerLimitCPEP[j] & numericCPEP <= reportingFrame$upperLimitCPEP[j])
    reportingFrame$n_IDs[j] <- nrow(sub)
    
    insulin_sub <- subset(sub, nPrescriptionsPerID>0 & timeToFirstInsulinYears<=n_years)
    reportingFrame$n_on_insulin[j] <- nrow(insulin_sub)
    
  }
  
  reportingFrame$ppv <- reportingFrame$n_on_insulin / reportingFrame$n_IDs
  reportingFrame$npv <- 1-reportingFrame$ppv
 # plot(reportingFrame$lowerLimitCPEP,reportingFrame$ppv, xlim=xlimVals,pch=16, cex=(sqrt(reportingFrame$n_IDs))/2)
#  print(reportingFrame)
  
  return(reportingFrame)
  
  
}

# subset of those not on insulin at the time of c-peptide measurement    
ppv_test_set <- subset(cpepMergeDT_firstCPEPval_firstIns,timeToFirstInsulinYears>=0)

predictiveFrame<-ppv_function(ppv_test_set, 1, 0.2, c(0,2))

# plot probability over time per bin
# need to comment out plot in ppv_function function before using this loop
totalFU_years <- 1
bin <- 1/(365.25*2)
n_bins <- totalFU_years/bin
plotFrame <- as.data.frame(matrix(0, ncol=2, nrow=n_bins)); colnames(plotFrame) <- c("time", "prob")
plotFrame$time <- seq(bin, totalFU_years, bin)


for (p in seq(1,5,1)) {
  
    for (t in seq(1,n_bins,1)) {
      
      outputData<-ppv_function(ppv_test_set, plotFrame$time[t], 0.2, c(0,2))
      plotFrame$prob[t] <- outputData$ppv[p]
      
    }
    if (p==1) { plot(log(plotFrame$time), plotFrame$prob, col=p, ylim=c(0,1)); lines(log(plotFrame$time), plotFrame$prob, col=p) }
  if (p>1) { points(log(plotFrame$time), plotFrame$prob, col=p); lines(log(plotFrame$time), plotFrame$prob, col=p) }
    
}


