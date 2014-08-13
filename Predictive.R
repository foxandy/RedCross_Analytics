closeAllConnections()
rm(list=ls())
setwd("C:/Users/Andrew/Desktop/Red Cross/Model Data")
library(sqldf)
library(geosphere)
library(ggmap)
library(lubridate)
library(MASS)
library(car)
library(tree)
library(gbm)

#cross-validation
CVInd <- function(n,K) {  #n is sample size; K is number of parts; returns K-length list of indices for each part
  m<-floor(n/K)  #approximate size of each part
  r<-n-m*K  
  I<-sample(n,n)  #random reordering of the indices
  Ind<-list()  #will be list of indices for all K parts
  length(Ind)<-K
  for (k in 1:K) {
    if (k <= r) kpart <- ((m+1)*(k-1)+1):((m+1)*k)  
    else kpart<-((m+1)*r+m*(k-r-1)+1):((m+1)*r+m*(k-r))
    Ind[[k]] <- I[kpart]  #indices for kth part of data
  }
  Ind
}

###DATA CLEANING###
calls <- read.csv("callout.csv", header=T, na.strings=c("","*","NA"))
volunteers <- read.csv("volunteer_locations.csv", header=T, na.strings=c("","*","NA"))
zip_stats <- read.csv("zips.csv", header=T, na.string=c("","*","NA"))
zip_stats <- zip_stats[,1:5]
#convert income and population to Ks
zip_stats$income <- zip_stats$income / 1000
zip_stats$population <- zip_stats$population / 1000
#constant for Meters to Miles distance conversion
kMetersToMiles <- 0.000621371

#incident geocoding only needs to be run once
incidents <- sqldf("SELECT DISTINCT drNum, incidentAddress, incidentCity, incidentState FROM calls")
#full address is easier to look up
incidents$fullAddress <- paste(incidents$incidentAddress, incidents$incidentCity, incidents$incidentState)
#look up incident coodinates via the Google Maps API
incidentCoordinates <- geocode(incidents$fullAddress, output = 'more')
incidents <- cbind(incidents, incidentCoordinates$lon, incidentCoordinates$lat, incidentCoordinates$postal_code)
incidents <- incidents[,c(1, 5:8)]
colnames(incidents)[3] = "incidentLong"
colnames(incidents)[4] = "incidentLat"
colnames(incidents)[5] = "incidentZip"
incidents <- incidents[which(substring(incidents$incidentZip, 1, 1) == '4' |
                             substring(incidents$incidentZip, 1, 1) == '6'),]

#write and read incident data since geocoding takes a while
write.csv(incidents, "incidents.csv")

incidents <- read.csv("incidents.csv", header=T, na.strings=c("","*","NA"))

#join all data into a regression data set
calls_interim <- sqldf("SELECT C.*, I.incidentLong, I.incidentLat, I.incidentZip, V.volLat, V.volLong
                      FROM calls AS C
                      LEFT OUTER JOIN incidents AS I
                        ON C.drNum = I.drNum
                      LEFT OUTER JOIN volunteers AS V
                        ON C.volName = V.volName")

calls_joined <- sqldf("SELECT C.*, Z.income, Z.population
                      FROM calls_interim AS C
                      LEFT OUTER JOIN zip_stats AS Z
                        ON C.incidentZip = Z.zip")

calls_joined$responseBinary <- ifelse(calls_joined$response == 'Responding'
                                      | calls_joined$response == 'No Longer Needed',
                                      1, 0)
calls_joined$response <- calls_joined$responseBinary
calls_joined$roleTrainee <- ifelse(calls_joined$volType == 'DAT-Trainee', 1, 0)
calls_joined$roleFull <- ifelse(calls_joined$volType == 'DAT-Full Responder', 1, 0)
calls_joined$schedule <- ifelse(calls_joined$onSchedule == TRUE, 1, 0)
#Extract hour from timestamp using lubridate package
calls_joined$hour <- hour(mdy_hm(calls_joined$incidentTime))
calls_joined$lateNight <- ifelse(calls_joined$hour < 6, 1, 0)
calls_joined$morning <- ifelse(calls_joined$hour >= 6 & calls_joined$hour < 12, 1, 0)
calls_joined$afternoon <- ifelse(calls_joined$hour >= 12 & calls_joined$hour < 18, 1, 0)
calls_joined$evening <- ifelse(calls_joined$hour > 18, 1, 0)
calls_joined$downtown <- ifelse(toupper(calls_joined$incidentCity) == 'CHICAGO', 1, 0)
#day of week - 1 is Sunday
calls_joined$dayOfWeek <- wday(mdy_hm(calls_joined$incidentDate))
calls_joined$weekend <- ifelse(calls_joined$dayOfWeek > 5 | calls_joined$dayOfWeek == 1, 1, 0)
calls_joined$distance <- distHaversine(cbind(calls_joined$incidentLong, calls_joined$incidentLat),
                                       cbind(calls_joined$volLong, calls_joined$volLat)) * kMetersToMiles
calls_joined$distance <- ifelse(calls_joined$distance > 500, NA, calls_joined$distance)

calls_joined_2013 <- calls_joined[which(mdy_hm(calls_joined$incidentDate) < mdy_hm("7/1/2013 00:00")),]
reg_variables <- c("response", "schedule", "roleTrainee", "roleFull",
                   "lateNight", "afternoon", "evening",
                   "downtown", "weekend", "distance", "income", "population")
dispatch <- calls_joined_2013[reg_variables]

#regression tiiiiime
dispatch <- dispatch[complete.cases(dispatch),]
dispatch$response<-factor(dispatch$response,levels=c(0,1))

#logistic - residual dev: 5920; AIC: 5937
logistic<-glm(response~.,data=dispatch,family=binomial())
logisticStep=stepAIC(logistic,direction="both")
summary(logisticStep)
vif(logisticStep)
sqrt(vif(logisticStep))

#logistic with some interactions
logistic_int <- glm(response~. + distance*lateNight + distance*evening
                    + distance*afternoon + distance*schedule + distance*roleFull
                    + distance*roleTrainee + distance*weekend, data=dispatch,
                    family=binomial())
summary(logistic_int)
logistic_int_step <- stepAIC(logistic_int, direction="both")
summary(logistic_int_step)
#res. dev: 5876, AIC: 5906
vif(logistic_int_step)

Ind<-CVInd(n=nrow(dispatch),10); K<-length(Ind)
y<-dispatch[[1]]; yhat<-y
for (k in 1:K) {
  out<-glm(response ~ schedule + roleTrainee + lateNight
           + afternoon + evening + downtown + weekend + distance
           + income + distance*lateNight + distance*evening
           + distance*roleFull + distance*schedule
           ,data=dispatch[-Ind[[k]],],family=binomial())
  phat<-as.numeric(predict(out,dispatch[Ind[[k]],])) 
  yhat[Ind[[k]]]<-as.numeric(phat >= 0.5)
}
sum(y != yhat)/length(y)  #CV misclassification rate

#CV for simple stepwise logistic: 37.0%
out<-glm(response ~ schedule + roleTrainee + roleFull
         + afternoon + downtown + weekend + income + distance

#CV for interactions: 36.5%
out<-glm(response ~ schedule + roleTrainee + lateNight
        + afternoon + evening + downtown + weekend + distance
        + income + distance*lateNight + distance*evening
        + distance*roleFull + distance*schedule       

#decsion tree
control = tree.control(nobs=nrow(dispatch), mincut = 5, minsize = 10, mindev = 0.0005)
dispatch.tree.prep<-tree(response ~ .,dispatch,control=control)
deviance(dispatch.tree.prep,detail=FALSE)
dispatch.tree<-prune.tree(dispatch.tree.prep)
plot(dispatch.tree)
plot(cv.tree(dispatch.tree.prep, , prune.tree))

#after pruning, the best regression tree has 8-10 nodes
dispatch.tree<-prune.tree(dispatch.tree.prep, best=8)
dispatch.tree

plot(dispatch.tree,type="u");text(dispatch.tree,digits=3)

#boosted tree
dispatch.boost <- gbm(response==1 ~ .,data=dispatch,
                      var.monotone=NULL,distribution="bernoulli",
                      n.trees=1000,shrinkage=0.05,interaction.depth=3,
                      bag.fraction=0.5,train.fraction=1,
                      n.minobsinnode = 10, cv.folds = 10,
                      keep.data=TRUE, verbose=FALSE)
                      
best.iter <- gbm.perf(dispatch.boost,method="cv")
best.iter
dispatch.boost$cv.error[best.iter]
summary(dispatch.boost,n.trees=best.iter)
dispatch.boost$var.names

Ind<-CVInd(n=nrow(dispatch),10); K<-length(Ind)
y<-dispatch[[1]]; yhat<-y
for (k in 1:K) {
  phat<-as.numeric(predict(dispatch.boost,newdata=dispatch[Ind[[k]],],
                           n.trees=best.iter,type="response",single.tree=FALSE))
  yhat[Ind[[k]]]<-as.numeric(phat >= 0.5)
}
sum(y != yhat)/length(y)  #CV misclassification rate


#out of time test, misclass = 36.3%
calls_joined_2014 <- calls_joined[which(mdy_hm(calls_joined$incidentDate) > mdy_hm("7/1/2013 00:00")),]
dispatch_oot <- calls_joined_2014[reg_variables]
dispatch_oot <- dispatch_oot[complete.cases(dispatch_oot),]
dispatch_oot$response<-factor(dispatch_oot$response,levels=c(0,1))
yhat<-predict(dispatch.boost,newdata=dispatch_oot,n.trees=best.iter,type="response",single.tree=FALSE)
y=dispatch_oot$response
yhat=round(yhat)
sum(y != yhat)/length(y)

#out of time test for best logistic, misclass = 39.0%
yhat <- predict(logistic_int_step, newdata=dispatch_oot,
                type = "response")
y <- dispatch_oot$response
yhat = round(yhat)
sum(y != yhat)/length(y)

hist(dispatch$distance[which(dispatch$response==1)], breaks=100, xlim=c(0,50))
hist(dispatch$distance[which(dispatch$response==0)], breaks=100, xlim=c(0,50))
