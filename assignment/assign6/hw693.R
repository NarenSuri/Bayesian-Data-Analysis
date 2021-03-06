# 9.3 Crime data Question
# a. Fit a regression model
library(MASS)
# in page # 159, its clearly given the algorithm to use the MontoCarlo simulation 
# using the code discussed by the professor in the class
UsCrimeData<-UScrime
# assining the values from the data to the variables or selectign the relevant columns from the US crime data
UsCrimeResultantY<-UsCrimeData$y
CrimeDataToAnalyze<-UsCrimeData[,-16]
CrimeDataToAnalyze<-cbind(rep(1,47),CrimeDataToAnalyze)
CrimeDataToAnalyze<-as.matrix(CrimeDataToAnalyze)

# loaded all the required data above
# now lets initialize and set all the parameters required to start the iterative MC process
# initialize as shown in the text book
#############################
nu0<-2
s20=1
SigmaSquare=1
SamplingIters<-10000
Beta<-S2<-NULL
##############################
xtxinv<-solve(t(CrimeDataToAnalyze)%*%CrimeDataToAnalyze)
linearModelingResult<-lm(UsCrimeResultantY~CrimeDataToAnalyze[,-1])
# using code shared in text book and by Professor
SSRg<-t(UsCrimeResultantY)%*%UsCrimeResultantY-length(UsCrimeData$y)/(length(UsCrimeData$y)+1)*t(UsCrimeResultantY)%*%predict(linearModelingResult)
for(i in 1:SamplingIters){
  s2<-1/rgamma(n=1,nu0+length(UsCrimeData$y)/2,(nu0*s20-SigmaSquare+SSRg)/2)
  ValueOfbeta<-mvrnorm(n=1,mu=length(UsCrimeData$y)/(length(UsCrimeData$y)+1)*linearModelingResult$coefficients,length(UsCrimeData$y)/(length(UsCrimeData$y)+1)*s2*xtxinv)
  Beta<-rbind(Beta,ValueOfbeta)
  S2<-rbind(S2,s2)
}

########################
## Now working on the posterior samples
posteriorsample<-data.frame(Beta,S2)
colnames(posteriorsample)<-c("Beta0","Beta1","Beta2","Beta3","Beta4","Beta5","Beta6","Beta7","Beta8","Beta9","Beta10","Beta11","Beta12","Beta13","Beta14","Beta15","Sigma2")
posteriorMean<-colMeans(posteriorsample)
quantiles<-matrix(NA,15,2)
quantiles<-apply(posteriorsample,2,function(CrimeDataToAnalyze)quantile(CrimeDataToAnalyze,c(0.025,0.975)))
print(quantiles)
print(posteriorMean)
coffecients<-linearModelingResult$coefficients
print(summary(linearModelingResult))

###################################################################################################################

#Question part 2 - b - 1
SizeOfSamlple<-0.5*length(UsCrimeData$y)
## set the seed to make your partition reproductible
set.seed(145)
trainDataindexValues <- sample(seq_len(length(UsCrimeData$y)), size = SizeOfSamlple)
train <- UsCrimeData[trainDataindexValues, ]
test <- UsCrimeData[-trainDataindexValues, ]
trainValueOfY<-as.matrix(train$y)
Xtrain<-as.matrix(train[,-16])
ytest<-as.matrix(test$y)
Xtest<-as.matrix(test[,-16])
Xtest<-cbind(rep(1,nrow(Xtest)),Xtest)
linear_model<-lm(trainValueOfY~Xtrain)
B<-linear_model$coefficients
Yhat<-Xtest%*%B
plot(Yhat,ytest)
abline(0,1)
perror<-(1/length(ytest))*sum((ytest-Yhat)^2)
print(perror)

######################################################################################

##Question part 2 - b - 2
train <- UsCrimeData[trainDataindexValues, ]
test <- UsCrimeData[-trainDataindexValues, ]
trainValueOfY<-as.matrix(train$y)
UsCrimeResultantY<-trainValueOfY
Xtrain<-as.matrix(train[,-16])
Xtrain<-cbind(rep(1,nrow(Xtrain)),Xtrain)
CrimeDataToAnalyze<-Xtrain
ytest<-as.matrix(test$y)
Xtest<-as.matrix(test[,-16])
Xtest<-cbind(rep(1,nrow(Xtest)),Xtest)

nu0<-2
s20 <- SigmaSquare<-1
SamplingIters<-10000
Beta<-S2<-NULL
n<-length(trainValueOfY)
xtxinv<-solve(t(CrimeDataToAnalyze)%*%CrimeDataToAnalyze)
linearModelingResult<-lm(UsCrimeResultantY~CrimeDataToAnalyze[,-1])
SSRg<-t(UsCrimeResultantY)%*%UsCrimeResultantY-length(trainValueOfY)/(length(trainValueOfY)+1)*t(UsCrimeResultantY)%*%predict(linearModelingResult)
for(i in 1:SamplingIters){
  s2<-1/rgamma(n=1,nu0+length(trainValueOfY)/2,(nu0*s20-SigmaSquare+SSRg)/2)
  ValueOfbeta<-mvrnorm(n=1,mu=length(trainValueOfY)/(length(trainValueOfY)+1)*linearModelingResult$coefficients,length(trainValueOfY)/(length(trainValueOfY)+1)*s2*xtxinv)
  Beta<-rbind(Beta,ValueOfbeta)
  S2<-rbind(S2,s2)
}
posteriorsample<-data.frame(Beta)
colnames(posteriorsample)<-c("Beta0","Beta1","Beta2","Beta3","Beta4","Beta5","Beta6","Beta7","Beta8","Beta9","Beta10","Beta11","Beta12","Beta13","Beta14","Beta15")
posteriorMean<-colMeans(posteriorsample)
posteriorMean<-as.matrix(posteriorMean)
Yhat<-Xtest%*%posteriorMean
plot(Yhat,ytest)
perror<-(1/length(ytest))*sum((ytest-Yhat)^2)
print(perror)
