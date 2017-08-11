ce= read.csv("D:/sem3/Baysien/assignment/assign4/data.csv", header= T,  sep=",")
typeof(ce)
x <- ce$glucose
h<-hist(x, xlab="Glucose", main="Histogram with Normal Curve") 
xfit<-seq(min(x),max(x),0.1) 
yfit<-dnorm(xfit,mean=mean(x),sd=sd(x)) 
yfit <- yfit*diff(h$mids[1:2])*length(x) 
lines(xfit, yfit, col="blue", lwd=2)
## Seems like normal, but the data is not completely normal.
qqnorm(ce$glucose)
qqline(ce$glucose, col = 2)
#####################################################################
