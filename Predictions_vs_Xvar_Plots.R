#load packages and data
library(randomForestSRC)
library(MASS)
data(Boston)

#Function to calculate mode of a vector (used later)
Mode <- function(x) {
  x_table <- table(x)
  return(names(x_table)[which.max(x_table)])
}

#chas variable should be factor (Does neighborhood border Charles River?)
Boston$chas <- as.factor(Boston$chas)

#Tune
TunePar <- tune(medv~., data=Boston, nodesizeTry = c(1:9, seq(10, 100, by = 5)))

#Grow RF
RF <- rfsrc(medv~., data=Boston, nodesize = TunePar$optimal[1])


#Create new dataframe for predictions
Newdf <- Boston #start with original df

#set factor variables equal to most frequently occurring level
Newdf$chas[1:nrow(Newdf)] <- Mode(Boston$chas) #[1:nrow(Newdf)] is needed to keep all levels of factor, otherwise we get error when making predictions

#Set numeric variables equal to median value
Newdf$crim <- median(Boston$crim)
Newdf$zn <- median(Boston$zn)
Newdf$indus <- median(Boston$indus)
Newdf$nox <- median(Boston$nox)
Newdf$rm <- median(Boston$rm)
Newdf$age <- median(Boston$age)
Newdf$dis <- median(Boston$dis)
Newdf$rad <- median(Boston$rad)
Newdf$tax <- median(Boston$tax)
Newdf$ptratio <- median(Boston$ptratio)
Newdf$black <- median(Boston$black)
Newdf$lstat <- median(Boston$lstat)
Newdf$medv <- median(Boston$medv)

#Enter values for variable we're interested in studying. 
#Change to from and to arguments to appropriate values
#Don't change the length.out argument
Newdf$age <- seq(from=1, to=100, length.out=nrow(Boston))

#Make predictions
Predictions <- predict(RF, newdata=Newdf)$predicted

#Calculate prediction interval, using OOB method assuming symmetry (modify for other methods if needed)
D <- RF$predicted.oob-Boston$medv
Margin <- quantile(abs(D), 0.95)
Lower <- Predictions - Margin
Upper <- Predictions + Margin

#Create dataframe with predictions and bounds
Preddf <- data.frame(Newdf$age, Predictions, Lower, Upper)
names(Preddf)[1] <- "Age"

#Create plot
p <- ggplot(data=Preddf, aes(x=Age, y=Predictions))
p + geom_line() + geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = .3, linetype=2) + ylim(c(0,30))
  

