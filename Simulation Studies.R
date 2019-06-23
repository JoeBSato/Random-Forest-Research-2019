# Simulation studies to investigate RF and LM prediction intervals

###################################################################
# Simulation 1: Linear Model (one variable)

set.seed(06232019)#set seed
n <- 2000  #number of observations (trainig and test set combined)
sig <- 5  #error standard deviation
x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
e <- rnorm(n, 0, sig) #generate error term 
y <- 2*x+3 + e  #generate response values 
data <- data.frame(x, y) #combine explanatory and response varables into data frame
train <- data[1:(n/2), ] #use first half of observations as training data
test <- data[(n/2+1):n, ] #use second half of observations as test data

###################################################################
# Simulation 2: Nonlinear Model (one variable)

set.seed(06232019) #set seed
n <- 2000  #number of observations (trainig and test set combined)
sig <- 5  #error standard deviation
x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
e <- rnorm(n, 0, sig) #generate error term 
y <- 0.1*(x-7)^2-3*cos(x) + 5*log(abs(x)) + 3 + e  #generate response values 
data <- data.frame(x, y) #combine explanatory and response varables into data frame
train <- data[1:(n/2), ] #use first half of observations as training data
test <- data[(n/2+1):n, ] #use second half of observations as test data

###################################################################
# Simulation 3: Linear Model (multivariate)

set.seed(06232019) #set seed
n <- 2000  #number of observations (trainig and test set combined)
xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
sig <- 5  #error standard deviation
e <- rnorm(n, 0, sig) #generate error term 
y <- 2*X[,1]+3*X[,4]+4*X[,6]-3*X[,7]+ X[,9] + e  #generate response values 
data <- data.frame(X, y) #combine explanatory and response varables into data frame
train <- data[1:(n/2), ] #use first half of observations as training data
test <- data[(n/2+1):n, ] #use second half of observations as test data

###################################################################
# Simulation 4: Nonlinear Model (multivariate)

set.seed(06232019) #set seed
n <- 2000  #number of observations (trainig and test set combined)
xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
sig <- 5  #error standard deviation
e <- rnorm(n, 0, sig) #generate error term 
y <- (X[,1]-6)^2 + 12*cos(X[,3]) + (X[,7]-5)*(X[,8]-3) + 0.02*(X[,10]-5)^5+ e  #generate response values 
data <- data.frame(X, y) #combine explanatory and response varables into data frame
train <- data[1:(n/2), ] #use first half of observations as training data
test <- data[(n/2+1):n, ] #use second half of observations as test data
