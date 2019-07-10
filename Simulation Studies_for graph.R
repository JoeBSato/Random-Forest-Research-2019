library(randomForestSRC)
##### CR: computes the coverage rate of the model under the given data #####
CR <- function(quantile, testdata){
  attempt <- quantile[,1] <= testdata$y  & testdata$y <= quantile[,2]
  s=0
  for (k in 1:(length(attempt))) {
    if (attempt[k] == "TRUE"){s <- s+1}
  }
  return(s/length(attempt))
}


##### Results: given data, outputs MSE, Avg.Width, Coverage Rate for LR & QRF & OOB & OOB:sym #####
Results <- function(traindata, testdata, nodesize){
  Results <- matrix(NA, nrow=4,ncol=3)
  #########################
  #Linear Regression Model#
  #########################
  
  #Fit a linear model with all predictors
  m <- lm(data=traindata,y~.)
  
  #Predction intervals by m for each simulated value of y
  PI <- data.frame(predict(m, newdata=testdata, interval="prediction", level=0.95))
  
  #Count how many values of "testdata" lie in the prediction intervals by PI1
  attempt <- PI$lwr <= testdata$y & testdata$y <= PI$upr
  s=0
  for (k in 1:(length(attempt))) {
    if (attempt[k] == "TRUE"){s <- s+1}
  }
  
  # MSE
  Results[1,1] <- mean((PI$fit - testdata$y)^2)
  # Avg. width of the prediction intervals
  Results[1,2] <- mean(PI$upr-PI$lwr)#sum(PI[,3]-PI[,2])/(dim(PI)[[1]])
  #Coverage Rate of the prediction intervals by LR
  Results[1,3] <- s/(length(attempt)) #CR(quantile=PI, testdata=testdata)
  
  ################
  #Random Forests#
  ################
  library(randomForestSRC)
  
  # Grow Forests on the training dataset "trainOne"
  ranfor <- quantreg(y~., data = traindata,  splitrule="mse", nodesize=nodesize)
  # Apply the grown model on the test set "test1"
  rfTest <- quantreg(object=ranfor, newdata=testdata)
  
  # 95% Prediction intervals for the test set
  Pinter <- get.quantile(rfTest, c(.025,.975))
  
  # MSE
  Results[2,1] <- mean((testdata$y - rfTest$predicted)^2)
  # Avg. width of the prediction intervals 
  Results[2,2] <- mean(Pinter[,2]-Pinter[,1])
  # Counts of the values of testdata that lie in the prediction intervals by RF
  Results[2,3] <- CR(quantile=Pinter, testdata=testdata)#z/(length(check))
  
  #####
  #OOB#
  #####
  
  # Error in each prediction of the training data
  oDiff <- ranfor$predicted.oob - traindata$y
  # 95% OOB PI
  o <- quantile(oDiff, c(0.025,0.975)) 
  # 95% OOB PI for each prediction 
  oobPI <- matrix(NA, nrow=length(oDiff), ncol=2)
  for (k in 1:length(oDiff)) {
    for(h in 1:2){
      oobPI[k,h] <- rfTest$predicted[k] + o[h]
    }
  }
  
  #MSE (same as QRF's)
  Results[3,1] <- Results[2,1]
  #Avg. PI Width
  Results[3,2] <- mean(oobPI[,2]-oobPI[,1])
  #Coverage Rate
  Results[3,3] <- CR(quantile=oobPI, testdata=testdata)
  
  ###############
  #OOB:symmetric#
  ###############
  
  # 95% OOB:sym PI
  oSym <- quantile(abs(oDiff), 0.95) 
  # 95% OOB:sym PI for each prediction 
  oobPIsym <- matrix(NA, nrow=length(oDiff), ncol=2)
  for (k in 1:length(oDiff)) {
    oobPIsym[k,1] <- rfTest$predicted[k] - oSym
    oobPIsym[k,2] <- rfTest$predicted[k] + oSym
  }
  
  #MSE (same as QRF & OOB's)
  Results[4,1] <- Results[3,1]
  #Avg.PI Width
  Results[4,2] <- mean(oobPIsym[,2]-oobPIsym[,1])
  #Coverage Rate
  Results[4,3] <- CR(quantile = oobPIsym, testdata=testdata)
  
  ########
  #output#
  ########
  return(Results)
}


##### "mTabs": gets Avg.(MSE, Avg.width of PI-s, CR) of the dataset simulated from the given equation #####
mTabs <- function(n, np, sig, equation, loops, ns){
  # np: number of predictors = # of columns in the predictor matrix
  
  #array to store MSE & Avg.Width & CR
  opt <- array(NA, dim=c(4,3, loops, length(ns)), 
               dimnames=list(c("LR","RF","OOB","OOB:sym"),c("MSE","avg.width","CR"),c(1:loops),ns))
  
  #Make a For Loop # Maybe need to change the order of the loop/where the data is generated
  for(t in 1:loops){
    #Make datasets from specified equation
    train <- makedata(n, np, sig, equation)$train
    test  <- makedata(n, np, sig, equation)$test
    for(k in 1:length(ns)){
    #for(t in 1:loops){
      #Make datasets from specified equation
      #train <- makedata(n, np, sig, equation)$train
      #test  <- makedata(n, np, sig, equation)$test
      #Get MSE & Avg.Wid & CR for each iteration
      opt[,,t,k] <- Results(traindata=train, testdata=test, nodesize = ns[k])
    }
  }
  #Average over iterations
  optTab <- array(NA, dim=c(4, 3, length(ns)), 
                   dimnames = list(c("LR","RF","OOB","OOB:sym"),c("Avg.MSE","Avg.Wid","Avg.CR"), ns))
  for (k in 1:length(ns)){
    for (i in 1:4) {
      for (j in 1:3) {
        optTab[i,j,k] <- mean(opt[i,j,,k])
      }
    }
  }
  #output
  return(list("results"=opt, "Avg.results"=optTab))
}

##### Define the equations ##### 

# Define functions to simulate response vector y 
# given predictor vector x and error vector e 

# eq1 
eq1 <- function(variable, error){
  y <- c(rep(NA, dim(variable)[[1]])) 
  
  for (k in 1:(dim(variable)[[1]])) {
    y[k] <- 2*variable[k,1] + 3 + error[k]
  }
  
  return(y)
}

# eq2
eq2 <- function(variable, error){
  y <- c(rep(NA, dim(variable)[[1]])) 
  
  for (k in 1:dim(variable)[[1]]) {
    y[k] <- 0.1*(variable[k,1]-7)^2-3*cos(variable[k,1]) + 5*log(abs(variable[k,1])) + 3 + error[k]
  }
  return(y)
}

# eq3 
eq3 <- function(variable, error){
  
  y <- c(rep(NA,dim(variable)[[1]]))
  
  for (k in 1:(dim(variable)[[1]])) {
    y[k] <- 2*variable[k,1]+3*variable[k,4]+4*variable[k,6]-3*variable[k,7]+ variable[k,9]   + error[k]
  }
  return(y)
  
}

# eq4 
eq4 <- function(variable, error){
  
  y <- c(rep(NA,dim(variable)[[1]]))
  
  for (k in 1:(dim(variable)[[1]])) {
    y[k] <- (variable[k,1]-6)^2 + 12*cos(variable[k,3]) + (variable[k,7]-5)*(variable[k,8]-3) + 0.02*(variable[k,10]-5)^5 + error[k]
  }
  return(y)
  
}

##### Function to generate a corresponding dataset #####
makedata <- function(n, np, sig, equation){
  xvec <- runif((n)*(np), 0, 10) #generate observations for 10 predictors
  X <- matrix(xvec, ncol=np, byrow=TRUE) #arrange predictors into matrix
  e <- rnorm(n , 0, sig) #generate error term 
  y <- equation(variable=X, error=e) #generate response values 
  
  data <- data.frame(X, y) #combine predictors and responses into data frame
  trainOne <- data[1:(n/2), ] #1st half of observations as training data
  test <- data[(n/2+1):n, ] #2nd half of observations as test data
  
  return(list("train"=trainOne, "test"=test))
}

##### Find optimal nodesizes, make the vector of nodesizes #####

### Function to find optimal nodesizes and their mean ###
optimize <- function(loops, equation, n, np){
  ns <- c()
  
  for (k in 1:loops) {
    ns[k] <- tune(y~., makedata(n=n, np=np, sig=5, equation=equation)$train)$optimal[[1]]
  }
  
  return(list("Avg.NS"=mean(ns), "NS"=ns))
}

### Simulation 1 ###
refer1 <- optimize(loops=3, equation=eq1, n=2000, np=1)
# 11 nodesizes (so mean(ns)=20 is centered)
ns1 <- c( 1, 3, 5,10,15,
         20,
         25,30,35,40,45)

### Simulation 2 ### 
refer2 <- optimize(loops=3, equation=eq2, n=2000, np=1)
ns2 <- c( 1, 3, 5, 7, 9,
         10,
         11,13,15,17,19)

### Simulation 3 ###
refer3 <- optimize(loops=3, equation=eq3, n=2000, np=10)
ns3 <- c(1:11) #optimal ns were 3,2,1 respectively

### Simulation 4 ### 
refer4 <- optimize(loops=3, equation=eq4, n=2000, np=10)
ns4 <- c(1:11) #optimal ns were 3,1,2 respectively

##### Compute the values using the 4 methods, with the given nodesizes #####
# Then summarize the results into a data frame

### Function that outputs the desired dataframe 
tada <- function(ns, data){
  gyoretu <- matrix(NA, nrow=(4*length(ns)), ncol=5, 
                    dimnames=list( c(  c(rep("LR", length(ns))), c(rep("QRF",    length(ns))),
                                       c(rep("OOB",length(ns))), c(rep("OOB:sym",length(ns)))  ),
                                   c("Nodesize", "Avg.MSPE","Avg.Wid","Avg.CR","Method")
                    ))
  
  #Fill the first column, "nodesize"
  gyoretu[,1]   <- c(rep(ns,4))
  
  #2nd to 4th column
  for (k in 1:4) {
    gyoretu[((k-1)*length(ns)+1):(k*length(ns)),2:4] <- t(data$Avg.results[k,,])
    
  }

  #Put it into a data frame 
  df <- data.frame(gyoretu)
  
  #5th column
  for (k in 1:4) {
   for (t in 1:length(ns)) {
    df[((k-1)*length(ns)+t),5] <- list("LR","QRF","OOB","OOB:sym")[[k]]
  }
  }
  
  #output
  return(df)
}

### Sim 1 ###
set.seed(1)
data1 <- mTabs(n=2000, np=1,  sig=5, equation=eq1, loops=10, ns=ns1)
df1 <- tada(ns=ns1, data=data1) 

### Sim 2 ###
set.seed(2)
data2 <- mTabs(n=2000, np=1,  sig=5, equation=eq2, loops=10, ns=ns2)
df2 <- tada(ns=ns2, data=data2)

### Sim 3 ###
set.seed(3)
data3 <- mTabs(n=2000, np=10, sig=5, equation=eq3, loops=10, ns=ns3)
df3 <- tada(ns=ns3, data=data3)

### Sim 4 ###
set.seed(4)
data4 <- mTabs(n=2000, np=10, sig=5, equation=eq4, loops=10, ns=ns4)
df4 <- tada(ns=ns4, data=data4)

##### Plot #####
library(ggformula)

### Sim 1 ### 
plot1 <- 
  list("Avg.MSPE" = ggplot(data=df1, aes(x=Nodesize, y=Avg.MSPE,color=Method)) + geom_point() + geom_line(),
        "Avg.Wid" = ggplot(data=df1, aes(x=Nodesize, y=Avg.Wid, color=Method)) + geom_point() + geom_line(),
         "Avg.CR" = ggplot(data=df1, aes(x=Nodesize, y=Avg.CR, color=Method)) + geom_point() + geom_line() 
          +geom_hline(yintercept=0.95) )

### Sim 2 ### 
plot2 <- list("Avg.MSPE" = ggplot(data=df2, aes(x=Nodesize, y=Avg.MSPE,color=Method)) + geom_point() + geom_line(),
              "Avg.Wid" = ggplot(data=df2, aes(x=Nodesize, y=Avg.Wid, color=Method)) + geom_point() + geom_line(),
              "Avg.CR" = ggplot(data=df2, aes(x=Nodesize, y=Avg.CR, color=Method)) + geom_point() + geom_line() 
              +geom_hline(yintercept=0.95))

### Sim 3 ### 
plot3 <- list("Avg.MSPE" = ggplot(data=df3, aes(x=Nodesize, y=Avg.MSPE,color=Method)) + geom_point() + geom_line(),
              "Avg.Wid" = ggplot(data=df3, aes(x=Nodesize, y=Avg.Wid, color=Method)) + geom_point() + geom_line(),
              "Avg.CR" = ggplot(data=df3, aes(x=Nodesize, y=Avg.CR, color=Method)) + geom_point() + geom_line() 
              +geom_hline(yintercept=0.95))

### Sim 4 ### 
plot4 <- list("Avg.MSPE" = ggplot(data=df4, aes(x=Nodesize, y=Avg.MSPE,color=Method)) + geom_point() + geom_line(),
              "Avg.Wid" = ggplot(data=df4, aes(x=Nodesize, y=Avg.Wid, color=Method)) + geom_point() + geom_line(),
              "Avg.CR" = ggplot(data=df4, aes(x=Nodesize, y=Avg.CR, color=Method)) + geom_point() + geom_line() 
              +geom_hline(yintercept=0.95))






