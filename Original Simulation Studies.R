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
mTabs <- function(n, np, sig, equation, loops, ns, error){
  # np: number of predictors = # of columns in the predictor matrix
  
  #array to store MSE & Avg.Width & CR
  opt <- array(NA, dim=c(4,3, loops, length(ns)), 
               dimnames=list(c("LR","RF","OOB","OOB:sym"),c("MSE","avg.width","CR"),c(1:loops),ns))
  #store the simulated datas
  storedata <- data.frame(matrix(NA, nrow=n, ncol=(np+1)))
  if (dim(storedata)[[2]]==2){colnames(storedata) <- c("X1","Y")}
  else       {colnames(storedata) <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10","Y")}
  
  
  #Make a For Loop 
  for(t in 1:loops){
    #Make datasets from specified equation
    preserve <- makedata(n, np, sig, equation, error)
    train <- preserve[[1]]
    test  <- preserve[[2]]
    
    if(t==1){
    storedata[(1:(n/2)),]      <- train
    storedata[(((n/2)+1):n),] <- test
    }
    
    for(k in 1:length(ns)){
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
  return(list("results"=opt, "Avg.results"=optTab, "storeddata"=storedata))
}

##### Define the error terms #####

### error1 ###
# Normal & heteroscedastic #
error1 <- function(variable,n){
  e <- rnorm(n=dim(variable)[[1]], mean=0, sd=abs(variable[,1]))
  
  #e <- c(rep(NA, dim(variable)[[1]])) 
  #for (k in 1:(dim(variable)[[1]])){ 
   # e[k] <- rnorm(n=1, mean=0, sd=abs(variable[k,1]))
  #}
  return(e)
}

### error2 ###
# Normal & homoscedastic #
error2 <- function(variable,n){
  e <- rnorm(n=dim(variable)[[1]], mean=0, sd=5)
  
   #e <- c(rep(NA, dim(variable)[[1]])) 
  #for (k in 1:(dim(variable)[[1]])){ 
   # e[k] <- rnorm(n=1, mean=0, sd=5)
  #}
  return(e)
}

### error3 ###
# non-normal & symmetric & homoscedastic #
# t-distribution with df=3 #
error3 <- function(variable,n){
  e <- rt(n=dim(variable)[[1]], df=3)
  
  #e <- c(rep(NA, dim(variable)[[1]])) 
  #for (k in 1:(dim(variable)[[1]])){ 
   # e[k] <- rt(n=1, df=3)
  #}
  return(e)
}

### error4 ###
# non-normal & symmetric & heteroscedastic #
# t-distribution + error1 # 
error4 <- function(variable,n){
  #e <- c(rep(NA, dim(variable)[[1]])) 
  #for (k in 1:(dim(variable)[[1]])){ 
   # e[k] <- rt(n=1, df=3) + rnorm(n=1, mean=0, sd=abs(variable[k,1]))
  #}
  e <- rt(n=dim(variable)[[1]], df=3) 
  + rnorm(n=dim(variable)[[1]], mean=0, sd=abs(variable[k,1]))
  
  return(e)
}

### error5 ###
# non-symmetric & homoscedastic # 
error5 <- function(variable,n){
  #e <- c(rep(NA, dim(variable)[[1]])) 
  #for (k in 1:(dim(variable)[[1]])){ 
   # e[k] <- rexp(n=1,rate=2)-0.5
  #}
  e <- rexp(n=dim(variable)[[1]],rate=2)-0.5
  
  return(e)
}

### error6 ###
# non-symmetric & heteroscedastic # 
error6 <- function(variable,n){
  #e <- c(rep(NA, dim(variable)[[1]])) 
  #for (k in 1:(dim(variable)[[1]])){ 
   # e[k] <- rexp(n=1,rate=2)-0.5 +runif(n=1, -abs(variable[k,1]), abs(variable[k,1]))
  #}
  e <- rexp(n=dim(variable)[[1]],rate=2)-0.5 
  +runif(n=dim(variable)[[1]], -abs(variable[k,1]), abs(variable[k,1]))
  
  return(e)
}

### Summarize these error functions ### 
errorAll <- list(error1,error2,error3,error4,error5,error6)

##### Define the equations ##### 

# Define functions to simulate response vector y 
# given predictor vector x and error vector e 

# eq1 
eq1 <- function(variable, error){
  #y <- c(rep(NA, dim(variable)[[1]])) 
  y <-2*variable[,1] + 3 + error
  #for (k in 1:(dim(variable)[[1]])) {
   # y[k] <- 2*variable[k,1] + 3 + error[k]
  #}
  
  return(y)
}

# eq2
eq2 <- function(variable, error){
  y <- 0.1*(variable[,1]-7)^2-3*cos(variable[,1]) + 5*log(abs(variable[,1])) + 3 + error

  #y <- c(rep(NA, dim(variable)[[1]])) 
  
  #for (k in 1:dim(variable)[[1]]) {
  #  y[k] <- 0.1*(variable[k,1]-7)^2-3*cos(variable[k,1]) + 5*log(abs(variable[k,1])) + 3 + error[k]
  #}
  return(y)
}

# eq3 
eq3 <- function(variable, error){
  #y <- c(rep(NA,dim(variable)[[1]]))
  y <- 2*variable[,1]+3*variable[,4]+4*variable[,6]-3*variable[,7]+ variable[,9] + error
  
  #for (k in 1:(dim(variable)[[1]])) {
   # y[k] <- 2*variable[k,1]+3*variable[k,4]+4*variable[k,6]-3*variable[k,7]+ variable[k,9]   + error[k]
  #}
  return(y)
}

# eq4 
eq4 <- function(variable, error){
  y <- (variable[,1]-6)^2 + 12*cos(variable[,3])  + (variable[,7]-5)*(variable[,8]-3)  
  + 0.02*(variable[,10]-5)^5 + error
  #c(rep(NA,dim(variable)[[1]]))
   
  #for (k in 1:(dim(variable)[[1]])) {
   # y[k] <- (variable[k,1]-6)^2 + 12*cos(variable[k,3]) 
  #+ (variable[k,7]-5)*(variable[k,8]-3) + 0.02*(variable[k,10]-5)^5 + error[k]
  #}
  return(y)
}

# eq5 
eq5 <- function(variable, error){
  #y <- c(rep(NA, dim(variable)[[1]])) 
  y <- 10*exp(-variable[,1])*sin(variable[,1]) + error
  
  #for (k in 1:(dim(variable)[[1]])){ 
   # y[k] <- 10*exp(-variable[k,1])*sin(variable[k,1]) + error[k]
  #}
  
  return(y)
}

### Summarize all these equations ###
eqAll <- list(eq1,eq2,eq3,eq4,eq5)

##### Function to generate a corresponding dataset #####
makedata <- function(n, np, sig, equation, error){
  xvec <- runif((n)*(np), 0, 10) #generate observations for 10 predictors
  X <- matrix(xvec, ncol=np, byrow=TRUE) #arrange predictors into matrix
  y <- equation(variable=X, error=error(variable=X, n=n)) #generate response values 
  
  data <- data.frame(X, y) #combine predictors and responses into data frame
  trainOne <- data[1:(n/2), ] #1st half of observations as training data
  test <- data[(n/2+1):n, ] #2nd half of observations as test data
  
  return(list("train"=trainOne, "test"=test))
}

##### Find optimal nodesizes, make the vector of nodesizes #####

### Function to find optimal nodesizes and their mean ###
#optimize <- function(loops, equation, n, np, error){
 # library(randomForestSRC)
  #ns <- c()
  
  #for (k in 1:loops) {
   # see <- makedata(n=n, np=np, sig=5, equation=equation, error=error)
    #ns[k]   <- tune(y~., see$train)$optimal[[1]]
    #ns[k+1] <- tune(y~., see$test )$optimal[[1]]
  #}
  
  #return(list("Avg.NS"=mean(ns), "NS"=ns))
#}

### optimal ns for each combination of base and error equation ###
nps <- c(1,1,10,10,1)
errors <- list("error1"=list(),"error2"=list(),"error3"=list(),"error4"=list(),"error5"=list(),"error6"=list())

#refer <- list("eq1"=errors,"eq2"=errors,"eq3"=errors,"eq4"=errors,"eq5"=errors)
#for (k in 1:length(eqAll)) {
 # for(t in 1:length(errorAll)){
  #  refer[[k]][[t]] <- optimize(loops=2, equation=eqAll[[k]], n=50, np=nps[k], error=errorAll[[t]])
  #}
#}

### ns vectors ###
errorsVec <- list("error1"=c(rep(NA,11)),"error2"=c(rep(NA,11)),"error3"=c(rep(NA,11)),
                  "error4"=c(rep(NA,11)),"error5"=c(rep(NA,11)),"error6"=c(rep(NA,11)))
ns <- list("eq1"=errorsVec,"eq2"=errorsVec,"eq3"=errorsVec,"eq4"=errorsVec,"eq5"=errorsVec)
for (k in 1:length(eqAll)) {
  for(t in 1:length(errorAll)){
    Flo <- floor(optimals[[k]][[t]][[2]])
    if(optimals[[k]][[t]][[2]] <= 6){
      ns[[k]][[t]] <- c(1:11)
    }
    else if(6    < optimals[[k]][[t]][[2]] & optimals[[k]][[t]][[2]] <= 6.5){
      ns[[k]][[t]] <- c(1,2,4,5,6,7,8,9,11,12,13)
    }
    else if (6.5 < optimals[[k]][[t]][[2]] & optimals[[k]][[t]][[2]] <= 7){
      ns[[k]][[t]] <- c(1,2,3,4,6,7,8,10,11,12,14)
    }
    else if (7   < optimals[[k]][[t]][[2]] & optimals[[k]][[t]][[2]] <= 7.5){
      ns[[k]][[t]] <- c(1,3,5,6,7,8,9,10,11,13,15)
    }
    else if (7.5 < optimals[[k]][[t]][[2]] & optimals[[k]][[t]][[2]] <= 8){
      ns[[k]][[t]] <- c(1,3,5,6,7,8,9,11,13,15,16)
    }
    else{
      ns[[k]][[t]][1] <- 1 
      ns[[k]][[t]][6] <- Flo 
      ns[[k]][[t]][11]<- 2*Flo 
      for (i in 1:4) {
        ns[[k]][[t]][i+1] <- floor((Flo-1)/4)*i
      }
      for (i in 1:4) {
        ns[[k]][[t]][i+6] <- Flo + floor((Flo-1)/4)*i
      }
      #"refer" and "optimals" have "nodesizes" and "nodesizesAvg" in a reverse order
    }
    
  }
}

##### Compute the values using the 4 methods, with the given nodesizes #####
# Then summarize the results into a data frame

### Function that outputs the desired dataframe 
tada <- function(ns, data){
  #matrix to be put into a dataframe
  gyoretu <- matrix(NA, nrow=(4*length(ns)), ncol=5, 
                    dimnames=list( c(  c(rep("LR", length(ns))), c(rep("QRF",    length(ns))),
                                       c(rep("OOB",length(ns))), c(rep("OOB:sym",length(ns)))  ),
                                   c("Nodesize", "Avg.MSPE","Avg.Wid","Avg.CR","Method"))         )
  
  #Fill the first column, "nodesize"
  gyoretu[,1] <- c(rep(ns,4))
  
  #2nd to 4th column
  for (k in 1:4) {
    gyoretu[((k-1)*length(ns)+1):(k*length(ns)),2:4] <- t(data[[2]][k,,])#$Avg.results
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

### 
set.seed(1)
dfs    <- list("eq1"=errors,"eq2"=errors,"eq3"=errors,"eq4"=errors,"eq5"=errors)
stored <- list("eq1"=errors,"eq2"=errors,"eq3"=errors,"eq4"=errors,"eq5"=errors)
for (k in 1:length(eqAll)) {
  for(t in 1:length(errorAll)){
    d <- mTabs(n=2000, np=nps[k],  sig=5, equation=eqAll[[k]], loops=10, ns=ns[[k]][[t]], error=errorAll[[t]])
    stored[[k]][[t]] <- d$storeddata
       dfs[[k]][[t]] <- tada(ns=ns[[k]][[t]], data=d)
  }
}

##### Plot #####

### Function for plotting ###
library(ggformula)
plotter <- function(data){
  return(list("Avg.MSPE" = ggplot(data=data, aes(x=Nodesize, y=Avg.MSPE,color=Method)) + geom_point() + geom_line(),
              "Avg.Wid"  = ggplot(data=data, aes(x=Nodesize, y=Avg.Wid, color=Method)) + geom_point() + geom_line(),
              "Avg.CR"   = ggplot(data=data, aes(x=Nodesize, y=Avg.CR,  color=Method)) + geom_point() + geom_line() 
              +geom_hline(yintercept=0.95) ))
}

### Function for plotting the stored data ###
Dplotter <- function(data){
  if(dim(data)[[2]]==2){
    return(list("X1/Y"=ggplot(data=data, aes(x=X1, y=Y)) + geom_point()
                )
           )}
  else                {
    return(list("X1/Y" =ggplot(data=data, aes(x=X1,  y=Y)) + geom_point(),
                "X2/Y" =ggplot(data=data, aes(x=X2,  y=Y)) + geom_point(),
                "X3/Y" =ggplot(data=data, aes(x=X3,  y=Y)) + geom_point(),
                "X4/Y" =ggplot(data=data, aes(x=X4,  y=Y)) + geom_point(),
                "X5/Y" =ggplot(data=data, aes(x=X5,  y=Y)) + geom_point(),
                
                "X6/Y" =ggplot(data=data, aes(x=X6,  y=Y)) + geom_point(),
                "X7/Y" =ggplot(data=data, aes(x=X7,  y=Y)) + geom_point(),
                "X8/Y" =ggplot(data=data, aes(x=X8,  y=Y)) + geom_point(),
                "X9/Y" =ggplot(data=data, aes(x=X9,  y=Y)) + geom_point(),
                "X10/Y"=ggplot(data=data, aes(x=X10, y=Y)) + geom_point()
                )
           )
  }
}

### plots  ###
plots <- list("eq1"=errors,"eq2"=errors,"eq3"=errors,"eq4"=errors,"eq5"=errors)
for (k in 1:length(eqAll)) {
  for(t in 1:length(errorAll)){
   plots[[k]][[t]] <- list( "Three Values"=plotter(data=dfs[[k]][[t]]), "data plot"=Dplotter(stored[[k]][[t]]))
  }
}


##### Compare Avg. Wid #####

### Function that returns ratios, their avg, min&max of the Avg.Width ### 
widAnalysis <- function(data){
  df <- data.frame(matrix(NA, nrow=dim(data)[[1]], ncol=dim(data)[[2]]))
  colnames(df) <-  c("methods", "nodesize","MSPE","Avg.Wid", "CR")
  df$methods  <- data$Method #methods
  df$nodesize <- data$Nodesize #nodesizes
  df[1:(dim(df)[[1]]/4), 3:4] <- data[1:(dim(df)[[1]]/4),2:3]
  for (k in 1:3) {
    df$MSPE   [((k*dim(data)[[1]]/4+1):((k+1)*dim(data)[[1]]/4))] <- data$Avg.MSPE[((k*dim(data)[[1]]/4+1)):((k+1)*dim(data)[[1]]/4)] / data$Avg.MSPE[1:(dim(data)[[1]]/4)]
    df$Avg.Wid[((k*dim(data)[[1]]/4+1):((k+1)*dim(data)[[1]]/4))] <- data$Avg.Wid [((k*dim(data)[[1]]/4+1)):((k+1)*dim(data)[[1]]/4)] / data$Avg.Wid[1:(dim(data)[[1]]/4)]
    #df$CR   [((k*dim(data)[[1]]/4+1):((k+1)*dim(data)[[1]]/4))] <- data$Avg.CR[((k*dim(data)[[1]]/4+1)):((k+1)*dim(data)[[1]]/4)] / data$Avg.CR[1:(dim(data)[[1]]/4)]
  }
  df$CR <- data$Avg.CR
  dfwanted <- df[order(df$nodesize),]
  
  return(dfwanted)
}

# Summarize: all the ratios are in this list, for each of the 30 combinations. Only remains to extract the rows with the optimal nodesize, once the data simulation finishes
width <- list("eq1"=lapply(X=dfs[[1]],widAnalysis),"eq2"=lapply(X=dfs[[2]],widAnalysis),"eq3"=lapply(X=dfs[[3]],widAnalysis),
              "eq4"=lapply(X=dfs[[4]],widAnalysis),"eq5"=lapply(X=dfs[[5]],widAnalysis)  )

# extract the ones with optimal nodesizes from each simulation #
table <- data.frame(matrix(NA, nrow=4*4, ncol=5))
colnames(table) <- colnames(width$Sim1)

table[1:4,  ] <- width$Sim1[width$Sim1$nodesize==ns1[(length(ns1)+1)/2],]
table[5:8,  ] <- width$Sim2[width$Sim2$nodesize==ns2[(length(ns2)+1)/2],]
table[9:12, ] <- width$Sim3[width$Sim3$nodesize==ns3[2],]
table[13:16,] <- width$Sim4[width$Sim4$nodesize==ns4[2],]

# ALSO APPLY MAKE THE PLOTS FOR THE 30 COMBINATIONS ONCE THE SIMULATION IS COMPLETE

