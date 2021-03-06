library(MASS)
library(glmnet)
library(ggplot2)
library(dplyr)

####### generate signal variables
signal_data <- function(numx, numxs, XY, XSY){
  if ( XY =="Low" && XSY == "Low"){
    beta <- matrix( c( rep(0.5, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
  } 
  if ( XY =="Low" && XSY == "Medium"){
    beta <- matrix( c( rep(0.5, numx), rep(1, numxs)), ncol = 1 ) #coefficient
  } 
  if ( XY =="Low" && XSY == "High"){
    beta <- matrix( c( rep(0.5, numx), rep(2, numxs)), ncol = 1 ) #coefficient
  } 
  if ( XY =="Medium" && XSY == "Low"){
    beta <- matrix( c( rep(1, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
  } 
  if ( XY =="Medium" && XSY == "Medium"){
    beta <- matrix( c( rep(1, numx), rep(1, numxs)), ncol = 1 ) #coefficient
  } 
  if ( XY =="Medium" && XSY == "High"){
    beta <- matrix( c( rep(1, numx), rep(2, numxs)), ncol = 1 ) #coefficient
  } 
  if ( XY =="High" && XSY == "Low"){
    beta <- matrix( c( rep(2, numx), rep(0.5, numxs)), ncol = 1 ) #coefficient
  } 
  if ( XY =="High" && XSY == "Medium"){
    beta <- matrix( c( rep(2, numx), rep(1, numxs)), ncol = 1 ) #coefficient
  } 
  if ( XY =="High" && XSY == "High"){
    beta <- matrix( c( rep(2, numx), rep(2, numxs)), ncol = 1 ) #coefficient
  } 
  mu <- rep(.1, (numx + numxs)) 
  sigma <- matrix(0.5, nrow =  (numx + numxs), 
                  ncol = (numx + numxs) )
  sigma <- sigma + diag(0.5, (numx + numxs)) 
  XX.star <- mvrnorm(500, mu, sigma)
  
  epsilon <- rnorm(500,0,.1) #generate noise
  
  Y <- XX.star %*% beta + epsilon #generate  return(list(X = X, Xstar = Xstar, Y = Y, XX.star = XX.star))
  df <- cbind(XX.star,Y)
  colnames(df) <- c(paste0("X", 1:numx),paste0("Xstar",1:numxs),"Y")
  return(df)
}

noise_data <- function(number){

  mu    <- rep(.1, number) 
  sigma <- matrix(0.2, ncol = number, nrow = number)
  sigma <- sigma + diag(0.8, number)
  
  noise <- mvrnorm(study*study_size, mu, sigma)
  
  return(data.frame(noise = noise))
}

multiple_sparse <- function(numx, numxs, numnoise, XY, XSY,
                     xxs.non.miss.num, noise.non.miss.num, sim.iter){
  MSE <- matrix(NA, nrow=sim.iter, ncol=5)
  colnames(MSE) <- c("Omit", "PL", "PP", "ML", "MP")
  num.impu.merge <- c()
  num.impu.pair  <- c()
  for(sim in 1:sim.iter){
    set.seed(sim)
    signal   <- signal_data(numx = numx, numxs = numxs, XY = XY, XSY = XSY)
    noise    <- noise_data(number = numnoise)
    data     <- as.matrix(cbind(noise, signal))
    
    index.full <- 1:500
    index1 <- sample(index.full, 100)
    index2 <- sample(index.full[-index1], 100)
    index3 <- sample(index.full[-c(index1, index2)], 100)
    index4 <- sample(index.full[-c(index1, index2, index3)], 100)
    index5 <- sample(index.full[-c(index1, index2, index3, index4)], 100)
    
    T1  <- as.matrix(data[index1,])
    T2  <- as.matrix(data[index2,])
    T3  <- as.matrix(data[index3,])
    T4  <- as.matrix(data[index4,])
    T5  <- as.matrix(data[index5,]) ### T5 as validation set
    
    T1.available <- c( (numnoise+1):(numnoise+2), 
                       sample( (numnoise+3):(numnoise+numx+numxs), xxs.non.miss.num), 
                       sample( c(1:numnoise), noise.non.miss.num) )
    T1.name      <- colnames(T1)[T1.available]
    T1.true      <- T1[, c(T1.available,ncol(T1))]
    colnames(T1.true)[ncol(T1.true)] <- "Y"
    
    T2.available <- c( (numnoise+1):(numnoise+2), 
                       sample( (numnoise+3):(numnoise+numx+numxs), xxs.non.miss.num), 
                       sample( c(1:numnoise), noise.non.miss.num) )
    T2.name      <- colnames(T2)[T2.available]
    T2.true      <- T2[, c(T2.available,ncol(T2))]
    colnames(T2.true)[ncol(T2.true)] <- "Y"
    
    T3.available <- c( (numnoise+1):(numnoise+2), 
                       sample( (numnoise+3):(numnoise+numx+numxs), xxs.non.miss.num), 
                       sample( c(1:numnoise), noise.non.miss.num) )
    T3.name      <- colnames(T3)[T3.available]
    T3.true      <- T3[, c(T3.available,ncol(T3))]
    colnames(T3.true)[ncol(T3.true)] <- "Y"
    
    T4.available <- c( (numnoise+1):(numnoise+2), 
                       sample( (numnoise+3):(numnoise+numx+numxs), xxs.non.miss.num), 
                       sample( c(1:numnoise), noise.non.miss.num) )
    T4.name      <- colnames(T4)[T4.available]
    T4.true      <- T4[, c(T4.available,ncol(T4))]
    colnames(T4.true)[ncol(T4.true)] <- "Y"
    
    name1.T <- T1.name
    name2.T <- T2.name
    name3.T <- T3.name
    name4.T <- T4.name
    
    Name    <- matrix( nrow=4, ncol  = length(name1.T),
                       c(name1.T, name2.T, name3.T, name4.T), byrow=TRUE)
    
    V   <- T5
    
    inter2 <- intersect(name1.T, name2.T)
    inter3 <- intersect(inter2, name3.T)
    inter4 <- intersect(inter3, name4.T)
    num.impu.merge <- append(num.impu.merge, length(inter4))
    
    inter4.y <- c(inter4, "Y")
    
    T.omit <- rbind(T1[,inter4.y], T2[,inter4.y],
                    T3[,inter4.y], T4[,inter4.y])
    #### Omitting
    X.train     <- as.matrix(T.omit[, inter4])
    Y.train     <- as.matrix(T.omit[, "Y"])
    cv_fit      <- cv.glmnet(X.train, Y.train, alpha = 1, 
                             lambda = 10^seq(3, -2, by = -.1), standardize = FALSE)
    opt_lambda  <- cv_fit$lambda.min
    fit.omit    <- glmnet(X.train, Y.train, alpha = 1, lambda = opt_lambda,standardize = FALSE)
    pred.omit   <- predict(fit.omit, as.matrix(V[,inter4]))
    MSE[sim,1]  <- mean((pred.omit -  V[,"Y"])^2)
    
    #### Imputation
    T1.impute   <- c()
    T2.impute   <- c()
    T3.impute   <- c()
    T4.impute   <- c()
    
    index.full  <- matrix(nrow = 4, ncol = 100,
                          c(index1, index2, index3, index4),
                          byrow = TRUE)
    name.full   <- Name
    combinations <- combn(4, 2)
    for(j in 1:ncol(combinations)){
      index.data1 <- index.full[combinations[,j][1],]
      index.data2 <- index.full[combinations[,j][2],]
      name1       <- name.full[combinations[,j][1],]
      name2       <- name.full[combinations[,j][2],]
      
      T.first     <- data[index.data1,][,name1]
      T.sec       <- data[index.data2,][,name2]
      names.inter <- intersect(name1, name2)
      num.impu.pair <- append(num.impu.pair, length(names.inter))
      
      name.first  <- name1[!(name1 %in% names.inter)]
      index       <- which(name1 %in% name.first)
      if(length(index) > 0){
        T.sec.impu  <- matrix(NA, ncol=length(index), nrow=100)
        colnames(T.sec.impu) <- name.first
        for(i in (1:length(index))){
          index.1         <- index[i]
          
          X.train         <- (T.first[, names.inter])
          Y.train         <- (T.first[, index.1])
          cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                       lambda = 10^seq(3, -2, by = -.1), 
                                       standardize = FALSE)
          opt_lambda      <- cv_fit$lambda.min
          fit.impu        <- glmnet(X.train, Y.train, alpha = 1, 
                                    lambda = opt_lambda,standardize = FALSE)
          T.sec.impu[,i]  <- predict(fit.impu, T.sec[, names.inter])
        }
        if(combinations[,j][2]==1){
          T1.impute   <- cbind(T1.impute, T.sec.impu)
        }
        if(combinations[,j][2]==2){
          T2.impute   <- cbind(T2.impute, T.sec.impu)
        }
        if(combinations[,j][2]==3){
          T3.impute   <- cbind(T3.impute, T.sec.impu)
        }
        if(combinations[,j][2]==4){
          T4.impute   <- cbind(T4.impute, T.sec.impu)
        }
      } 
      
      name.sec    <- name2[!(name2 %in% names.inter)]
      index       <- which(name2 %in% name.sec)
      if(length(index) > 0){
        T.first.impu      <- matrix(NA, ncol=length(index), nrow=100)
        colnames(T.first.impu) <- name.sec
        for(i in 1:length(index)){
          index.1         <- index[i]
          X.train         <- (T.sec[, names.inter])
          Y.train         <- (T.sec[, index.1])
          cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                       lambda = 10^seq(3, -2, by = -.1), 
                                       standardize = FALSE)
          opt_lambda      <- cv_fit$lambda.min
          fit.impu        <- glmnet(X.train, Y.train, alpha = 1, 
                                    lambda = opt_lambda,standardize = FALSE)
          T.first.impu[,i]  <- predict(fit.impu, T.first[, names.inter])
        }
        if(combinations[,j][1]==1){
          T1.impute   <- cbind(T1.impute, T.first.impu)
        }
        if(combinations[,j][1]==2){
          T2.impute   <- cbind(T2.impute, T.first.impu)
        }
        if(combinations[,j][1]==3){
          T3.impute   <- cbind(T3.impute, T.first.impu)
        }
        if(combinations[,j][1]==4){
          T4.impute   <- cbind(T4.impute, T.first.impu)
        }
      }
    }
    
    name.star  <- paste0("Xstar",1:numxs)
    name.x     <- paste0("X", 1:numx)
    name.noise <- paste0("noise.", 1:numnoise)
    
    name1.impute <- colnames(T1.impute)
    name2.impute <- colnames(T2.impute)
    name3.impute <- colnames(T3.impute)
    name4.impute <- colnames(T4.impute)
    
    T1.impute2 <- c()
    T2.impute2 <- c()
    T3.impute2 <- c()
    T4.impute2 <- c()
    
    for(name in name.star){
      if(length(which(name1.impute %in% name))>1){
        temp       <- matrix(rowMeans(T1.impute[,which(name1.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2, temp)
      } else if(length(which(name1.impute %in% name))==1){
        temp       <- matrix(T1.impute[,which(name1.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2,temp)
      }
      
      if(length(which(name2.impute %in% name))>1){
        temp       <- matrix(rowMeans(T2.impute[,which(name2.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2, temp)
      } else if(length(which(name2.impute %in% name))==1){
        temp       <- matrix(T2.impute[,which(name2.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2,temp)
      }
      
      if(length(which(name3.impute %in% name))>1){
        temp       <- matrix(rowMeans(T3.impute[,which(name3.impute %in% name)]),ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2, temp)
      } else if(length(which(name3.impute %in% name))==1){
        temp       <- matrix(T3.impute[,which(name3.impute %in% name)],ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2,temp)
      }
      
      if(length(which(name4.impute %in% name))>1){
        temp       <- matrix(rowMeans(T4.impute[,which(name4.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2, temp)
      } else if(length(which(name4.impute %in% name))==1){
        temp       <- matrix(T4.impute[,which(name4.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2,temp)
      }
      
    }
    
    for(name in name.x){
      if(length(which(name1.impute %in% name))>1){
        temp       <- matrix(rowMeans(T1.impute[,which(name1.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2, temp)
      } else if(length(which(name1.impute %in% name))==1){
        temp       <- matrix(T1.impute[,which(name1.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2,temp)
      }
      
      if(length(which(name2.impute %in% name))>1){
        temp       <- matrix(rowMeans(T2.impute[,which(name2.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2, temp)
      } else if(length(which(name2.impute %in% name))==1){
        temp       <- matrix(T2.impute[,which(name2.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2,temp)
      }
      
      if(length(which(name3.impute %in% name))>1){
        temp       <- matrix(rowMeans(T3.impute[,which(name3.impute %in% name)]),ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2, temp)
      } else if(length(which(name3.impute %in% name))==1){
        temp       <- matrix(T3.impute[,which(name3.impute %in% name)],ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2,temp)
      }
      
      if(length(which(name4.impute %in% name))>1){
        temp       <- matrix(rowMeans(T4.impute[,which(name4.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2, temp)
      } else if(length(which(name4.impute %in% name))==1){
        temp       <- matrix(T4.impute[,which(name4.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2,temp)
      }
      
    }
    
    for(name in name.noise){
      if(length(which(name1.impute %in% name))>1){
        temp       <- matrix(rowMeans(T1.impute[,which(name1.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2, temp)
      } else if(length(which(name1.impute %in% name))==1){
        temp       <- matrix(T1.impute[,which(name1.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2,temp)
      }
      
      if(length(which(name2.impute %in% name))>1){
        temp       <- matrix(rowMeans(T2.impute[,which(name2.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2, temp)
      } else if(length(which(name2.impute %in% name))==1){
        temp       <- matrix(T2.impute[,which(name2.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2,temp)
      }
      
      if(length(which(name3.impute %in% name))>1){
        temp       <- matrix(rowMeans(T3.impute[,which(name3.impute %in% name)]),ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2, temp)
      } else if(length(which(name3.impute %in% name))==1){
        temp       <- matrix(T3.impute[,which(name3.impute %in% name)],ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2,temp)
      }
      
      if(length(which(name4.impute %in% name))>1){
        temp       <- matrix(rowMeans(T4.impute[,which(name4.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2, temp)
      } else if(length(which(name4.impute %in% name))==1){
        temp       <- matrix(T4.impute[,which(name4.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2,temp)
      }
      
    }
    
    
    T1.true2 <- cbind(T1.true, T1.impute2)
    T2.true2 <- cbind(T2.true, T2.impute2)
    T3.true2 <- cbind(T3.true, T3.impute2)
    T4.true2 <- cbind(T4.true, T4.impute2)
    
    T1.true3 <- T1.true2[,sort(colnames(T1.true2))]
    T2.true3 <- T2.true2[,sort(colnames(T2.true2))]
    T3.true3 <- T3.true2[,sort(colnames(T3.true2))]
    T4.true3 <- T4.true2[,sort(colnames(T4.true2))]
    
    
    #### Imputation
    df.imput <- rbind(T1.true3, T2.true3, T3.true3, T4.true3)
    variables <- colnames(df.imput)[-length(colnames(df.imput))]
    
    X.train     <- as.matrix(df.imput[,variables])
    Y.train     <- as.matrix(df.imput[,length(colnames(df.imput))])
    cv_fit      <- cv.glmnet(X.train, Y.train, alpha = 1, 
                             lambda = 10^seq(3, -2, by = -.1), standardize = FALSE)
    opt_lambda  <- cv_fit$lambda.min
    fit.impu    <- glmnet(X.train, Y.train, alpha = 1, lambda = opt_lambda,standardize = FALSE)
    pred.impu   <- predict(fit.impu, as.matrix(V[,variables]))
    MSE[sim,2]  <- mean((pred.impu -  V[,"Y"])^2)
    
    #### Polynomial KT
    T1.impute   <- c()
    T2.impute   <- c()
    T3.impute   <- c()
    T4.impute   <- c()
    
    combinations <- combn(4, 2)
    for(j in 1:ncol(combinations)){
      index.data1 <- index.full[combinations[,j][1],]
      index.data2 <- index.full[combinations[,j][2],]
      name1       <- name.full[combinations[,j][1],]
      name2       <- name.full[combinations[,j][2],]
      
      T.first     <- data[index.data1,][,name1]
      T.sec       <- data[index.data2,][,name2]
      names.inter <- intersect(name1, name2)
      
      name.first  <- name1[!(name1 %in% names.inter)]
      index       <- which(name1 %in% name.first)
      if(length(index) > 0){
        T.sec.impu  <- matrix(NA, ncol=length(index), nrow=100)
        colnames(T.sec.impu) <- name.first
        for(i in (1:length(index))){
          index.1         <- index[i]
          X.train         <- (cbind(T.first[, names.inter], T.first[, names.inter]^2))
          colnames(X.train) <- c(names.inter,paste0( names.inter,"sq"))
          Y.train         <- (T.first[, index.1])
          cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                       lambda = 10^seq(3, -2, by = -.1), 
                                       standardize = FALSE)
          opt_lambda      <- cv_fit$lambda.min
          fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                    lambda = opt_lambda,standardize = FALSE)
          X.for.pred      <- (cbind(T.sec[, names.inter], T.sec[, names.inter]^2))
          colnames(X.for.pred) <- c(names.inter,paste0( names.inter,"sq"))
          T.sec.impu[,i]  <- predict(fit.poly, X.for.pred)
        }
        if(combinations[,j][2]==1){
          T1.impute   <- cbind(T1.impute, T.sec.impu)
        }
        if(combinations[,j][2]==2){
          T2.impute   <- cbind(T2.impute, T.sec.impu)
        }
        if(combinations[,j][2]==3){
          T3.impute   <- cbind(T3.impute, T.sec.impu)
        }
        if(combinations[,j][2]==4){
          T4.impute   <- cbind(T4.impute, T.sec.impu)
        }
      } 
      
      name.sec    <- name2[!(name2 %in% names.inter)]
      index       <- which(name2 %in% name.sec)
      if(length(index) > 0){
        T.first.impu      <- matrix(NA, ncol=length(index), nrow=100)
        colnames(T.first.impu) <- name.sec
        for(i in 1:length(index)){
          index.1         <- index[i]
          X.train         <- (cbind(T.sec[, names.inter], T.sec[, names.inter]^2))
          colnames(X.train) <- c(names.inter,paste0( names.inter,"sq"))
          Y.train         <- (T.sec[, index.1])
          cv_fit          <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                       lambda = 10^seq(3, -2, by = -.1), 
                                       standardize = FALSE)
          opt_lambda      <- cv_fit$lambda.min
          fit.poly        <- glmnet(X.train, Y.train, alpha = 1, 
                                    lambda = opt_lambda,standardize = FALSE)
          X.for.pred      <- (cbind(T.first[, names.inter], T.first[, names.inter]^2))
          colnames(X.for.pred) <- c(names.inter,paste0( names.inter,"sq"))
          T.first.impu[,i]<- predict(fit.poly, X.for.pred)
        }
        if(combinations[,j][1]==1){
          T1.impute   <- cbind(T1.impute, T.first.impu)
        }
        if(combinations[,j][1]==2){
          T2.impute   <- cbind(T2.impute, T.first.impu)
        }
        if(combinations[,j][1]==3){
          T3.impute   <- cbind(T3.impute, T.first.impu)
        }
        if(combinations[,j][1]==4){
          T4.impute   <- cbind(T4.impute, T.first.impu)
        }
      }
    }
    
    name.star  <- paste0("Xstar",1:numxs)
    name.x     <- paste0("X", 1:numx)
    name.noise <- paste0("noise.", 1:numnoise)
    
    name1.impute <- colnames(T1.impute)
    name2.impute <- colnames(T2.impute)
    name3.impute <- colnames(T3.impute)
    name4.impute <- colnames(T4.impute)
    
    T1.impute2 <- c()
    T2.impute2 <- c()
    T3.impute2 <- c()
    T4.impute2 <- c()
    
    for(name in name.star){
      if(length(which(name1.impute %in% name))>1){
        temp       <- matrix(rowMeans(T1.impute[,which(name1.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2, temp)
      } else if(length(which(name1.impute %in% name))==1){
        temp       <- matrix(T1.impute[,which(name1.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2,temp)
      }
      
      if(length(which(name2.impute %in% name))>1){
        temp       <- matrix(rowMeans(T2.impute[,which(name2.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2, temp)
      } else if(length(which(name2.impute %in% name))==1){
        temp       <- matrix(T2.impute[,which(name2.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2,temp)
      }
      
      if(length(which(name3.impute %in% name))>1){
        temp       <- matrix(rowMeans(T3.impute[,which(name3.impute %in% name)]),ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2, temp)
      } else if(length(which(name3.impute %in% name))==1){
        temp       <- matrix(T3.impute[,which(name3.impute %in% name)],ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2,temp)
      }
      
      if(length(which(name4.impute %in% name))>1){
        temp       <- matrix(rowMeans(T4.impute[,which(name4.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2, temp)
      } else if(length(which(name4.impute %in% name))==1){
        temp       <- matrix(T4.impute[,which(name4.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2,temp)
      }
      
    }
    
    for(name in name.x){
      if(length(which(name1.impute %in% name))>1){
        temp       <- matrix(rowMeans(T1.impute[,which(name1.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2, temp)
      } else if(length(which(name1.impute %in% name))==1){
        temp       <- matrix(T1.impute[,which(name1.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2,temp)
      }
      
      if(length(which(name2.impute %in% name))>1){
        temp       <- matrix(rowMeans(T2.impute[,which(name2.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2, temp)
      } else if(length(which(name2.impute %in% name))==1){
        temp       <- matrix(T2.impute[,which(name2.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2,temp)
      }
      
      if(length(which(name3.impute %in% name))>1){
        temp       <- matrix(rowMeans(T3.impute[,which(name3.impute %in% name)]),ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2, temp)
      } else if(length(which(name3.impute %in% name))==1){
        temp       <- matrix(T3.impute[,which(name3.impute %in% name)],ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2,temp)
      }
      
      if(length(which(name4.impute %in% name))>1){
        temp       <- matrix(rowMeans(T4.impute[,which(name4.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2, temp)
      } else if(length(which(name4.impute %in% name))==1){
        temp       <- matrix(T4.impute[,which(name4.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2,temp)
      }
      
    }
    
    for(name in name.noise){
      if(length(which(name1.impute %in% name))>1){
        temp       <- matrix(rowMeans(T1.impute[,which(name1.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2, temp)
      } else if(length(which(name1.impute %in% name))==1){
        temp       <- matrix(T1.impute[,which(name1.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T1.impute2 <- cbind(T1.impute2,temp)
      }
      
      if(length(which(name2.impute %in% name))>1){
        temp       <- matrix(rowMeans(T2.impute[,which(name2.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2, temp)
      } else if(length(which(name2.impute %in% name))==1){
        temp       <- matrix(T2.impute[,which(name2.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T2.impute2 <- cbind(T2.impute2,temp)
      }
      
      if(length(which(name3.impute %in% name))>1){
        temp       <- matrix(rowMeans(T3.impute[,which(name3.impute %in% name)]),ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2, temp)
      } else if(length(which(name3.impute %in% name))==1){
        temp       <- matrix(T3.impute[,which(name3.impute %in% name)],ncol=1)
        colnames(temp) <- name
        T3.impute2 <- cbind(T3.impute2,temp)
      }
      
      if(length(which(name4.impute %in% name))>1){
        temp       <- matrix(rowMeans(T4.impute[,which(name4.impute %in% name)]), ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2, temp)
      } else if(length(which(name4.impute %in% name))==1){
        temp       <- matrix(T4.impute[,which(name4.impute %in% name)], ncol=1)
        colnames(temp) <- name
        T4.impute2 <- cbind(T4.impute2,temp)
      }
      
    }
    
    
    T1.true2 <- cbind(T1.true, T1.impute2)
    T2.true2 <- cbind(T2.true, T2.impute2)
    T3.true2 <- cbind(T3.true, T3.impute2)
    T4.true2 <- cbind(T4.true, T4.impute2)
    
    T1.true3 <- T1.true2[,sort(colnames(T1.true2))]
    T2.true3 <- T2.true2[,sort(colnames(T2.true2))]
    T3.true3 <- T3.true2[,sort(colnames(T3.true2))]
    T4.true3 <- T4.true2[,sort(colnames(T4.true2))]
    
    #### Polynomial KT
    df.poly <- rbind(T1.true3, T2.true3, T3.true3, T4.true3)
    variables <- colnames(df.poly)[-length(colnames(df.poly))]
    
    X.train     <- as.matrix(df.poly[,variables])
    Y.train     <- as.matrix(df.poly[,length(colnames(df.poly))])
    cv_fit      <- cv.glmnet(X.train, Y.train, alpha = 1, 
                             lambda = 10^seq(3, -2, by = -.1), standardize = FALSE)
    opt_lambda  <- cv_fit$lambda.min
    fit.poly    <- glmnet(X.train, Y.train, alpha = 1, lambda = opt_lambda,standardize = FALSE)
    pred.poly   <- predict(fit.poly, V[,variables])
    MSE[sim,3]  <- mean((pred.poly -  V[,"Y"])^2)
    
    #### merge
    T1.true.df <- data.frame(T1.true)
    T2.true.df <- data.frame(T2.true)
    T3.true.df <- data.frame(T3.true)
    T4.true.df <- data.frame(T4.true)
    
    T1234.merge <- bind_rows(T1.true.df, T2.true.df, T3.true.df, T4.true.df)
    index.Y     <- which(colnames(T1234.merge)=="Y")
    T1234.merge <- cbind(T1234.merge[,-index.Y], T1234.merge[,index.Y])
    colnames(T1234.merge)[ncol(T1234.merge)] <- "Y"
    index.imput <-  1: (ncol(T1234.merge)-1)
    T1234.merge.linear <- T1234.merge
    T1234.merge.poly   <- T1234.merge
    
    for(kk in index.imput){
      row.index.miss <- which(is.na(T1234.merge[,kk])==TRUE)
      if(length(row.index.miss)!=0){
        #### Linear imputation
        X.train      <- as.matrix(T1234.merge[-row.index.miss, inter4 ])
        Y.train      <- T1234.merge[-row.index.miss, kk ]
        cv_fit       <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                  lambda = 10^seq(3, -2, by = -.1), 
                                  standardize = FALSE)
        opt_lambda   <- cv_fit$lambda.min
        fit.impu     <- glmnet(X.train, Y.train, alpha = 1, 
                               lambda = opt_lambda,standardize = FALSE)
        T1234.merge.linear[row.index.miss, kk]  <- predict(fit.impu, 
                                                           as.matrix(T1234.merge[row.index.miss, inter4]))
        
        
        #### polynomial imputation
        X.train      <- cbind(as.matrix(T1234.merge[-row.index.miss, inter4 ]), 
                              as.matrix(T1234.merge[-row.index.miss, inter4 ]^2))
        colnames(X.train) <- c(inter4,paste0( inter4,"sq"))
        Y.train      <- T1234.merge[-row.index.miss, kk ]
        cv_fit       <- cv.glmnet(X.train, Y.train, alpha = 1, 
                                  lambda = 10^seq(3, -2, by = -.1), 
                                  standardize = FALSE)
        opt_lambda   <- cv_fit$lambda.min
        fit.impu     <- glmnet(X.train, Y.train, alpha = 1, 
                               lambda = opt_lambda,standardize = FALSE)
        X.for.pred   <- cbind(as.matrix(T1234.merge[row.index.miss, inter4]), 
                              as.matrix(T1234.merge[row.index.miss, inter4]^2))
        colnames(X.for.pred) <- c(inter4,paste0( inter4,"sq"))
        T1234.merge.poly[row.index.miss, kk]  <- predict(fit.impu, 
                                                         X.for.pred)
        
      }
    }
    #### predicting on merged with linear imputation
    X.train     <- as.matrix(T1234.merge.linear[,-ncol(T1234.merge.linear)])
    Y.train     <- as.matrix(T1234.merge.linear[,ncol(T1234.merge.linear)])
    cv_fit      <- cv.glmnet(X.train, Y.train, alpha = 1, 
                             lambda = 10^seq(3, -2, by = -.1), standardize = FALSE)
    opt_lambda  <- cv_fit$lambda.min
    fit.merged.linear    <- glmnet(X.train, Y.train, alpha = 1, lambda = opt_lambda,standardize = FALSE)
    pred.merged.linear   <- predict(fit.merged.linear, V[,colnames(X.train)])
    MSE[sim,4]  <- mean((pred.merged.linear -  V[,"Y"])^2)
    #### predicting on merged with polynomial imputation
    X.train     <- as.matrix(T1234.merge.poly[,-ncol(T1234.merge.linear)])
    Y.train     <- as.matrix(T1234.merge.poly[,ncol(T1234.merge.linear)])
    cv_fit      <- cv.glmnet(X.train, Y.train, alpha = 1, 
                             lambda = 10^seq(3, -2, by = -.1), standardize = FALSE)
    opt_lambda  <- cv_fit$lambda.min
    fit.merged.poly    <- glmnet(X.train, Y.train, alpha = 1, lambda = opt_lambda,standardize = FALSE)
    pred.merged.poly  <- predict(fit.merged.poly, V[,colnames(X.train)])
    MSE[sim,5]  <- mean((pred.merged.poly -  V[,"Y"])^2)
    
  }
  return(list(MSE=MSE,
              num.impu.pair = num.impu.pair,
              num.impu.merge = num.impu.merge))
}
