#################################################################
## Counterexample code to show better models can perform worse ##
## Samuel Emerson, 6 Oct 2020                                  ##
#################################################################

#Note the package 'MASS' must be installed.

# Fuction to generate the true risk score f:
f <- function(xs,xa){
  return((1+exp(-xs-xa))^(-1))
}

# Function which samples the response given f:
f.sample <-  function(xs,xa){
  return(rbinom(1,1,f(xs,xa)))
}

#Function which transforms the actionable covariate xa based on a risk score p:
g.a <- function(p,xa){
  return((xa+3)*(1-p) + (xa-3)*p)
}

# Function which approximates the function m_{tilde{f_0}} given an estimate risk score function at epoch 0 (p.0):
m.tildef.0 <- function(p.0,N = 1000){
  S <- mvrnorm(N,mu = rep(0,2),Sigma = diag(2)) #We first sample N values from X_{0}(0)
  absfun <- function(v){ #This function simply gives the distance between tilde{f_0} and rho at a particular value of X_{0}(0)
    return(abs(f(v[1],v[2]) - p.0(v[1],v[2])))
  }
  V <- apply(S,1,absfun)
  return(mean(V)) #By returning the mean we give an Monte Carlo estimate of m_{tilde{f_e}}
}

# Function which approximates the function m_{tilde{f_e}} (e>0) given an estimate risk score function at epoch e (p.current) and epoch e-1 (p.previous):
m.tildef.e <- function(p.current,p.previous,N = 1000){
  S <- mvrnorm(N,mu = rep(0,2),Sigma = diag(2)) #We first sample N values from X_{e}(0)
  absfun <- function(v){ #This function simply gives the distance between tilde{f_e} and rho at a particular value of X_{e}(0)
    return(abs(f(v[1],g.a(p.previous(v[1],v[2]),v[2])) - p.current(v[1],v[2])))
  }
  V <- apply(S,1,absfun)
  return(mean(V)) #By returning the mean we give an Monte Carlo estimate of m_{tilde{f_e}}
}


# Function for generation of training data and the subsequent risk score at epoch 0 (i.e. p0) given this training data:
initial.set.up <- function(n=100){
  D.initial <- data.frame(mvrnorm(n,mu = rep(0,2),Sigma = diag(2)))
  D.initial$y <- rep(0,n)
  for(i in 1:n){
    D.initial$y[i] <- f.sample(D.initial$X1[i],D.initial$X2[i])
  } # This forms the initial training data sample at epoch 0
  p0.initial <- function(xs,xa){
    if(xs+xa>0){
      return(mean(D.initial[which(D.initial[,1] + D.initial[,2]>0),]$y))
    } else{
      return(mean(D.initial[which(D.initial[,1] + D.initial[,2]<=0),]$y))
    }
  } #From this training data p0.initial is the resulting risk score function
  return(list(D.initial,p0.initial,c(mean(D.initial[which(D.initial[,1] + D.initial[,2]>0),][,3]),
                                     mean(D.initial[which(D.initial[,1] + D.initial[,2]<=0),][,3])))) #We return; the simulated training data set at epoch 0, the risk function at epoch 0, and a vector of the two risk scores (one for if xs+xa>0 and one for if xs+xa<=0) given in the risk function at epoch 0
}

# Function for approximation of the expectation, with respect to the training data, of m_{tilde{f_0}}(p0):
training.expectation.p0 <- function(N=1000,n=100){
  m.p0 <- rep(0,N)
  # In order to gain an approximation we provide N training data samples and average over the resulting N m_{tilde{f_0}}(p0) values
  for(j in 1:N){
    D <- data.frame(mvrnorm(n,mu = rep(0,2),Sigma = diag(2)))
    D$y <- rep(0,n) #To obtain one of the samples for the approximation, we first sample a training data set X^{*}_{0}
    for(i in 1:n){
      D$y[i] <- f.sample(D$X1[i],D$X2[i])
    }
    p0 <- function(xs,xa){ #We then obtain a risk score function at epoch 0 from this data set
      if(xs+xa>0){
        return(mean(D[which(D[,1] + D[,2]>0),]$y))
      } else{
        return(mean(D[which(D[,1] + D[,2]<=0),]$y))
      }
    }
    m.p0[j] <- m.tildef.0(p0) #Then we calculate m_{tilde{f_0}} given this function
  }
  return(mean(m.p0)) #We return the a Monte Carlo estimate for m_{tilde{f_0}} of the risk function at epoch 0
}

#Function for approximation of the expectation, with respect to the training data, of m_{tilde{f_0}}(p1) when no interventions are made:
training.expectation.p1.NI <- function(N=1000,n=100){
  m.p1.NI <- rep(0,N)
  # In order to gain an approximation we provide N training data samples and average over the resulting N m_{tilde{f_0}}(p1) values
  for(j in 1:N){
    D <- data.frame(mvrnorm(n,mu = rep(0,2),Sigma = diag(2)))
    D$y <- rep(0,n) #To obtain one of the samples for the approximation, we first sample a training data set X^{*}_{1}
    for(i in 1:n){
      D$y[i] <- f.sample(D$X1[i],D$X2[i])
    }
    beta <- as.vector(glm(y~X1+X2,family = binomial,data = D)$coefficients) #The logistic regression coefficients are then obtained for the model at epoch 1
    p1.NI <- function(xs,xa){
      return((1+exp(-as.vector(c(1,xs,xa)%*%beta)))^(-1))
    }
    m.p1.NI[j] <- m.tildef.0(p1.NI) #By creating a risk score function at epoch 1, we can now calculate m_{tilde{f_0}} given this function
  }
  return(mean(m.p1.NI)) #We return the a Monte Carlo estimate for m_{tilde{f_0}} of the risk function at epoch 1
}

#Function for approximation of the expectation, with respect to the training data, of m_{tilde{f_1}}(p1) when interventions are made:
training.expectation.p1.I <- function(N=1000,n=100,p0.initial){ #Here p0.initial would be the second list element in the output of initial.set.up(n)
  m.p1.I <- rep(0,N)
  # In order to gain an approximation we provide N training data samples, intervene and then average over the resulting N m_{tilde{f_1}}(p1) values
  for(j in 1:N){
    D <- data.frame(mvrnorm(n,mu = rep(0,2),Sigma = diag(2))) #To obtain one of the samples for the approximation, we first sample a training data set X^{*}_{1}
    D$X3 <- rep(0,n)
    for(i in 1:n){
      D$X3[i] <- g.a(p0.initial(D$X1[i],D$X2[i]),D$X2[i])
    } #Here we intervene of the actionable covariates of the training data using the initial risk score p0.initial to form the temporary X^{*}_{1}(1)
    D$y <- rep(0,n)
    for(i in 1:n){
      D$y[i] <- f.sample(D$X1[i],D$X3[i])
    } #Our response training data at epoch 1, Y_{1}, is then obtained through sampling using X^{*}_{1}(1) as an input
    beta <- as.vector(glm(y~X1+X2,family = binomial,data = D)$coefficients) #The logistic regression coefficients are then obtained for the model at epoch 1
    p1.I <- function(xs,xa){
      return((1+exp(-as.vector(c(1,xs,xa)%*%beta)))^(-1))
    }
    m.p1.I[j] <- m.tildef.e(p1.I,p0.initial) #By creating a risk score function at epoch 1, we can now calculate m_{tilde{f_1}} given this function and the previous risk function
  }
  return(mean(m.p1.I)) #We return the a Monte Carlo estimate for m_{tilde{f_1}} of the risk function at epoch 1
}

#Function for approximation of the expectation, with respect to the training data, of m_{tilde{f_0}}(p1) when interventions are made:
training.expectation.p1.I.mftilde0 <- function(N=1000,n=100,p0.initial){ #Here p0.initial would be the second list element in the output of initial.set.up(n)
  m.p1.I <- rep(0,N)
  # In order to gain an approximation we provide N training data samples, intervene and then average over the resulting N m_{tilde{f_0}}(p1) values
  for(j in 1:N){
    D <- data.frame(mvrnorm(n,mu = rep(0,2),Sigma = diag(2))) #To obtain one of the samples for the approximation, we first sample a training data set X^{*}_{1}
    D$X3 <- rep(0,n)
    for(i in 1:n){
      D$X3[i] <- g.a(p0.initial(D$X1[i],D$X2[i]),D$X2[i])
    } #Here we intervene of the actionable covariates of the training data using the initial risk score p0.initial to form the temporary X^{*}_{1}(1)
    D$y <- rep(0,n)
    for(i in 1:n){
      D$y[i] <- f.sample(D$X1[i],D$X3[i])
    } #Our response training data at epoch 1, Y_{1}, is then obtained through sampling using X^{*}_{1}(1) as an input
    beta <- as.vector(glm(y~X1+X2,family = binomial,data = D)$coefficients) #The logistic regression coefficients are then obtained for the model at epoch 1
    p1.I <- function(xs,xa){
      return((1+exp(-as.vector(c(1,xs,xa)%*%beta)))^(-1))
    }
    m.p1.I[j] <- m.tildef.0(p1.I) #By creating a risk score function at epoch 1, we can now calculate m_{tilde{f_0}} given this function and the previous risk function
  }
  return(mean(m.p1.I)) #We return the a Monte Carlo estimate for m_{tilde{f_0}} of the risk function at epoch 1
}
