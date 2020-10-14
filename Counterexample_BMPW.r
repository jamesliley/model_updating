#################################################################
## Counterexample code to show better models can perform worse ##
## Sam Emerson, 6 Oct 2020                                     ##
#################################################################

f <- function(xs,xa,xl){
  return((1+exp(-xs-xa-0.01*xl))^(-1))
}

f.sample <-  function(xs,xa,xl){
  return(rbinom(1,1,f(xs,xa,xl)))
}

m <- function(p,N = 1000){
  S <- mvrnorm(N,mu = rep(0,3),Sigma = diag(3))
  absfun <- function(v){
    return(abs(f(v[1],v[2],v[3]) - p(v[1],v[2])))
  }
  V <- apply(S,1,absfun)
  return(mean(V))
}

training.expectation <- function(N=1000,n=100){
  D.initial <- data.frame(mvrnorm(n,mu = rep(0,3),Sigma = diag(3)))
  D.initial$y <- rep(0,n)
  for(i in 1:n){
    D.initial$y[i] <- f.sample(D.initial$X1[i],D.initial$X2[i],D.initial$X3[i])
  }
  p0.initial <- function(xs,xa){
    if(xs+xa>0){
      return(mean(D.initial[which(D.initial[,1] + D.initial[,2]>0),]$y))
    } else{
      return(mean(D.initial[which(D.initial[,1] + D.initial[,2]<=0),]$y))
    }
  }
  intervention <- function(xs,xa){
    return((xa+3)*(1-p0.initial(xs,xa)) + (xa-3)*p0.initial(xs,xa))
  }
  m.p0 <- rep(0,N)
  m.p1.NI <- rep(0,N)
  m.p1.I <- rep(0,N)
  for(j in 1:N){
    D <- data.frame(mvrnorm(n,mu = rep(0,3),Sigma = diag(3)))
    D$y <- rep(0,n)
    for(i in 1:n){
      D$y[i] <- f.sample(D$X1[i],D$X2[i],D$X3[i])
    }
    p0 <- function(xs,xa){
      if(xs+xa>0){
        return(mean(D[which(D[,1] + D[,2]>0),]$y))
      } else{
        return(mean(D[which(D[,1] + D[,2]<=0),]$y))
      }
    }
    m.p0[j] <- m(p0)
    beta <- as.vector(glm(y~X1+X2,family = binomial,data = D)$coefficients)
    p1.NI <- function(xs,xa){
      return((1+exp(-as.vector(c(1,xs,xa)%*%beta)))^(-1))
    }
    m.p1.NI[j] <- m(p1.NI)
    D$X4 <- rep(0,n)
    for(i in 1:n){
      D$X4[i] <- intervention(D$X1[i],D$X2[i])
    }
    D$y1 <- rep(0,n)
    for(i in 1:n){
      D$y1[i] <- f.sample(D$X1[i],D$X4[i],D$X3[i])
    }
    beta.1 <- as.vector(glm(y1~X1+X2,family = binomial,data = D)$coefficients)
    p1.I <- function(xs,xa){
      return((1+exp(-as.vector(c(1,xs,xa)%*%beta.1)))^(-1))
    }
    m.p1.I[j] <- m(p1.I)
  }
  return(list(D.initial,
              c(mean(D.initial[which(D.initial[,1] + D.initial[,2]>0),][,4]),
                mean(D.initial[which(D.initial[,1] + D.initial[,2]<=0),][,4])),
              c(mean(m.p0),mean(m.p1.NI),mean(m.p1.I))))
}