library("SQUAREM")
library("Rcpp")
## C++ implementation of SFS function
sourceCpp("code/sampleSFS.cpp")
#################################################################
## Finding SVD solutions alpha from samples of matrices P and E
#################################################################
samplesToSVD = function(Presults, Eresults, N){
  D = t(Presults[1:N,])%*%Eresults[1:N,]
  svd.D = svd(t(D), nu = N, nv = N)
  svdV = svd.D$v
  svdU = svd.D$u
  PT.points = matrix(0, nrow = nrow(Presults), ncol = N-1)
  TE.points = matrix(0, nrow = nrow(Eresults), ncol = N-1)
  for(i in 1:(nrow(Presults)/N)){
    p = Presults[(i*N-(N-1)):(i*N),]
    Tmat = p%*%svdV  # find Tmat
    Tmat = Tmat/Tmat[,1]
    PT.points[(i*3-2):(i*3),] = Tmat[,c(2,3)]
    e = Eresults[(i*N-(N-1)):(i*N),]
    Tmat = e%*%svdU  # find Tmat
    Tmat = Tmat/Tmat[,1]
    TE.points[(i*3-2):(i*3),] = Tmat[,c(2,3)]
  }
  Output = list()
  Output$PT.points = PT.points
  Output$TE.points = TE.points
  return(Output)
}

#######################################################################
## Function for calculating the generalised Kullback Leibler divergence
#######################################################################
gkl.dev <- function(y, mu){
  r <- mu
  p <- which(y > 0)
  r[p] <- (y * (log(y)- log(mu)) - y + mu)[p]
  return(sum(r))
}

###############################################################################
## Function factorizing M into two matrices P and E of
## dimension ncol(M) x N and N x nrow(M). 
## The objective function is the generalized Kullbach-Leibler divergence(GKLD).
##
## Input: M     - non-negative data matrix of size
##        N     - Small dimension of the two new matrices
##        tol   - Change of GKLD when algorithm is stopped 
##        seed  - Vector of random seeds to initialize the matrices
##
## Output: P   - Non-negative matrix of dimension ncol(V) x K, with columns summing to one
##         E   - Non-negative matrix of dimension K x nrow(V), where rows sum to one  
##         gkl - Smallest Value of the generalized kullbach leibler
###############################################################################
NMFPoisEM <- function(M,N,tol,seed, arrange = TRUE){
  K <- dim(M)[1]  # patients
  G <- dim(M)[2]  # mutations
  
  div <- rep(0,length(seed)) # vector of different GKLD values
  Plist <- list()            # list of P matrices
  Elist <- list()            # list of E matrices
  
  
  for(i in 1:length(seed)){ 
    set.seed(seed[i])
    
    P <- matrix(runif(K*N), nrow = K, ncol = N)  # Initialize P
    E <- matrix(runif(N*G), nrow = N, ncol = G)  # Initialize E
    
    GKLold <- 0
    
    repeat{
      PE <- P%*%E
      P <- P * ((M/PE) %*% t(E))      # update of signatures
      P <- P %*% diag(1/colSums(P))   # make sure the columns sum to one
      
      PE <- P%*%E
      E <- E * (t(P) %*% (M/PE))      # update of exposures
      
      
      GKL <- gkl.dev(as.vector(M),as.vector(P%*%E)) # GKLD value
      
      print(GKL)                                    # print GKLD value
      
      if(abs(GKLold - GKL) < tol){break}            # stop iterating if GKLD change less than tol
      GKLold <- GKL
    }
    
    Plist[[i]] <- P # signatures
    Elist[[i]] <- E # exposures
    div[i] <- GKL   # final generalized Kullback-Leibler divergence
  }
  
  best <- which.min(div) # Smallest GKLD value
  P = Plist[[best]]
  E = Elist[[best]]
  
  if(arrange == TRUE){
    idx = order(rowSums(E),decreasing = TRUE)
    P = P[,idx]
    E = E[idx,]
  }
  
  Output <- list()
  Output$P <-  P
  Output$E <-  E
  Output$gkl <- div[best]
  
  return(Output)
}

#################################################################################
## Function factorizing M into two matrices P and E of
## dimension ncol(M) x N and N x nrow(M) with the acceleration of SQUAREM.
## The objective function is the generalized Kullbach-Leibler divergence(GKLD).
## 
##
## Input: M     - non-negative data matrix of size
##        N     - Small dimension of the two new matrices
##        tol   - Maximum change of P and E when stopping 
##        seed  - Vector of random seeds to initialize the matrices
##
## Output: P   - Non-negative matrix of dimension ncol(V) x K, with columns summing to one
##         E   - Non-negative matrix of dimension K x nrow(V), where rows sum to one  
##         gkl - Smallest Value of the generalized kullbach leibler
####################################################################################
NMFPoisEMsquarem = function(M,N,seed, arrange = TRUE, tol = 1e-5){
  K <- dim(M)[1]  # mutations
  G <- dim(M)[2]  # patients
  
  div <- rep(0,length(seed)) # vector of different GKLD values
  Plist <- list()            # list of P matrices
  Elist <- list()            # list of E matrices
  reslist <- list()
  
  poisson.em = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)
    
    PE <- P%*%E
    P <- P * ((M/PE) %*% t(E))      # update of signatures
    P <- P %*% diag(1/colSums(P))   # make sure the columns sum to one
    
    PE <- P%*%E
    E <- E * (t(P) %*% (M/PE))      # update of exposures
    
    par = c(as.vector(P),as.vector(E))
    par[par <= 0] = 1e-10
    return(log(par))
  }
  
  gklobj = function(x){
    x = exp(x)
    P = matrix(x[1:(K*N)], nrow = K, ncol = N)
    E = matrix(x[-c(1:(K*N))], nrow = N, ncol = G)
    
    GKL <- gkl.dev(as.vector(M),as.vector(P%*%E)) # GKLD value
    
    return(GKL)
  }
  
  for(i in 1:length(seed)){ 
    set.seed(seed[i])
    
    P <- matrix(runif(K*N), nrow = K, ncol = N)  # Initialize P
    E <- matrix(runif(N*G), nrow = N, ncol = G)  # Initialize E
    
    init = log(c(as.vector(P),as.vector(E)))
    sres = squarem(init, fixptfn = poisson.em, objfn = gklobj, control = list(tol = tol))
    
    P = matrix(exp(sres$par[1:(K*N)]), nrow = K, ncol = N)
    E = matrix(exp(sres$par[-c(1:(K*N))]), nrow = N, ncol = G)
    
    Plist[[i]] <- P # signatures
    Elist[[i]] <- E # exposures
    div[i] <- gklobj(sres$par)   # final generalized Kullback-Leibler divergence
    reslist[[i]] = sres
  }
  
  best <- which.min(div) # Smallest GKLD value
  P = Plist[[best]]
  E = Elist[[best]]
  
  if(arrange == TRUE){
    idx = order(rowSums(E),decreasing = TRUE)
    P = P[,idx]
    E = E[idx,]
  }
  
  Output <- list()
  Output$P <-  P
  Output$E <-  E
  Output$gkl <- div[best]
  Output$results <- reslist
  
  return(Output)
}

#####################################################
## Cosine similarity between a vector x and y
#####################################################
similarity <- function(x,y){
  dot.prod <- sum(x*y) 
  norm.x <- sqrt(sum(x^2))
  norm.y <- sqrt(sum(y^2))
  frac <- dot.prod / (norm.x * norm.y)
  return(as.numeric(frac))
}

#############################################################
## similarity between the rows of H1 and H2,
## where the rows that are closest are matched
##############################################################
cos.sim <- function(H1,H2){
  K <- nrow(H1)
  d <- numeric(K)
  m <- numeric(K)
  dist <- sapply(1:K, function(y) sapply(1:K,function(x) similarity(H1[x,],H2[y,])))
  dist <- as.matrix(dist)
  residual <- 0
  for(s in 1:K){
    max.dist <- max(dist)
    remove = which(dist == max.dist, arr.ind = TRUE)
    dist[remove[1,1],] <- 0
    dist[,remove[1,2]] <- 0
    d[remove[1,1]] <- max.dist
    m[remove[1,1]] <- remove[1,2]
    residual = residual + max.dist
  }
  
  Output <- list()
  Output$total <- residual/K # The average cosine similarity
  Output$sig   <- d          # Individual cosine similarity numbered after signatures in H1
  Output$match <- m
  return(Output)
}
###############################################################################
## Change because of variance. Illustrated by poisson parametric bootstrapping
##
## Input :  P            - resulting signatures
##          E            - resulting exposures
##          iter         - number of iterations in bootstrap
##          same.init    - logical to decide if the initialization are the same
##                         across the bootstrap samples.
##
## Output: Presults - 'iter' results of P stacked in a matrix ('iter'*ncol(P) x nrow(P))
##         Eresults - 'iter' results of E stacked in a matrix ('iter'*nrow(E) x ncol(E))
##         
###############################################################################
boot_pois = function(P,E, iter, same.init = FALSE){
  N = ncol(P)
  K = nrow(P)
  G = ncol(E)
  
  P.probs = matrix(0, nrow = iter*N, ncol = K)
  E.expos = matrix(0, nrow = iter*N, ncol = G)
  
  PE = P%*%E
  no.seed = c(1,200,400,600,800,1000,1200,1400,1600,1800)

  for(i in 1:iter){
    set.seed(i)
    M.sim = matrix(rpois(K*G, lambda = PE), nrow = K) # compute random sample
    res = NMFPoisEMsquarem(M.sim, N, seed = ifelse(same.init,no.seed,i*c(1,200,400,600,800,1000,1200,1400,1600,1800)), tol = 1e-3)
    dist = cos.sim(t(P),t(res$P))
    P.probs[(1+N*(i-1)):(N+N*(i-1)),] = t(res$P[,dist$match])
    E.expos[(1+N*(i-1)):(N+N*(i-1)),] = res$E[dist$match,]
  }
  
  Output = list()
  
  Output$Presults = P.probs
  Output$Eresults = E.expos
  return(Output)
}

#################################################################
## Finding SVD solutions alpha from samples of matrices P and E 
## from parametric bootstrapping
#################################################################
poissonToSVD = function(P, E, samples = 1000, same.init = FALSE){
  N = ncol(P)
  D = P%*%E
  res = boot_pois(P = P,E = E, iter = samples, same.init = same.init)
  svd.D = svd(t(D), nu = N, nv = N)
  svdV = svd.D$v
  svdU = svd.D$u
  Psample = res$P_lastCheckResults
  Esample = res$E_lastCheckResults
  PT.points = matrix(0, nrow = nrow(Psample), ncol = N-1)
  TE.points = matrix(0, nrow = nrow(Esample), ncol = N-1)
  for(i in 1:samples){
    p = Psample[(i*N-(N-1)):(i*N),]
    Tmat = p%*%svdV  # find Tmat
    Tmat = Tmat/Tmat[,1]
    PT.points[(i*3-2):(i*3),] = Tmat[,c(2,3)]
    e = Esample[(i*N-(N-1)):(i*N),]
    Tmat = e%*%svdU  # find Tmat
    Tmat = Tmat/Tmat[,1]
    TE.points[(i*3-2):(i*3),] = Tmat[,c(2,3)]
  }
  Output = list()
  Output$PT.points = PT.points
  Output$TE.points = TE.points
  return(Output)
}

