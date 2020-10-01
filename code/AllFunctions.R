library("SQUAREM")
library("Rcpp")
library("ggplot2")
library("ggpubr")

## C++ implementation of SFS function
sourceCpp("code/sampleSFS.cpp")

#' Finding SVD representation from solutions of matrices P and E
#'
#' @param Presults Matrix of results of P transposed stacked on top of each other. Dimension is (N*results x nrow(P)).
#' @param Eresults Matrix of results of E stacked on top of each other. Dimension is (N*results x ncol(E))
#' @param Mfit The initial factorization of P and E to use as a reference for the eigenvectors. 
#'  Default is the factorization of the first matrix in Presults and Eresults.
#' @param N The rank of the factorization
#' 
#' @return The SVD representation of the set of feasible solutions
#' \itemize{
#'  \item P.points - Matrix of P results as SVD solution (results x (N-1)).
#'  \item E.points - Matrix of E results as SVD solution (results x (N-1)).
#'  \item plotP - Plot of P.points
#'  \item plotE - Plot of E.points
#'  }
samplesToSVD = function(Presults, Eresults, N, Mfit = t(Presults[1:N,])%*%Eresults[1:N,]){
  svd.Mfit = svd(t(Mfit), nu = N, nv = N)
  svdV = svd.Mfit$v
  svdU = svd.Mfit$u
  P.points = matrix(0, nrow = nrow(Presults), ncol = N-1)
  E.points = matrix(0, nrow = nrow(Eresults), ncol = N-1)
  for(i in 1:(nrow(Presults)/N)){
    p = Presults[(i*N-(N-1)):(i*N),]
    Tmat = p%*%svdV  # find Tmat
    Tmat = Tmat/Tmat[,1]
    P.points[(i*3-2):(i*3),] = Tmat[,c(2,3)]
    e = Eresults[(i*N-(N-1)):(i*N),]
    Tmat = e%*%svdU  # find Tmat
    Tmat = Tmat/Tmat[,1]
    E.points[(i*3-2):(i*3),] = Tmat[,c(2,3)]
  }
  
  # plot svd for P 
  dat.P = data.frame(x = P.points[,1], y = P.points[,2])
  dat.E = data.frame(x = E.points[,1], y = E.points[,2])
  
  gP = ggplot(dat.P, aes(x = x, y = y))+
    geom_point(size = 1, alpha = 0.2)+
    labs(x = expression(alpha[1]), y = expression(alpha[2]))+
    theme_bw()+
    theme(legend.position = "none")
  
  gE = ggplot(dat.E, aes(x = x, y = y))+
    geom_point(size = 1, alpha = 0.2)+
    labs(x = expression(alpha[1]), y = expression(alpha[2]))+
    theme_bw()+
    theme(legend.position = "none")
  
  
  Output = list()
  Output$P.points = P.points
  Output$E.points = E.points
  Output$plotP = gP
  Output$plotE = gE
  return(Output)
}


#' Function for calculating the generalised Kullback Leibler divergence
#' 
#' Internal function used in NMFPois
#' 
#' @param y Observation
#' @param mu Estimate
#' 
#' @return Generalized Kullback-Leibler
gkl.dev <- function(y, mu){
  r <- mu
  p <- which(y > 0)
  r[p] <- (y * (log(y)- log(mu)) - y + mu)[p]
  return(sum(r))
}

#' @title Non-negative matrix factorization algorithm for Poisson data
#' 
#' Factorizing M into two matrices P and E of
#' dimension ncol(M) x N and N x nrow(M) with the acceleration of SQUAREM.
#' The objective function is the generalized Kullback-Leibler divergence(GKLD).
#' 
#' @param M Non-negative data matrix of size
#' @param N Small dimension of the two new matrices
#' @param seed  Vector of random seeds to initialize the matrices
#' @param arrange Arranging columns in P and rows in E after largest row sums of E 
#' @param tol Maximum change of P and E when stopping 
#'
#' @return A list of the matrices derived by the factorization and the corresponding generalized Kullback-Leibler
#'  \itemize{
#'  \item P   - Non-negative matrix of dimension ncol(V) x K, with columns summing to one
#'  \item E   - Non-negative matrix of dimension K x nrow(V), where rows sum to one 
#'  \item gkl - Smallest Value of the Generalized Kullback-Leibler
#'  }
NMFPois = function(M,N,seed, arrange = TRUE, tol = 1e-5){
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
    P <- P %*% diag(1/rowSums(E))   
    
    PE <- P%*%E
    E <- E * (t(P) %*% (M/PE))      # update of exposures
    E <- diag(1/colSums(P)) %*% E
    
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
    
    P <- matrix(stats::runif(K*N), nrow = K, ncol = N)  # Initialize P
    E <- matrix(stats::runif(N*G), nrow = N, ncol = G)  # Initialize E
    
    init = log(c(as.vector(P),as.vector(E)))
    sres = squarem(init, fixptfn = poisson.em, objfn = gklobj, control = list(tol = tol))
    
    P = matrix(exp(sres$par[1:(K*N)]), nrow = K, ncol = N)
    E = matrix(exp(sres$par[-c(1:(K*N))]), nrow = N, ncol = G)
    E = diag(colSums(P)) %*% E # normalizing 
    P = P %*% diag(1/colSums(P))
    
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


#' Cosine similarity between two vectors x and y
similarity <- function(x,y){
  dot.prod <- sum(x*y) 
  norm.x <- sqrt(sum(x^2))
  norm.y <- sqrt(sum(y^2))
  frac <- dot.prod / (norm.x * norm.y)
  return(as.numeric(frac))
}


#' similarity and match between the rows of P1 and P2, using hierarchical clustering of cosine similarity
cos.sim <- function(P1,P2){
  K <- nrow(P1)
  d <- numeric(K)
  m <- numeric(K)
  dist <- sapply(1:K, function(y) sapply(1:K,function(x) similarity(P1[x,],P2[y,])))
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


#' Change of P and E because of poisson parametric bootstrapping
#'
#' @param P resulting factor P from NMF
#' @param E resulting factor E from NMF
#' @param B number of bootstrap samples
#' @param same.init logical to decide if the initialization are the same across the bootstrap samples.
#'
#' @return The results of P and E from the parametric bootstrapping
#' \itemize{
#'  \item Presults - 'iter' results of P stacked in a matrix (B*ncol(P) x nrow(P))
#'  \item Eresults - 'iter' results of E stacked in a matrix ('iter'*nrow(E) x ncol(E))
#'  }
boot_pois = function(P,E, iter, same.init = FALSE){
  N = ncol(P)
  K = nrow(P)
  G = ncol(E)
  
  P.probs = matrix(0, nrow = iter*N, ncol = K)
  E.expos = matrix(0, nrow = iter*N, ncol = G)
  
  PE = P%*%E
  no.seed = sample(1:1000,10)

  for(i in 1:iter){
    set.seed(i)
    M.sim = matrix(rpois(K*G, lambda = PE), nrow = K) # compute random sample
    res = NMFPois(M.sim, N, seed = ifelse(same.init,no.seed,sample(1:1000,10)), tol = 1e-3)
    dist = cos.sim(t(P),t(res$P))
    P.probs[(1+N*(i-1)):(N+N*(i-1)),] = t(res$P[,dist$match])
    E.expos[(1+N*(i-1)):(N+N*(i-1)),] = res$E[dist$match,]
  }
  
  Output = list()
  
  Output$Presults = P.probs
  Output$Eresults = E.expos
  return(Output)
}
