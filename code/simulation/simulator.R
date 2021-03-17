
library(MASS)
library(Matrix)
library(poweRlaw)                                                                                     ### for power-law distribution

######################################## parameter setup #############################################

path = "F:/networkTS/netVAR_package/data/rda/"                                                       ### set the default path
Nrep = 1000                                                                                           ### number of replications
NSize = c(100, 500, 1000)                                                                             ### network sizes
BetaList = list(c(0.3, 0, 0.5),
                c(0, 0.1, -0.2),
                c(0.3, -0.1, 0.5))

TrueParaList = list(c(10, 30, 100),                                                                   ### setup for TSize (Dyad), block size (Block), and alpha (power-law)
                    c(5, 10, 20),
                    c(1.2, 2, 3))
gamma0 = c(-0.5,0.3,0.8,0,0)                                                                          ### true parameter for gamma
ZSigma = 0.5^abs(outer(1:5,1:5,"-"))                                                                  ### Sigma_z (covariance for Z)


######################################## useful functions #############################################

mvrnorm.Autoreg <- function(n, p, rho = 0.5)
{
  X=rnorm(n*p,0,1);X=array(X,c(n,p));    	                                                             ### autoregressive covariance structure
  for (j in 2:p) {X[,j]=rho*X[,j-1]+sqrt(1-rho^2)*rnorm(n,0,1)}
  return(X)
}



getY<-function(Y0, Beta0, Beta, W, sig = 1, Time = 10)                                                 ### function for simulationg Y series given Y0, beta0, beta = (beta1, beta2), W, sigma, and T
{
  Ymat = matrix(0, nrow = length(Beta0), ncol = Time+1)                                                ### use Ymat to store the simulated data
  Ymat[,1] = as.vector(Y0)                                                                             ### the first column of Ymat is assigned Y0
  for (i in 1:Time)                 
  {
    Ymat[,i+1] = as.vector(Beta0 + Beta[1]*W%*%Ymat[,i] + Beta[2]*Ymat[,i] 
                           + rnorm(length(Beta0), sd = sig))                                           ### follow the NAR model to simulate Y time series
  }
  return(Ymat)
}

getMu<-function(G,Beta0)                                                                               ### given G and Beta0, calculte Mu (stationary mean)
{
  mu = solve(Diagonal(length(Beta0))-G)%*%Beta0                                                        ### See Thoerem 1, stationary mean
  return(mu)
}

getGamma0Approx<-function(G, sigma)                                                                    ### given G and sigma, calculate approximated Gamma(0) using Taylor expansion
{
  Gamma0 = (Diagonal(nrow(G)) + tcrossprod(G) + tcrossprod(G%*%G) +tcrossprod(G%*%G%*%G))*sigma^2       ### calculate approximated Gamma(0) using Taylor expansion
  return(Gamma0)
}

simu.data<-function(W, beta0, Beta, Time = Time, G = G, Z = Z, sig = 1)
{
  mu = getMu(G, beta0 + Z%*%gamma0)
  Gamma0 = getGamma0Approx(G, sigma = sig)
  Y0 = mvrnorm(n = 1, mu = mu, Sigma = Gamma0)
  
  Ymat = getY(Y0 = Y0,#rep(0,N), 
              beta0+ Z%*%gamma0, Beta, W, Time = Time)
  return(Ymat)
}



getDyadW<-function(N, N1 = 10, delta = 1.2, normalize = T)                                             ### simulate Dyad network W; N1: number of mutual pairs 2N1*N, delta: N((0,1)) = N((1,0)) = 0.5*N^{delta}
{
  A = matrix(0, nrow = N, ncol = N)                                                                    ### use A to store network structure
  
  ######################################### mutual follower ###########################################
  ind = which(upper.tri(A), arr.ind = T)                                                               ### record all the index of upper triangular matrix
  indM = ind[sample(1:nrow(ind), N*N1),]                                                               ### sample N*N1 as mutual pairs in the upper triangular matrix
  A[indM] = 1                                                                                          ### the corresponding links are set to be 1
  A[indM[,2:1]] = 1                                                                                    ### the following matrix is set to be symmetric, as a result, mutual pairs are 2N1*N
  
  ######################################### single relationship #######################################
  ind1 = which(A==0&upper.tri(A), arr.ind = T)                                                         ### record all the zero pairs in the upper triangular matrix
  indS = ind1[sample(1:nrow(ind1), N^delta),]                                                          ### choose N^delta index as single relations
  tmp = sample(1:nrow(indS), floor(N^delta/2))                                                         ### randomly choose 0.5*N^delta to be in the lower triangular matrix
  indS[tmp,] = indS[tmp, 2:1]                                                                          ### change the corresponding index to be inverse
  A[indS] = 1                                                                                          ### the single following relation is set to be 
  diag(A) = 0                                                                                          ### aii = 0
  A = as(A, "dgCMatrix")
  if (!normalize)
    return(A)
  W = A/rowSums(A)                                                                                     ### W is row-normalized
  return(W)
}


getBlockW<-function(N, Nblock, normalize = T)                                                          ### get block network
{
  if (N%%Nblock==0){                                                                                   ### if N mod Nblock is integer
    isDiagList = rep(list(matrix(1, nrow = N/Nblock, ncol = N/Nblock)), Nblock)                        ### obtain the diagnal block list
    mList = rep(list(matrix(rbinom((N/Nblock)^2, size = 1, prob = 0.3*N^{-0.3}),                       ### generate following relations within the blocks
                            nrow = N/Nblock, ncol = N/Nblock)), Nblock)
  }
  else
  {
    isDiagList = rep(list(matrix(1, nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)        ### if N mod Nblock is not integer
    isDiagList[[length(Nblock)]] = matrix(1, nrow = N%%Nblock, ncol = N%%Nblock)
    
    mList = rep(list(matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.3*N^{-0.3}),                  ### generate following relations within the blocks
                            nrow = floor(N/Nblock), ncol = floor(N/Nblock))), Nblock-1)
    mList[[Nblock]] = matrix(rbinom(floor(N/Nblock)^2, size = 1, prob = 0.3*N^{-0.3}),                 ### generate following relations within the blocks
                             nrow = floor(N/Nblock), ncol = floor(N/Nblock))
  }
  isDiag = bdiag(isDiagList)                                                                           ### combine the blocks in matrix
  offDiag = which(isDiag == 0, arr.ind = T)                                                            ### to calculate the index of the off digonal indexes
  mList = lapply(mList, function(M){
    ind = which(rowSums(M)==0)
    if (length(ind)>0)
      M[cbind(ind, sample(1:nrow(M), length(ind)))] = 1
    return(M)
  })
  bA = bdiag(mList)
  bA[offDiag] = rbinom(nrow(offDiag), size = 1, prob = 0.3/N)                                          ### people between blocks have 0.3 prob to follow
  bA = as.matrix(bA)
  upperInd = which(upper.tri(bA), arr.ind = T)
  
  ################ transform bA to be a symmetric matrix ##############################################
  bA[upperInd[,2:1]] = bA[upper.tri(bA)]
  diag(bA) = 0
  
  
  ind = which(rowSums(bA)==0)                                                                          ### in case some row sums are zero
  for (i in ind)
  {
    bA[i, sample(setdiff(1:N,i), 3)] = 1                                                               ### for those node, randomly select 3 followees
  }
  bA = as(bA, "dgCMatrix")
  if (!normalize)
    return(bA)
  W = bA/rowSums(bA)                                                                                   ### row normalize bA
  return(as.matrix(W))
}



getPowerLawW<-function(N, alpha, normalize = T)                                                        ### get power-law network W
{
  Nfollowers = rpldis(N, 1, alpha)                                                                     ### generate N random numbers following power-law(1, alpha): k1-kN
  A = sapply(Nfollowers, function(n) {                                                                 ### for node i, randomly select ki nodes to follow it
    vec = rep(0, N)
    vec[sample(1:N, min(n,N))] = 1
    return(vec)
  })
  diag(A) = 0
  ind = which(rowSums(A)==0)                                                                           ### in case some row sums are zero
  for (i in ind)
  {
    A[i, sample(setdiff(1:N,i), 3)] = 1                                                                ### for those node, randomly select 3 followees
  }
  A = as(A, "dgCMatrix")
  if (!normalize)
    return(A)
  W = A/rowSums(A)
  return(W)
}
