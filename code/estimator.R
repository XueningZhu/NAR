
################################### Estimator Function for NAR model ##########################################

betaOLS<-function(Ymat, W, Z)                                                                                  ### OLS estimation for theta: eq (2.8)
{
  Ymat1 = W%*%Ymat                                                                                             ### obtain WY
  Time = ncol(Ymat)-1                                                                                          ### Time doesn't count Y0
  X = cbind(rep(1, nrow(Ymat)*Time),                                                                           ### the intercept
            as.vector(Ymat1[,-ncol(Ymat)]),                                                                    ### WY_{t-1}
            as.vector(Ymat[,-ncol(Ymat)]),                                                                     ### Y_{t-1}
            do.call("rbind", rep(list(Z), Time)))                                                              ### nodal covariates
  invXX = solve(crossprod(X))                                                                                  ### {t(X)X}^{-1}
  Yvec = as.vector(Ymat[,-1])                                                                                  ### the response vector
  thetaEst = invXX%*%colSums(X*Yvec)                                                                           ### estimation equation (2.8)
  sigmaHat2 = mean((Yvec - X%*%thetaEst)^2)                                                                    ### estimation for hat sigma^2
  covHat = invXX*sigmaHat2                                                                                     ### covariance for hat theta
  return(list(theta = thetaEst,
              covHat = covHat, 
              sigmaHat = sqrt(sigmaHat2)))                                                                     ### return the result
}