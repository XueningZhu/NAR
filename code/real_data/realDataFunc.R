
library(plyr)
library(Matrix)
library(car)
specify_decimal = function(x, k) gsub('\\s+','',format(round(x, k), nsmall=k))        ### format function for keeping 2 digits after
lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")                        ### set the time format

setM = function (D, vid, score)                                                       ### transform to sparse matrix
{
  M = Matrix(data = 0, nrow = length(vid), ncol = length(vid))
  uidx = match(D[[1]], vid, nomatch = 0)
  sidx = match(D[[2]], vid, nomatch = 0)
  i = uidx > 0 & sidx > 0
  M = as(M, "dgTMatrix")
  M@i = as.integer(uidx[i] - 1)
  M@j = as.integer(sidx[i] - 1)
  M@x = as.numeric(score[i])
  M = as(M, "dgCMatrix")
  #colnames(M) = vid
  return(M)
}


getNetStruc<-function(file)                                                            ### function for get the network structure: Adjency mat
{
  A = read.csv(file)                                                                   ### read the network structure files (edges*2)
  IDs = sort(union(A$id1, A$id2))                                                  ### keep the nodes who not only follow others but also followed by others                                                  
  A = A[is.element(A$id1, IDs)&is.element(A$id2, IDs),]                                ### only keep the edges related to IDs
  A = A[order(A$id1, A$id2),]                                                          ### reorder A
  Amat = setM(D = A, vid = IDs, score = rep(1, nrow(A)))                               ### obtain the sparse adjency matrix
  return(list(IDs = IDs, Amat = Amat))                                                 ### return the IDs and adjency matrix
}



getRegDf<-function(logYmatNchar,  W , Nodal)                                                 ### obtain regression dataset
{
  Time = ncol(logYmatNchar); tt = 1:Time
  logYmatNchar1 = W%*%logYmatNchar
  Zmat = do.call(rbind, rep(list(Nodal), Time-1))
  
  X = cbind(network = as.vector(logYmatNchar1[,tt[-length(tt)]]),                                        ### obtain all variables in X
            momentum = as.vector(logYmatNchar[,tt[-length(tt)]]),
            Zmat)
  Yvec = as.vector(logYmatNchar[,tt[-1]])                                                                ### the response Y vector
  regDf = data.frame(Y = Yvec, X)
  return(regDf)
}

