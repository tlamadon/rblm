#' Model creation
#'
#' Creates a 2 period mini model. Can be set to be linear or not, with serial correlation or not
#' @export
m2.mini.new <- function(nf,linear=FALSE,serial=F,fixb=F) {
  m = list()
  m$nf = nf

  # generating intercepts and interaction terms
  m$A1 = 0.5*rnorm(nf)
  m$B1 = exp(rnorm(nf)/3)
  m$A2 = 0.5*rnorm(nf)
  m$B2 = exp(rnorm(nf)/3)

  # generate the variance of eps
  m$eps1_sd  = exp(rnorm(nf)/2)
  m$eps2_sd  = m$eps1_sd
  m$eps_cor  = runif(1)

  # generate the E(alpha|l) and sd
  m$Em   = sort(rnorm(nf))
  m$Esd  = exp(rnorm(nf)/2)

  # set normalization
  m$B2 = m$B2/m$B2[1]
  m$B1 = m$B1/m$B2[1]
  m$A1[1] = 0

  if (linear) {
    m$B1[] = 1
    m$B2[] = 1
  }

  if (serial==F) {
    m$eps_cor = 0
  }

  if (fixb) {
    m$B1=m$B2
  }

  # sort
  I = order(m$A1 + m$B1*m$Em)
  m$A1 = m$A1[I]
  m$A2 = m$A2[I]
  m$B1 = m$B1[I]
  m$B1 = m$B2[I]
  m$Em = m$Em[I]

  # generate the E(alpha|l,l') and sd
  m$EEm  = 0.8*spread(m$Em,2,nf) +  0.2*spread(m$Em,1,nf) + 0.3*array(rnorm(nf^2),c(nf,nf))
  m$EEsd = exp(array(rnorm(nf^2),c(nf,nf))/2)

  return(m)
}

lm.wfitc <- function(XX,YY,rw,C1,C0,meq) {

  XXw      = diag(rw) %*% XX
  Dq       = t(XXw) %*% XX
  dq       = t(YY %*% XXw)

  # do quadprod
  fit      = solve.QP(Dq,dq,t(C1),C0,meq)
  return(fit)
}

#' simulate a dataset according to the model
#' @export
m2.mini.simulate <- function(model) {

  # compute decomposition
  stayer_share = sum(model$Ns)/(sum(model$Ns)+sum(model$Nm))
  model$vdec   = m2.mini.vdec(model,1e6,stayer_share,"y1")

  sdata = m2.mini.simulate.stayers(model,model$Ns)
  jdata = m2.mini.simulate.movers(model,model$Nm)

  # randomly assign firm IDs
  sdata <- sdata[,f1:=paste("F",j1 + model$nf*(sample.int(.N/100,.N,replace=T)-1),sep=""),j1]
  sdata <- sdata[,j1b:=j1]
  sdata <- sdata[,j1true := j1]
  jdata <- jdata[,j1true := j1][,j2true := j2]
  jdata <- jdata[,j1c:=j1]
  jdata <- jdata[,f1:=sample( unique(sdata[j1b==j1c,f1]) ,.N,replace=T),j1c]
  jdata <- jdata[,j2c:=j2]
  jdata <- jdata[,f2:=sample( unique(sdata[j1b==j2c,f1])  ,.N,replace=T),j2c]
  jdata$j2c=NULL
  jdata$j1c=NULL
  sdata$j1b=NULL
  sdata$x=1

  # combine the movers and stayers, ad stands for all data:
  ad = list(sdata=sdata,jdata=jdata)

  return(ad)
}


#' simulate movers according to the model
#' @export
m2.mini.simulate.movers <- function(model,NNm) {

  J1 = array(0,sum(NNm))
  J2 = array(0,sum(NNm))
  Y1 = array(0,sum(NNm))
  Y2 = array(0,sum(NNm))
  e1 = array(0,sum(NNm))
  e2 = array(0,sum(NNm))
  K  = array(0,sum(NNm))

  A1  = model$A1
  B1  = model$B1
  A2  = model$A2
  B2  = model$B2
  EEm = model$EEm
  EEsd= model$EEsd
  eps_sd1=model$eps1_sd
  eps_sd2=model$eps2_sd
  eps_cor= 0 # model$eps_cor, not relevant anymore
  nf  = model$nf

  i =1
  for (l1 in 1:nf) for (l2 in 1:nf) {
    if (NNm[l1,l2]==0) next;
    I = i:(i+NNm[l1,l2]-1)
    ni = length(I)
    jj = l1 + nf*(l2 -1)
    J1[I] = l1
    J2[I] = l2

    # draw alpha
    K[I] = EEm[l1,l2] + EEsd[l1,l2]*rnorm(ni)

    # draw epsilon1 and epsilon2 corolated
    eps1 = eps_sd1[l1] * rnorm(ni)
    eps2 = eps_cor/eps_sd1[l1]*eps_cor*eps1    +    sqrt(  (1-eps_cor^2) *  eps_sd2[l2]^2 ) * rnorm(ni)

    Y1[I]  = A1[l1] + B1[l1]*K[I] + eps1
    Y2[I]  = A2[l2] + B2[l2]*K[I] + eps2
    e1[I]  = eps1
    e2[I]  = eps2

    i = i + NNm[l1,l2]
  }

  jdatae = data.table(alpha=K,y1=Y1,y2=Y2,j1=J1,j2=J2,e1=e1,e2=e2)
  return(jdatae)
}

#' simulate movers according to the model
#' @export
m2.mini.simulate.stayers <- function(model,NNs) {

  J1 = array(0,sum(NNs))
  Y1 = array(0,sum(NNs))
  Y2 = array(0,sum(NNs))
  e1 = array(0,sum(NNs))
  K  = array(0,sum(NNs))

  A1  = model$A1
  B1  = model$B1
  A2  = model$A2
  B2  = model$B2
  Em  = model$Em
  Esd = model$Esd
  eps_sd1=model$eps1_sd
  eps_sd2=model$eps2_sd

  nf  = model$nf

  i =1
  for (l1 in 1:nf) {
    I = i:(i+NNs[l1]-1)
    ni = length(I)
    J1[I] = l1

    # draw alpha
    K[I] = Em[l1] + Esd[l1]*rnorm(ni)

    # draw Y2, Y3
    e1 = rnorm(ni) * eps_sd1[l1]
    e2 = rnorm(ni) * eps_sd2[l1]
    Y1[I]  = A1[l1] + B1[l1]*K[I] + e1
    Y2[I]  = A2[l1] + B2[l1]*K[I] + e2

    i = i + NNs[l1]
  }

  sdatae = data.table(alpha=K,y1=Y1,y2=Y2,j1=J1,j2=J1)
  return(sdatae)
}

#' this function uses LIML to estimate a model with a given normalization
#' and we want to do it in a non-stationary way
m2.mini.liml.int <- function(Y1,Y2,J1,J2,norm=1) {

  L = max(J1)
  N = length(Y1)

  # ----- prepare the different regressors ----- #
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  D2 = array(0,c(N,L))
  I = (1:N) + N*(J2-1)
  D2[I]=1

  # contrust full matrix
  X1 = cbind(D1*spread(Y1,2,L),D2*spread(Y2,2,L))  # construct the matrix for the interaction terms
  X2 = cbind(D1,D2)                                # construct the matrix for the intercepts

  # set the normalizations
  Y  = -X1[,norm]
  X1 =  X1[,setdiff(1:(2*L),norm)]
  X2 =  X2[,1:(2*L-1)]

  # construct the Z matrix of instruemts (j1,j2), dummies of (j1,j2) interactions
  Z = array(0,c(N,L^2))
  I = (1:N) + N*(J1-1 + L*(J2-1))
  Z[I]=1

  # make matrices sparse
  Z  = Matrix(Z,sparse=T)
  R  = cbind(Y,X1)
  X2 = Matrix(X2,sparse=T)
  X1 = Matrix(X1,sparse=T)

  # LIML
  Wz = t(R) %*% R  -  ( t(R) %*% Z)  %*% solve( t(Z)  %*% Z  ) %*% ( t(Z) %*% R)
  Wx = t(R) %*% R  -  ( t(R) %*% X2) %*% solve( t(X2) %*% X2 ) %*% ( t(X2) %*% R)

  WW = Wx %*% solve(Wz)
  #WW = Wx %*% ginv(as.matrix(Wz))
  lambdas = eigen(WW)$values
  lambda  = min(lambdas)

  XX = cBind(X1,X2)
  RR = (1-lambda)*t(XX) %*% XX + lambda * ( t(XX) %*% Z)  %*% solve( t(Z)  %*% Z  ) %*% ( t(Z) %*% XX)
  RY = (1-lambda)*t(XX) %*% Y  + lambda * ( t(XX) %*% Z)  %*% solve( t(Z)  %*% Z  ) %*% ( t(Z) %*% Y)

  # --------- extract the results ---------- #
  b_liml = as.numeric(solve(RR) %*% RY)
  tau = rep(0,L)
  for (i in 1:(L-1)) { tau[i+ (i>=norm)] = b_liml[i ]}
  tau[norm]=1
  B1 = 1/tau
  B2 = - 1/b_liml[L:(2*L-1)]
  A1 = - b_liml[(2*L):(3*L-1)] * B1
  A2 = c(b_liml[(3*L):(4*L-2)],0) * B2

  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=b_liml))
}

#' this function uses LIML to estimate a model with a given normalization
#' and we want to do it in a non-stationary way.
m2.mini.liml.int2 <- function(Y1,Y2,J1,J2,norm=1,coarse=0,tik=0) {

  L = max(J1)
  N = length(Y1)

  # ----- prepare the different regressors ----- #
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  D2 = array(0,c(N,L))
  I = (1:N) + N*(J2-1)
  D2[I]=1

  # construct the Z matrix of instruemts (j1,j2)
  if (coarse==0) {
    Z = array(0,c(N,L^2))
    I = (1:N) + N*(J1-1 + L*(J2-1))
    Z[I]=1
  } else {
    J1c = ceiling(J1/coarse)
    J2c = ceiling(J2/coarse)
    Lc = max(J1c)

    Z1 = array(0,c(N,L*Lc))
    I1 = (1:N) + N*(J1-1 + L*(J2c-1))
    Z1[I1]=1

    Z2 = array(0,c(N,L*Lc))
    I2 = (1:N) + N*(J2-1 + L*(J1c-1))
    Z2[I2]=1
    Z = cbind(Z1,Z2)
  }

  # create regressors, make sparse
  Z  = Matrix(Z,sparse=T)
  X  = Matrix(cbind(D2*spread(Y2,2,L), -D1*spread(Y1,2,L),D1[,2:L],-D2)  ,sparse=T)

  # as.numeric(t(Z) %*% X %*% c(model$B2,model$B1,model$A1[2:L],model$A2)
  t = Matrix:::t

  # ------ LIML procedure ---------- #
  Wz  = t(X) %*% Z %*% solve(t(Z) %*% Z + tik * diag(L^2)) %*% t(Z) %*% X
  Wx  = t(X) %*% X
  WW  = solve(Wx) %*% Wz
  dec = eigen(WW)
  b_liml = as.numeric(dec$vectors[,4*L-1])

  # --------- extract the results ---------- #
  b_liml = b_liml/b_liml[norm]
  B2 = 1/b_liml[1:L]
  B1 = 1/b_liml[(L+1):(2*L)]
  A1 = c(0,b_liml[(2*L+1):(3*L-1)]) * B1
  A2 = b_liml[(3*L):(4*L-1)] * B2

  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=b_liml))
}

#' this function uses LIML to estimate a model with a given normalization
#' and we want to do it in a non-stationary way.
m2.mini.liml.int.fixb <- function(Y1,Y2,J1,J2,norm=1,prof=F) {

  L = max(J1)
  N = length(Y1)

  # ----- prepare the different regressors ----- #
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  D2 = array(0,c(N,L))
  I = (1:N) + N*(J2-1)
  D2[I]=1

  # construct the Z matrix of instruemts (j1,j2)
  Z = array(0,c(N,L^2))
  I = (1:N) + N*(J1-1 + L*(J2-1))
  Z[I]=1

  # remove empty interactions
  I = colSums(Z)!=0
  Z=Z[,I]

  # create regressors, make sparse
  Z  = Matrix(Z,sparse=T)
  X  = Matrix(cbind(D2*spread(Y2,2,L) - D1*spread(Y1,2,L),D1[,2:L],-D2)  ,sparse=T)

  # profiling
  if (prof==TRUE) {
    D2Y2 = cbind(D2*spread(Y2,2,L),-D2[,2:L])
    D1Y1 = cbind(-D1*spread(Y1,2,L),D1)
    X = D2Y2 - D1Y1 %*% (ginv(D1Y1) %*% D2Y2)
    X = Matrix(X,sparse=T)
  }

  # as.numeric(t(Z) %*% X %*% c(model$B2,model$B1,model$A1[2:L],model$A2)
  t = Matrix:::t

  # ------ LIML procedure ---------- #
  Wz  = t(X) %*% Z %*% solve( t(Z) %*% Z ) %*% t(Z) %*% X
  Wx  = t(X) %*% X
  WW  = solve(Wx) %*% Wz
  dec = eigen(WW)
  b_liml = Re(as.numeric(dec$vectors[,3*L-1]))

  # --------- extract the results ---------- #
  b_liml = b_liml/b_liml[norm]
  B2 = 1/b_liml[1:L]
  B1 = B2
  A1 = c(0,b_liml[(L+1):(2*L-1)]) * B1
  A2 = b_liml[(2*L):(3*L-1)] * B2

  if (any(is.na(A1*A2))) warnings("A1 or A2 contains NA values");
  if (any(is.na(B1*B2))) warnings("A1 or A2 contains NA values");

  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=b_liml))
}

#' this function uses LIML to estimate a model with a given normalization.
#' In estimation the interaction terms B are profiled out in period 2 first,
#' then they are used to recover the intercepts in both periods.
m2.mini.liml.int.prof <- function(Y1,Y2,J1,J2,norm=1) {

  L = max(J1)
  N = length(Y1)

  # ----- prepare the different regressors ----- #
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  D2 = array(0,c(N,L))
  I = (1:N) + N*(J2-1)
  D2[I]=1

  # construct the Z matrix of instruemts (j1,j2)
  Z = array(0,c(N,L^2))
  I = (1:N) + N*(J1-1 + L*(J2-1))
  Z[I]=1

  # remove empty interactions
  I = colSums(Z)!=0
  Z=Z[,I]

  # create regressors, make sparse
  Z  = Matrix(Z,sparse=T)

  # profiling, we project on D2Y2 and intercepts
  D2Y2 = cbind(D2*spread(Y2,2,L))
  D1Y1 = cbind(-D1*spread(Y1,2,L),D1[,2:L],-D2)
  X = D2Y2 - D1Y1 %*% (ginv(D1Y1) %*% D2Y2)
  X = Matrix(X,sparse=T)

  # extract the correct transpose function
  t = Matrix:::t

  # ------ LIML procedure ---------- #
  Wz  = t(X) %*% Z %*% solve( t(Z) %*% Z ) %*% t(Z) %*% X
  Wx  = t(X) %*% X
  WW  = solve(Wx) %*% Wz
  dec = eigen(WW)
  b_liml = Re(as.numeric(dec$vectors[,L]))

  # finally extract the intercepts, this is a linear regression
  X = cbind(D1[,2:L],-D2)
  Y = (D1*spread(Y1,2,L)) %*% b_liml - (D2*spread(Y2,2,L))%*%b_liml
  fit2 = lm.fit(X,Y)
  b_liml = c(b_liml,as.numeric(fit2$coefficients))

  # --------- extract the interaction terms ---------- #
  b_liml = b_liml/b_liml[norm]
  B2 = 1/b_liml[1:L]
  B1 = B2
  A1 = c(0,b_liml[(L+1):(2*L-1)]) * B1
  A2 = b_liml[(2*L):(3*L-1)] * B2

  if (any(is.na(A1*A2))) warnings("A1 or A2 contains NA values");
  if (any(is.na(B1*B2))) warnings("A1 or A2 contains NA values");

  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=b_liml))
}

#' this function uses performs an estimation like the LIML but forcing
#' all the Bs=1. This is like an AKM estimator.
m2.mini.linear.int <- function(Y1,Y2,J1,J2,norm=1) {

  L = max(c(J1,J2))
  N = length(Y1)

  # ----- prepare the different regressors ----- #
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  D2 = array(0,c(N,L))
  I = (1:N) + N*(J2-1)
  D2[I]=1

  X = cbind(D1[,2:L],-D2)
  fit = lm.fit(X,Y1-Y2)$coefficients

  # --------- extract the interaction terms ---------- #
  B2 = rep(1,L)
  B1 = rep(1,L)
  A1 = c(0,fit[1:(L-1)])
  A2 = fit[L:(2*L-1)]

  if (any(is.na(A1*A2))) warnings("A1 or A2 contains NA values");
  if (any(is.na(B1*B2))) warnings("B1 or B2 contains NA values");

  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=fit))
}

m2.mini.linear.int.stationary <- function(Y1,Y2,J1,J2,norm=1) {

  L = max(c(J1,J2))
  N = length(Y1)

  # ----- prepare the different regressors ----- #
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  I = (1:N) + N*(J2-1)
  D1[I]=D1[I] - 1

  X = cbind(D1[,2:L])
  fit = lm.fit(X,Y1-Y2)$coefficients

  # --------- extract the interaction terms ---------- #
  B2 = rep(1,L)
  B1 = rep(1,L)
  A1 = c(0,fit[1:(L-1)])
  A2 = A1

  if (any(is.na(A1*A2))) warnings("A1 or A2 contains NA values");
  if (any(is.na(B1*B2))) warnings("B1 or B2 contains NA values");

  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=fit))
}



#' this function uses LIML to estimate a model with a given normalization
#' and we want to do it in a non-stationary way.
m2.mini.liml.int3 <- function(Y1,Y2,J1,J2,norm=1,C=1) {

  L = max(J1)
  N = length(Y1)

  # ----- prepare the different regressors ----- #
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  D2 = array(0,c(N,L))
  I = (1:N) + N*(J2-1)
  D2[I]=1

  # construct the Z matrix of instruemts (j1,j2)
  Z = array(0,c(N,L^2))
  I = (1:N) + N*(J1-1 + L*(J2-1))
  Z[I]=1

  # create regressors, make sparse
  Z  = Matrix(Z,sparse=T)
  X  = Matrix(cbind(D2*spread(Y2,2,L), -D1*spread(Y1,2,L),D1[,2:L],-D2)  ,sparse=T)

  # as.numeric(t(Z) %*% X %*% c(model$B2,model$B1,model$A1[2:L],model$A2)

  # ------ LIML procedure ---------- #
  Pz = Z %*% solve( t(Z) %*% Z ) %*% t(Z)
  diag(Pz)<-0
  Wz  = t(X) %*% Pz %*% X
  Wx  = t(X) %*% X
  WW  = solve(Wx) %*% Wz
  dec = eigen(WW)
  lambda = min(dec$values)

  # we use HFUL of Hausman Newey et Al (2012) QE
  lambda2 = (lambda - (1-lambda)*C/N)/(1 - (1-lambda)*C/N)

  # solve Wx^-1 x Wz v = lambda v where we normalize v[1]=1
  # extract the
  M = WW - lambda2*diag(nrow(WW))
  MX = M[,2:(4*L-1)]
  b_liml = solve(t(MX)%*%MX,-t(MX)%*%M[,1])

  # --------- extract the results ---------- #
  b_liml= c(1,as.numeric(b_liml))
  B2 = 1/b_liml[1:L]
  B1 = 1/b_liml[(L+1):(2*L)]
  A1 = c(0,b_liml[(2*L+1):(3*L-1)]) * B1
  A2 = b_liml[(3*L):(4*L-1)] * B2

  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=b_liml))
}


#' this function uses LIML to estimate a model with a given normalization
#' and we want to do it in a non-stationary way.
m2.mini.nls.int2 <- function(Y1,Y2,J1,J2,norm=1,coarse=0,tik=0) {

  L = max(J1)
  N = length(Y1)

  # ----- prepare the different regressors ----- #
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  D2 = array(0,c(N,L))
  I = (1:N) + N*(J2-1)
  D2[I]=1

  # create regressors, make sparse
  Z  = Matrix(Z,sparse=T)
  X  = Matrix(cbind(D2*spread(Y2,2,L), -D1*spread(Y1,2,L),D1[,2:L],-D2)  ,sparse=T)

  # as.numeric(t(Z) %*% X %*% c(model$B2,model$B1,model$A1[2:L],model$A2)

  # ------ LIML procedure ---------- #
  Wz  = t(X) %*% Z %*% solve(t(Z) %*% Z + tik * diag(L^2)) %*% t(Z) %*% X
  Wx  = t(X) %*% X
  WW  = solve(Wx) %*% Wz
  dec = eigen(WW)
  b_liml = as.numeric(dec$vectors[,4*L-1])
  browser()
  # --------- extract the results ---------- #
  b_liml = b_liml/b_liml[norm]
  B2 = 1/b_liml[1:L]
  B1 = 1/b_liml[(L+1):(2*L)]
  A1 = c(0,b_liml[(2*L+1):(3*L-1)]) * B1
  A2 = b_liml[(3*L):(4*L-1)] * B2

  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=b_liml))
}

# ===== MAIN ESTIMATOR ======

#' estimates interacted model
#'
#' @param method default is "ns" for non-stationary, "fixb" is for
#' stationary, "prof" is for profiling out B in period 2 first and "linear"
#' to run an AKM type estimation
#' @family step2 interacted
#' @export
m2.mini.estimate <- function(jdata,sdata,norm=1,model0=c(),method="ns",withx=FALSE,bigk=0) {

  flog.info("Beginning m2.mini.estimate with method %s",method)

  # --------- use LIML on movers to get A1,B1,A2,B2 ---------- #
  Y1=jdata$y1;Y2=jdata$y2;J1=jdata$j1;J2=jdata$j2;
  nf = max(c(J1,J2))
  if (method=="fixb") {
    rliml = m2.mini.liml.int.fixb(Y1,Y2,J1,J2,norm=norm)
  }  else if (method=="prof") {
    rliml = m2.mini.liml.int.prof(Y1,Y2,J1,J2,norm=norm)
  }  else if (method=="linear") {
    rliml = m2.mini.linear.int(Y1,Y2,J1,J2,norm=norm)
  }  else if (method=="linear.ss") {
    rliml = m2.mini.linear.int.stationary(Y1,Y2,J1,J2,norm=norm)
  } else {
    rliml = m2.mini.liml.int2(Y1,Y2,J1,J2,norm=norm)
  }
  B1 = rliml$B1
  B2 = rliml$B2
  N = length(Y1)

  flog.info("getting worker compositions in movers")
  # ---------  use stayers to get E[alpha|l] ------------- #
  EEm = jdata[, (mean(y1)- rliml$A1[j1])/rliml$B1[j1] + (mean(y2)- rliml$A2[j2])/rliml$B2[j2] ,list(j1,j2)]
  EEm = 0.5*mcast(EEm,"V1","j1","j2",c(nf,nf),0) # what do you think of this 0.5 right here?

  if (!withx) {
    Em  = 0.5*( sdata[, (mean(y1)- rliml$A1[j1])/rliml$B1[j1],j1][order(j1),V1]  +
                  sdata[, (mean(y2)- rliml$A2[j2])/rliml$B2[j2],j2][order(j2),V1])
  } else {
    Em1  = sdata[, (mean(y1)- rliml$A1[j1])/rliml$B1[j1],list(x,j1)]
    Em2  = sdata[, (mean(y2)- rliml$A2[j2])/rliml$B2[j2],list(x,j2)]
    Em   = 0.5*(acast(Em1,j1~x) + acast(Em2,j2~x) )
  }

  # ---------- MOVERS: use covariance restrictions  -------------- #
  # we start by computing Var(Y1), Var(Y2) and Cov(Y1,Y2)
  # we can only get this when there are at least 2 movers in the combination
  flog.info("getting residual wage variances in movers")
  if (bigk==0) {
    setkey(jdata,j1,j2)
    YY1 = c(mcast(jdata[,mvar(y1),list(j1,j2)],"V1","j2","j1",c(nf,nf),0))
    YY2 = c(mcast(jdata[,mcov(y1,y2),list(j1,j2)],"V1","j2","j1",c(nf,nf),0)) #jdata[,mcov(y1,y2),list(j1,j2)][,V1]
    YY3 = c(mcast(jdata[,mvar(y2),list(j1,j2)],"V1","j2","j1",c(nf,nf),0)) #jdata[,mvar(y2),list(j1,j2)][,V1]
    XX1 = array(0,c(nf^2, nf^2 + 2*nf))
    XX2 = array(0,c(nf^2, nf^2 + 2*nf))
    XX3 = array(0,c(nf^2, nf^2 + 2*nf))
    W = c(mcast(jdata[,.N,list(j1,j2)],"N","j2","j1",c(nf,nf),0)) #jdata[,.N,list(j1,j2)][,N]

    for (l1 in 1:nf) for (l2 in 1:nf) {
      ll = l2 + nf*(l1 -1)
      if (W[ll]>0) {
        XX1[ll,ll                 ] = B1[l1]^2
        XX1[ll,nf^2 + l1          ] = 1
        XX2[ll,ll                 ] = B1[l1]*B2[l2]
        XX3[ll,ll                 ] = B2[l2]^2
        XX3[ll,nf^2 + nf + l2     ] = 1
      } else { # force parameter to 0 when there is no info
        XX1[ll,ll                 ] = 0.5 # make sure to force both to 0
        XX3[ll,ll                 ] = 1.5 # make sure to force both to 0
        XX2[ll,ll]                  = 1
        XX1[ll,nf^2 + l1          ] = 1
        XX3[ll,nf^2 + nf + l2     ] = 1
        W[ll]=1e-6 # we use a very small weight, to make bias as small as possible
      }
    }

    Wm = mcast(jdata[,.N,list(j1,j2)],"N","j1","j2",c(nf,nf),0)
    XX = rbind(XX1,XX2,XX3)

    res     = lm.wfitnn( XX, c(YY1,YY2,YY3), c(W,W,W))$solution
    EEsd    = t(sqrt(pmax(array((res)[1:nf^2],c(nf,nf)),0)))
    eps1_sd = sqrt(pmax((res)[(nf^2+1):(nf^2+nf)],0))
    eps2_sd = sqrt(pmax((res)[(nf^2+nf+1):(nf^2+2*nf)],0))
    if (any(is.na(eps1_sd*eps2_sd))) warning("NAs in Var(nu1|k1,k2)")
    if (any(is.na(EEsd))) warning("NAs in Var(alpha|k1,k2)")
  } else {
    # we use a simpler approach here, and we ignore the EEsd and focus on the
    # eps_sd
    YY1 = mcast(jdata[,mvar(y1),list(j1,j2)],"V1","j2","j1",c(nf,nf),0)
    YY2 = mcast(jdata[,mcov(y1,y2),list(j1,j2)],"V1","j2","j1",c(nf,nf),0)
    YY3 = mcast(jdata[,mvar(y2),list(j1,j2)],"V1","j2","j1",c(nf,nf),0)
    Wm   = mcast(jdata[,.N,list(j1,j2)],"N","j2","j1",c(nf,nf),0)
    EEsd = array(0,c(nf,nf))
    eps1_sd = array(0,c(nf))
    eps2_sd = array(0,c(nf))

    for (l1 in 1:nf) {
      v1      =  wt.mean( YY1[l1,] - rliml$B1[l1]/rliml$B2 * YY2[l1,], Wm[l1,] + 1e-10)
      eps1_sd[l1] = sqrt(pmax(v1,0))
      v2      =  wt.mean( YY3[,l1] - rliml$B2[l1]/rliml$B1 * YY2[,l1], Wm[,l1] + 1e-10)
      eps2_sd[l1] = sqrt(pmax(v2,0))
    }

    for (l1 in 1:nf) for (l2 in 1:nf) {
      v1 = 1/(rliml$B1[l1]*rliml$B2[l2]) * YY2[l1,l2]
      EEsd[l1,l2] = sqrt(pmax(v1,0))
    }
  }
  # ---------- STAYERS: use covariance restrictions  -------------- #
  flog.info("getting residual wage variances in stayers")
  # function which computes the variance terms
  getEsd <- function(sdata) {
    setkey(sdata,j1)
    YY1 = sdata[,mvar(y1),list(j1)][,V1] - eps1_sd^2
    YY2 = sdata[,mvar(y2),list(j1)][,V1] - eps2_sd^2
    XX1 = diag(B1^2)
    XX2 = diag(B2^2)
    W = sdata[,.N,list(j1)][,N]

    XX = rbind(XX1,XX2)
    res2 = lm.wfitnn( XX, c(YY1,YY2), c(W,W))$solution
    Esd = sqrt(pmax(res2,0))
    return(Esd)
  }

  if (withx) {
    nx = max(sdata$x)
    Esd = array(0,c(nf,nx))
    for (xl in 1:nx) {
      Esd[,xl] = getEsd(sdata[x==xl])
    }
  } else {
    Esd = getEsd(sdata)
  }

  if (any(is.na(Esd))) warning("NAs in Var(alpha|k1)")

  # extract variances
  W = sdata[,.N,list(j1)][,N]
  model = list(nf=nf,A1=rliml$A1,B1=B1,A2=rliml$A2,B2=B2,EEsd=EEsd,Em=Em,eps1_sd=eps1_sd,
               eps2_sd=eps2_sd,Esd = Esd,Ns=W,Nm=Wm,EEm=EEm,withx=withx)

  # compute decomposition
  flog.info("computing wage variance decomposition")
  NNm = model$Nm; NNs = model$Ns
  NNm[!is.finite(NNm)]=0
  NNs[!is.finite(NNs)]=0
  stayer_share = sum(NNs)/(sum(NNs)+sum(NNm))

  if (bigk==0) {
    model$vdec   = m2.mini.vdec(model,1e6,stayer_share,"y1",do_interacted_reg=1)
  } else {
    model$vdec   = m2.mini.vdec(model,1e6,stayer_share,"y1",do_interacted_reg=0)
  }

  if (length(model0)>0) {
    rr = addmom(Em,model0$Em,"E(alpha|l1,l2)")
    rr = addmom(model$A1,model0$A1,"A1",rr)
    rr = addmom(model$A2,model0$A2,"A2",rr)
    rr = addmom(B1,model0$B1,"B1",rr)
    rr = addmom(B2,model0$B2,"B2",rr)
    rr = addmom(EEsd,model0$EEsd,"EEsd",rr)
    rr = addmom(eps1_sd,model0$eps1_sd,"eps_sd",rr)
    rr = addmom(Esd,model0$Esd,"Esd",rr)
    print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    catf("cor with true model:%f",cor(rr$val1,rr$val2))
  }

  return(model)
}

#' imputes data according to the model
#'
#' The results are stored in k_imp, y1_imp, y2_imp
#' this function simulates with the X or without depending on what model is passed.
#'
#' @export
m2.mini.impute.stayers <- function(model,sdatae) {
  A1      = model$A1
  B1      = model$B1
  A2      = model$A2
  B2      = model$B2
  Em      = model$Em
  Esd     = model$Esd
  eps1_sd  = model$eps1_sd
  eps2_sd  = model$eps2_sd


  # ------  impute K, Y1, Y4 on jdata ------- #
  sdatae.sim = copy(sdatae)
  if ("k_imp" %in% names(sdatae.sim)) sdatae.sim[,k_imp := NULL];
  if (model$withx) {
    sdatae.sim[, c('k_imp','y1_imp','y2_imp') := {
      ni = .N
      Ki  = Em[j1,x] + Esd[j1,x]*rnorm(ni)
      Y1  = A1[j1] + B1[j1]*Ki + eps1_sd[j1] * rnorm(ni)
      Y2  = A2[j1] + B2[j1]*Ki + eps2_sd[j1] * rnorm(ni)
      list(Ki,Y1,Y2)
    },list(j1,x)]
  } else {
    sdatae.sim[, c('k_imp','y1_imp','y2_imp') := {
      ni = .N
      Ki  = Em[j1] + Esd[j1]*rnorm(ni)
      Y1  = A1[j1] + B1[j1]*Ki + eps1_sd[j1] * rnorm(ni)
      Y2  = A2[j2] + B2[j2]*Ki + eps2_sd[j2] * rnorm(ni)
      list(Ki,Y1,Y2)
    },list(j1,j2)]
  }

  return(sdatae.sim)
}

#' imputes data for movers according to the model. The results are stored in k_imp, y1_imp, y2_imp
#' @export
m2.mini.impute.movers <- function(model,jdatae) {
  A1      = model$A1
  B1      = model$B1
  A2      = model$A2
  B2      = model$B2
  Em      = model$Em
  EEm     = model$EEm
  EEsd    = model$EEsd
  Esd     = model$Esd
  eps1_sd  = model$eps1_sd
  eps2_sd  = model$eps2_sd
  nf      = model$nf

  # ------  impute K, Y1, Y4 on jdata ------- #
  jdatae.sim = copy(jdatae)
  jdatae.sim[, c('k_imp','y1_imp','y2_imp') := {
    ni = .N
    jj = j1 + nf*(j2-1)
    Ki  = EEm[j1,j2] + rnorm(ni)*EEsd[j1,j2]
    # draw Y1, Y4
    Y1 = rnorm(ni)*eps1_sd[j1] + A1[j1] + B1[j1]*Ki
    Y2 = rnorm(ni)*eps2_sd[j2] + A2[j2] + B2[j2]*Ki
    list(Ki,Y1,Y2)
  },list(j1,j2)]

  return(jdatae.sim)
}


#' Computes the variance decomposition by simulation
#' @export
m2.mini.vdec <- function(model,nsim,stayer_share=1,ydep="y1",do_interacted_reg=1) {

  # simulate movers/stayers, and combine
  NNm = model$Nm
  NNs = model$Ns
  NNm[!is.finite(NNm)]=0
  NNs[!is.finite(NNs)]=0

  NNs = round(NNs*nsim*stayer_share/sum(NNs))
  NNm = round(NNm*nsim*(1-stayer_share)/sum(NNm))
  flog.info("computing var decomposition with ns=%i nm=%i",sum(NNs),sum(NNm))

  # we simulate from the model both movers and stayers
  sdata.sim = m2.mini.simulate.stayers(model,NNs)
  jdata.sim = m2.mini.simulate.movers(model,NNm)
  sdata.sim = rbind(sdata.sim[,list(j1,k=alpha,y1,y2)],jdata.sim[,list(j1,k=alpha,y1,y2)])
  proj_unc  = lin.proja(sdata.sim,ydep,"k","j1",do_interacted_reg=do_interacted_reg);

  return(proj_unc)
}



model.show <- function(model) {

  # compute the cavariance terms

  v1 = wtd.var(model$A1, model$Ns)
  v2 = wtd.mean(model$eps1_sd^2, model$Ns)
  v3 = wtd.var(model$B1 * model$Em, model$Ns ) + wtd.mean(model$B1^2 * model$Esd^2, model$Ns  )
  v4 = wtd.var(rbind(model$A1,model$B1*model$Em),model$Ns )
  catf(" V(alpha)=%f \t V(psi)=%f \t 2*cov")


}



m2.mini.test <- function() {

  model = m2.mini.new(10,serial = F)
  NNm   = array(200000/(model$nf^2),c(model$nf,model$nf))
  NNs  = array(300000/model$nf,model$nf)
  jdata = m2.mini.simulate.movers(model,NNm)
  sdata = m2.mini.simulate.stayers(model,NNs)




  # check moment restrictionB
  ggplot(jdata[,list(var(alpha)*B1[j1]^2+var(e1),var(y1),.N),list(j1,j2)],aes(x=V1,y=V2)) + geom_point() + geom_abline() + theme_bw()
  ggplot(jdata[,list(var(alpha)*B2[j2]^2+var(e2),var(y2),.N),list(j1,j2)],aes(x=V1,y=V2)) + geom_point() + geom_abline() + theme_bw()
  ggplot(jdata[,list(var(alpha)*B2[j2]*B1[j1],cov(y1,y2),.N),list(j1,j2)],aes(x=V1,y=V2)) + geom_point() + geom_abline() + theme_bw()

  jdata[(j1==6) & (j2==4),list(var(alpha)*B1[j1]+var(e1),var(y1),.N),list(j1,j2)]
  jdata[(j1==6) & (j2==4),sd(alpha)]

  # using estimated model
  load("../figures/src/minimodel.dat",verbose=T)



}

minimodel.exo.vardec <- function(res, alpha_l=res$alpha_l, NN=res$NN , nk = 10) {

  # we just simulate a cross-section and compute a linear variance decomposition

  # draw an l from p_l
  # draw an alpha from alpha | l
  # we draw from the mean: a(l) + b(l) alpha

  L  = length(res$a_l)
  Y  = array(0,c(sum(NN)))
  J  = array(0,c(sum(NN)))
  AL = array(0,c(sum(NN)))

  i =1
  for (l1 in 1:length(NN))  {
    I = i:(i+NN[l1] -1)
    alpha = res$alpha_l[l1] + sqrt(abs(res$alpha_var_l[l1]))*rnorm(NN[l1])
    Y[I]  = (res$a_l[l1] + res$b_l[l1]*alpha)
    J[I] = l1
    AL[I] = alpha
    i = i + NN[l1]
  }

  sim.data = data.table(j=J,y=Y,alpha=AL)
  sim.data[, k := ceiling( nk * ( rank(alpha)/.N))]

  # rr = sim.data[,mean(y),list(k,j)]
  # ggplot(rr,aes(x=j,y=V1,color=factor(k))) + geom_line() +theme_bw()

  fit = lm( y ~ alpha + factor(j),sim.data)
  sim.data$res = residuals(fit)

  pred = predict(fit,type = "terms")
  sim.data$k_hat = pred[,1]
  sim.data$l_hat = pred[,2]

  rr = list()

  rr$cc   = sim.data[,cov.wt( data.frame(y,k_hat,l_hat,res))$cov]
  rr$ccw  = rr$cc
  rr$ccm  = rr$cc

  fit2 = lm(y ~  factor(k) * factor(j), sim.data )
  rr$rsq1 = summary(fit)$r.squared
  rr$rsq2 = summary(fit2)$r.squared

  rr$pi_k_l = acast(sim.data[, .N , list(k,j)][ , V1:= N/sum(N),j],k~j,value.var="V1")
  rr$wage   = acast(sim.data[, mean(y) , list(k,j)],k~j,value.var="V1")

  return(rr)
}

#' extract Var(alpha|k1,k2) var
get.variances.movers <- function(jdata,sdata,model) {

  # USE MOVERS TO GET Var(epsilon|l1)
  setkey(jdata,j1,j2)
  YY1 = jdata[, var(y1)    ,list(j1,j2)][,V1]
  YY2 = jdata[, cov(y1,y2) ,list(j1,j2)][,V1]
  YY3 = jdata[, var(y2)    ,list(j1,j2)][,V1]
  W   = jdata[,.N,list(j1,j2)][,N]

  nf  = model$nf
  XX1 = array(0,c(nf^2,nf^2 + 2*nf))
  XX2 = array(0,c(nf^2,nf^2 + 2*nf))
  XX3 = array(0,c(nf^2,nf^2 + 2*nf))
  L = nf

  B1 = model$B1
  B2 = model$B2

  for (l1 in 1:nf) for (l2 in 1:nf) {
    ll = l2 + nf*(l1 -1)
    Ls = 0;

    # the Var(alpha|l1,l2)
    XX1[ll,ll] = B1[l1]^2
    XX2[ll,ll] = B1[l1]*B2[l2]
    XX3[ll,ll] = B2[l2]^2
    Ls= Ls+L^2

    # the var(nu|l)
    XX1[ll,Ls+l1] = 1;Ls= Ls+L
    XX3[ll,Ls+l2] = 1
  }

  XX = rbind(XX1,XX2,XX3)
  YY =     c(YY1,YY2,YY3)
  fit1 = lm.wfitnn(XX,YY,c(W,W,W))
  res = fit1$solution
  obj1 = t( YY - XX %*% res ) %*% diag(c(W,W,W)) %*% ( YY - XX %*% res )

  EEsd    = t(sqrt(pmax(array((res)[1:nf^2],c(nf,nf)),0)))
  nu1_sd = sqrt(pmax((res)[(nf^2+1):(nf^2+nf)],0))
  nu2_sd = sqrt(pmax((res)[(nf^2+nf+1):(nf^2+2*nf)],0))

  return(list(EEsd=EEsd,nu1_sd=nu1_sd,nu2_sd=nu2_sd,fit1=obj1))
}

#' plotting mini model
#' @export
m2.mini.plotw <- function(model,qt=6,getvals=F) {

  rr = data.table(l=1:length(model$A1),Em = model$Em,Esd = model$Esd,
                  N = model$Ns,A1=model$A1,B1=model$B1,A2=model$A2,B2=model$B2)

  alpha_m  = rr[, wtd.mean(Em,N)]
  alpha_sd = sqrt(rr[, wtd.mean(Esd^2,N) + wtd.var(Em,N) ])

  qts = qnorm( (1:qt)/(qt+1))

  rr2 = rr[, list( y1= (qts*alpha_sd + alpha_m)* B1+ A1 ,y2= (qts*alpha_sd + alpha_m)* B2+ A2, k=1:qt),l ]
  rr2 = melt(rr2,id.vars = c('l','k'))
  if (getvals==TRUE) {
    return(rr2)
  }
  ggplot(rr2,aes(x=l,y=value,color=factor(k))) + geom_line() + geom_point() + theme_bw() + facet_wrap(~variable)

}




#
extract.variances <- function(jdata,sdata,B1,B2,rr=list(),na.rm=F) {

  # compute the variance of espilon conditional on l1, l2
  alpha_var_l1l2 = acast(jdata[, list(V1 = mcov(y1,y2) /( B1[j1] * B2[j2])) , list(j1,j2)],j1~j2)

  # compute the variance of espilon conditional on l1, l2
  espilon_var_l1
  l2 = jdata[, list( mvar(y1) - mcov(y1,y2) * B1[j1]/B2[j2] ,.N) , list(j1,j2)]
  if (na.rm==TRUE) espilon_var_l1l2[ is.na(V1), V1:=0 ];

  # compute the variance of espilon conditional on l1, l2
  epsilon_var_l  = espilon_var_l1l2[ , weighted.mean(V1,N) , j1][, V1]

  # compute the variance of alpha conditional on l
  alpha_var_l = (sdata[,mvar(y),j][,V1] - epsilon_var_l)/b_l^2

  rr$espilon_var_l1l2 = acast(espilon_var_l1l2,j1~j2,value.var = "V1")
  rr$epsilon_var_l    = epsilon_var_l
  rr$alpha_var_l      = alpha_var_l
  rr$alpha_var_l1l2  = alpha_var_l1l2

  return(rr)
}

est.nls <- function() {

  # objective function given X
  eval_f <- function( x ) {

    B1  = x[1:L]
    B2  = x[(L+1):(2*L)]
    A1  = c(0,x[(2*L+1):(3*L-1)])
    A2  = x[(3*L):(4*L-1)]
    EEm = array( x[(4*L):(4*L + L^2-1)], c(L,L))

    # we sum square terms
    obj  = jdata[, (y1 - A1[j1] - B1[j1]*EEm[j1,j2])^2 + (y2 - A2[j2] - B2[j2]*EEm[j1,j2])^2,list(j1,j2)][,sum(V1^2)]

    gB1  = jdata[, -2*EEm[j1,j2]*(y1 - A1[j1] - B1[j1]*EEm[j1,j2]),list(j1,j2)][,sum(V1),j1]
    gB2  = jdata[, -2*EEm[j1,j2]*(y2 - A2[j2] - B2[j2]*EEm[j1,j2]),list(j1,j2)][,sum(V1),j2]
    gA1  = jdata[, -2*(y1 - A1[j1] - B1[j1]*EEm[j1,j2]),list(j1,j2)][,sum(V1),j1]
    gA2  = jdata[, -2*(y2 - A2[j2] - B2[j2]*EEm[j1,j2]),list(j1,j2)][,sum(V1),j2]
    gEEm = jdata[,  -2*B1[j1]*(y1 - A1[j1] - B1[j1]*EEm[j1,j2]) -2*B2[j2]*(y1 - A1[j1] - B1[j1]*EEm[j1,j2]),list(j1,j2)][,sum(V1),list(j1,j2)]

    return(list( "objective" = obj,
                 "gradient" = c(gB1,gB2,gA1[2:L],gA2)))
  }










}


#' estimates interacted model using quasi-linkelood on the data moments
#'
#' @export
m2.miniql.estimate <- function(jdata,sdata,norm=1,model0=c(),maxiter=1000,ncat=50,linear=F,init0=F) {

  nf = length(unique(sdata$j1))

  # we start by extracting the moments from the data
  Nm = acast(jdata[,.N,list(j1,j2)],j1~j2,fill=0,value.var = "N")
  M1 = acast(jdata[,mean(y1),list(j1,j2)],j1~j2,fill=0,value.var = "V1")
  M2 = acast(jdata[,mean(y2),list(j1,j2)],j1~j2,fill=0,value.var = "V1")
  V11 = acast(jdata[,var(y1),list(j1,j2)],j1~j2,fill=0,value.var = "V1")
  V22 = acast(jdata[,var(y2),list(j1,j2)],j1~j2,fill=0,value.var = "V1")
  V12 = acast(jdata[,cov(y1,y2),list(j1,j2)],j1~j2,fill=0,value.var = "V1")

  # then we initialize the parameters
  EEm  = array(rnorm(nf^2),c(nf,nf))
  EEsd = array(1,c(nf,nf))
  eps1_sd = rep(1,nf)
  eps2_sd = rep(1,nf)
  A  = rep(0,nf)
  B  = 0.5 + runif(nf)

  if (init0) {
    EEm  = model0$EEm
    EEsd = model0$EEsd
    #A = model0$A1
    B = model0$B1
    eps1_sd = model0$eps1_sd
    eps2_sd = model0$eps2_sd
  }

  SIG  = array(1,c(nf,nf))

  YY = rep(0,4)
  XX = array(0,c(4,2))
  XX[1,1]=1
  XX[4,2]=1
  SS = array(0,c(4,4))
  SS[1,1]=1

  DD1 = array(0,c(2*nf,2*nf))
  DD2 = rep(0,2*nf)
  lik0 = -Inf

  for (iter in 1:maxiter) {

    # E-step
    SIG    = ( spread(B^2/eps1_sd^2,2,nf)  + spread(B^2/eps2_sd^2,1,nf) + 1/EEsd^2 )^-1
    C0     = SIG * (   - spread(A*B/eps1_sd^2,2,nf)  - spread(A*B/eps2_sd^2,1,nf) + EEm/EEsd^2 )
    C1     = SIG * spread(B/eps1_sd^2,2,nf)
    C2     = SIG * spread(B/eps2_sd^2,1,nf)

    # evaluate the likelihood
    lik = 0
    for (k1 in 1:nf) for (k2 in 1:nf) {
      MM = c( M1[k1,k2] - A[k1] - B[k1]*EEm[k1,k2] , M2[k1,k2] - A[k2] - B[k2]*EEm[k1,k2] )
      OO = array( c(  B[k1]^2*EEsd[k1,k2]^2 + eps1_sd[k1]^2  ,
                      B[k1]*B[k2]*EEsd[k1,k2]^2,B[k1]*B[k2]*EEsd[k1,k2]^2,
                      B[k2]^2*EEsd[k1,k2]^2 + eps2_sd[k2]^2) ,c(2,2))
      VV = array(c(V11[k1,k2],V12[k1,k2],V12[k1,k2],V22[k1,k2]),c(2,2))
      lkk = - log(2*pi) - 1/2*log(det(OO)) -1/2* t(MM) %*% solve(OO) %*% MM -1/2* sum(diag(solve(OO) %*% VV))
      lik = lik + lkk*Nm[k1,k2]
    }
    dlik= lik0 - lik
    lik0 = lik
    if (iter%%ncat==0) flog.info("[%i] lik=%4.4f dlik=%4.4f",iter,lik,dlik);

    ntot = sum(Nm)

    # we construct the quad prob problem
    for (k1 in 1:nf) for (k2 in 1:nf) {
      # period 1 stuff
      YY[1]   = M1[k1,k2]
      YY[2]   = 1
      YY[3]   = 0

      XX[1,2] = C1[k1,k2]*M1[k1,k2] + C2[k1,k2]*M2[k1,k2] + C0[k1,k2]
      XX[2,2] = C1[k1,k2]
      XX[3,2] = -C2[k1,k2]

      SS[2,2] = V11[k1,k2]
      SS[2,3] = V12[k1,k2]
      SS[3,2] = V12[k1,k2]
      SS[3,3] = V22[k1,k2]
      SS[4,4] = SIG[k1,k2]

      I = 2*(k1-1)+1:2
      DD1[I,I] = DD1[I,I] + (t(XX) %*% SS %*% XX) * (Nm[k1,k2]/ntot)/eps1_sd[k1]^2
      DD2[I]   = DD2[I]   + (t(XX) %*% SS %*% YY) * (Nm[k1,k2]/ntot)/eps1_sd[k1]^2

      # period 2 stuff
      YY[1]   = M2[k1,k2]
      YY[2]   = 0
      YY[3]   = 1
      XX[2,2] = -C1[k1,k2]
      XX[3,2] = C2[k1,k2]

      I = 2*(k2-1)+1:2
      DD1[I,I] = DD1[I,I] + (t(XX) %*% SS %*% XX) * (Nm[k1,k2]/ntot)/eps2_sd[k2]^2
      DD2[I]   = DD2[I]   + (t(XX) %*% SS %*% YY) * (Nm[k1,k2]/ntot)/eps2_sd[k2]^2
    }

    # solve the problem
    Amat = array(0,c(2,2*nf)); Amat[1,1]=1; Amat[2,2]=1; bvec = c(0,1);
    if (linear) {
      Amat = array(0,c(nf+1,2*nf)); Amat[1,1]=1;
      for (k in 1:nf) Amat[k+1,2*k]=1;
      bvec=c(0,rep(1,nf))
    }

    fit = solve.QP(DD1,DD2,Amat = t(Amat),bvec =bvec,meq = nrow(Amat))
    R = array(fit$solution,c(2,nf))
    A = R[1,]
    B = R[2,]

    #Xtmp = C1*M1 + C2*M2 + C0
    #A = rowSums( (M1 - Xtmp)*Nm/spread(eps1_sd^2,2,nf) ) +  colSums( (M2 - Xtmp)*Nm/spread(eps2_sd^2,1,nf) )
    #A = A/( rowSums(Nm/spread(eps1_sd^2,2,nf)) +  colSums(Nm/spread(eps2_sd^2,1,nf)) )
    #A[1]=0

    # updating eps_sd
    #eps1_sd = as.numeric(sqrt(rowSums(Nm[k1,k2]*(  (1- spread(B,2,nf)*C1)^2*V11 -2*(1-spread(B,2,nf)*C1)*spread(B,2,nf)*C2*V12 + spread(B,2,nf)^2*C2^2*V22) )/rowSums(Nm)))
    #eps2_sd = as.numeric(sqrt(colSums(Nm[k1,k2]*(  (1- spread(B,1,nf)*C2)^2*V22 -2*(1-spread(B,1,nf)*C2)*spread(B,1,nf)*C1*V12 + spread(B,1,nf)^2*C1^2*V11) )/colSums(Nm)))

    # update mu and sigmas
    # EEm  = C0 + C1*M1 + C2*M2
    # EEsd = sqrt(SIG + C1^2*V11 + 2*C1*C2*V12 + C2^2*V22)
  } # we loop
  flog.info("[%i][final] lik=%4.4f dlik=%4.4f",iter,lik,dlik)

  model= list(A1=A,A2=A,B1=B,B2=B,EEm=EEm,EEsd=EEsd,eps1_sd=eps1_sd,eps2_sd=eps2_sd)


  if (length(model0)>0) {
    rr = addmom(model$EEm,model0$EEm,"E(alpha|l1,l2)")
    rr = addmom(model$A1,model0$A1,"A1",rr)
    rr = addmom(model$A2,model0$A2,"A2",rr)
    rr = addmom(model$B1,model0$B1,"B1",rr)
    rr = addmom(model$B2,model0$B2,"B2",rr)
    rr = addmom(model$EEsd,model0$EEsd,"EEsd",rr)
    rr = addmom(model$eps1_sd,model0$eps1_sd,"eps_sd",rr)
    #rr = addmom(Esd,model0$Esd,"Esd",rr)
    print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    catf("cor with true model:%f",cor(rr$val1,rr$val2))
  }

  return(model)
}

## RANDOM COEFFICIENT ESTIMATOR

#' Random coefficient estimator for the mini-model
#'
#' @param theta  value of the parameters
#' @param type   what part of the parameters to keep fixed (0: all parameters, 1: A,B, 2: A,eps, 3:eps, 4: homoskedastic)
#' @param model0 value for the parameters that are kept fixed
#' @param md  moments from the data (M1,M2,V11,V12,V22,Nm)
#' @param norm default is c(1,1) whicch gives the index for the normalization of A and B
#'
#' @export
m2.mini.rc.lik <- function(theta,type=0,model0=NA,md,norm=c(1,1)) {

  nf = dim(md$M1)[1]

  if (any(is.na(model0))) {
    model0 = list(A=rep(0,nf),B=rep(1,nf),eps_sd=rep(1,nf))
  }

  if (type==0) {
    A = rep(0,nf)
    A[setdiff(1:nf,norm[1])] = theta[1:(nf-1)]
    B = rep(1,nf)
    B[setdiff(1:nf,norm[2])] = theta[nf:(2*nf-2)]
    eps_sd = exp(theta[(2*nf-1):(3*nf-2)])
  }
  if (type==1) {
    A = rep(0,nf)
    A[setdiff(1:nf,norm[1])] = theta[1:(nf-1)]
    B = rep(1,nf)
    B[setdiff(1:nf,norm[2])] = theta[nf:(2*nf-2)]
    eps_sd = model0$eps_sd
  }
  if (type==2) {
    A = rep(0,nf)
    A[setdiff(1:nf,norm[1])] = theta[1:(nf-1)]
    B      = model0$B
    eps_sd = exp(theta[10:19])
  }
  if (type==3) {
    A      = model0$A
    B      = model0$B
    eps_sd = exp(theta[1:10])
  }
  if (type==4) {
    A = rep(0,nf)
    A[setdiff(1:nf,norm[1])] = theta[1:(nf-1)]
    B = rep(1,nf)
    B[setdiff(1:nf,norm[2])] = theta[nf:(2*nf-2)]
    eps_sd = rep(exp(theta[[19]]),10)
  }

  M1 = md$M1
  M2 = md$M2
  V11 = md$V11
  V22 = md$V22
  V12 = md$V12
  Nm = md$Nm

  lik_tot = 0
  SIG = diag(2)
  SIGir = diag(2)
  lik_all = array(0,c(nf,nf))

  for (k1 in 1:nf) for (k2 in 1:nf) {
    bv  = c(B[k1],B[k2])
    av  = c(A[k1],A[k2])
    SIGir[1,1] = eps_sd[k1]^-1 # inverse square root
    SIGir[2,2] = eps_sd[k2]^-1

    # data
    VV = array(c(V11[k1,k2],V12[k1,k2],V12[k1,k2],V22[k1,k2]),c(2,2))
    MM = c( M1[k1,k2], M2[k1,k2])

    # projector
    den = as.numeric( t(bv) %*% SIGir^2 %*% bv)
    WW = diag(2) - (SIGir %*% ( bv %*% t(bv) ) %*% SIGir) / den # @fixme: fix the scaling issue when the sigma are very large
    OO = SIGir %*% WW %*% SIGir

    # means
    lik_k1k2 = -0.5*log(2*pi) -0.5* t(MM - av) %*% OO %*% (MM - av)

    # variances
    lik_k1k2 = lik_k1k2 + 0.5*log(sum(diag(SIGir^2 %*% WW))) -1/2*sum(diag( VV %*% OO ))

    # add
    lik_all[k1,k2] = Nm[k1,k2]*as.numeric(lik_k1k2)
  }

  #flog.info("lik=%f",-lik)
  -sum(lik_all)
}

#' gradient function
#' @export
m2.mini.rc.lik.grad <- function(theta,type=0,model0=NA,md,norm=c(1,1)) {

  nf = dim(md$M1)[1]

  if (any(is.na(model0))) {
    model0 = list(A=rep(0,nf),B=rep(1,nf),eps_sd=rep(1,nf))
  }

  if (type==0) {
    A = rep(0,nf)
    A[setdiff(1:nf,norm[1])] = theta[1:(nf-1)]
    B = rep(1,nf)
    B[setdiff(1:nf,norm[2])] = theta[nf:(2*nf-2)]
    eps_sd = exp(theta[(2*nf-1):(3*nf-2)])
  }
  if (type==1) {
    A = rep(0,nf)
    A[setdiff(1:nf,norm[1])] = theta[1:(nf-1)]
    B = rep(1,nf)
    B[setdiff(1:nf,norm[2])] = theta[nf:(2*nf-2)]
    eps_sd = model0$eps_sd
  }
  if (type==2) {
    A = rep(0,nf)
    A[setdiff(1:nf,norm[1])] = theta[1:(nf-1)]
    B      = model0$B
    eps_sd = exp(theta[10:19])
  }
  if (type==3) {
    A      = model0$A
    B      = model0$B
    eps_sd = exp(theta[1:10])
  }
  if (type==4) {
    A = rep(0,nf)
    A[setdiff(1:nf,norm[1])] = theta[1:(nf-1)]
    B = rep(1,nf)
    B[setdiff(1:nf,norm[2])] = theta[nf:(2*nf-2)]
    eps_sd = rep(exp(theta[[19]]),10)
  }

  M1 = md$M1
  M2 = md$M2
  V11 = md$V11
  V22 = md$V22
  V12 = md$V12
  Nm = md$Nm

  lik_tot = 0
  SIG = diag(2)
  SIGir = diag(2)
  lik_all = array(0,c(nf,nf))

  diff1 = rep(0,3)
  diff2 = rep(0,3)
  grad  = array(0,c(3,nf))

  for (k1 in 1:nf) for (k2 in 1:nf) {
    a1 = A[k1]
    a2 = A[k2]
    b1 = B[k1]
    b2 = B[k2]
    s1 = eps_sd[k1]
    s2 = eps_sd[k2]

    # data
    v11 = V11[k1,k2]
    v12 = V12[k1,k2]
    v22 = V22[k1,k2]
    m1  = M1[k1,k2]
    m2  = M2[k1,k2]

    diff1[]=0
    diff2[]=0

    # term  -0.5* t(MM - av) %*% OO %*% (MM - av)
    diff1[1] =  ((a1/2 - m1/2)*(b1^2/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - 1))/s1^2 + ((a1 - m1)*(b1^2/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - 1))/(2*s1^2) + (b1*b2*(a2/2 - m2/2))/(s1^2*s2^2*(b1^2/s1^2 + b2^2/s2^2)) + (b1*b2*(a2 - m2))/(2*s1^2*s2^2*(b1^2/s1^2 + b2^2/s2^2))
    diff1[2] = (a1 - m1)*((((2*b1)/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - (2*b1^3)/(s1^4*(b1^2/s1^2 + b2^2/s2^2)^2))*(a1/2 - m1/2))/s1^2 + (b2*(a2/2 - m2/2))/(s1^2*s2^2*(b1^2/s1^2 + b2^2/s2^2)) - (2*b1^2*b2*(a2/2 - m2/2))/(s1^4*s2^2*(b1^2/s1^2 + b2^2/s2^2)^2)) - (a2 - m2)*((2*b1^2*b2*(a1/2 - m1/2))/(s1^4*s2^2*(b1^2/s1^2 + b2^2/s2^2)^2) - (b2*(a1/2 - m1/2))/(s1^2*s2^2*(b1^2/s1^2 + b2^2/s2^2)) + (2*b1*b2^2*(a2/2 - m2/2))/(s1^2*s2^4*(b1^2/s1^2 + b2^2/s2^2)^2))
    diff1[3] = (a2 - m2)*((2*b1^2*b2^2*(a2/2 - m2/2))/(s1^3*s2^4*(b1^2/s1^2 + b2^2/s2^2)^2) - (2*b1*b2*(a1/2 - m1/2))/(s1^3*s2^2*(b1^2/s1^2 + b2^2/s2^2)) + (2*b1^3*b2*(a1/2 - m1/2))/(s1^5*s2^2*(b1^2/s1^2 + b2^2/s2^2)^2)) - (a1 - m1)*((2*(a1/2 - m1/2)*(b1^2/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - 1))/s1^3 + (((2*b1^2)/(s1^3*(b1^2/s1^2 + b2^2/s2^2)) - (2*b1^4)/(s1^5*(b1^2/s1^2 + b2^2/s2^2)^2))*(a1/2 - m1/2))/s1^2 + (2*b1*b2*(a2/2 - m2/2))/(s1^3*s2^2*(b1^2/s1^2 + b2^2/s2^2)) - (2*b1^3*b2*(a2/2 - m2/2))/(s1^5*s2^2*(b1^2/s1^2 + b2^2/s2^2)^2))

    diff2[1] = ((a2/2 - m2/2)*(b2^2/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - 1))/s2^2 + ((a2 - m2)*(b2^2/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - 1))/(2*s2^2) + (b1*b2*(a1/2 - m1/2))/(s1^2*s2^2*(b1^2/s1^2 + b2^2/s2^2)) + (b1*b2*(a1 - m1))/(2*s1^2*s2^2*(b1^2/s1^2 + b2^2/s2^2))
    diff2[2] = (a2 - m2)*((((2*b2)/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - (2*b2^3)/(s2^4*(b1^2/s1^2 + b2^2/s2^2)^2))*(a2/2 - m2/2))/s2^2 + (b1*(a1/2 - m1/2))/(s1^2*s2^2*(b1^2/s1^2 + b2^2/s2^2)) - (2*b1*b2^2*(a1/2 - m1/2))/(s1^2*s2^4*(b1^2/s1^2 + b2^2/s2^2)^2)) - (a1 - m1)*((2*b1^2*b2*(a1/2 - m1/2))/(s1^4*s2^2*(b1^2/s1^2 + b2^2/s2^2)^2) - (b1*(a2/2 - m2/2))/(s1^2*s2^2*(b1^2/s1^2 + b2^2/s2^2)) + (2*b1*b2^2*(a2/2 - m2/2))/(s1^2*s2^4*(b1^2/s1^2 + b2^2/s2^2)^2))
    diff2[3] = (a1 - m1)*((2*b1^2*b2^2*(a1/2 - m1/2))/(s1^4*s2^3*(b1^2/s1^2 + b2^2/s2^2)^2) - (2*b1*b2*(a2/2 - m2/2))/(s1^2*s2^3*(b1^2/s1^2 + b2^2/s2^2)) + (2*b1*b2^3*(a2/2 - m2/2))/(s1^2*s2^5*(b1^2/s1^2 + b2^2/s2^2)^2)) - (a2 - m2)*((2*(a2/2 - m2/2)*(b2^2/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - 1))/s2^3 + (((2*b2^2)/(s2^3*(b1^2/s1^2 + b2^2/s2^2)) - (2*b2^4)/(s2^5*(b1^2/s1^2 + b2^2/s2^2)^2))*(a2/2 - m2/2))/s2^2 + (2*b1*b2*(a1/2 - m1/2))/(s1^2*s2^3*(b1^2/s1^2 + b2^2/s2^2)) - (2*b1*b2^3*(a1/2 - m1/2))/(s1^2*s2^5*(b1^2/s1^2 + b2^2/s2^2)^2))

    # term  0.5*log(sum(diag(SIGir^2 %*% WW))) -1/2*sum(diag( VV %*% OO ))
    diff1[2] = diff1[2] + (((2*b1)/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - (2*b1^3)/(s1^4*(b1^2/s1^2 + b2^2/s2^2)^2))/s1^2 - (2*b1*b2^2)/(s1^2*s2^4*(b1^2/s1^2 + b2^2/s2^2)^2))/(2*((b1^2/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - 1)/s1^2 + (b2^2/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - 1)/s2^2)) + (v11*((2*b1)/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - (2*b1^3)/(s1^4*(b1^2/s1^2 + b2^2/s2^2)^2)))/(2*s1^2) + (b2*v12)/(s1^2*s2^2*(b1^2/s1^2 + b2^2/s2^2)) - (2*b1^2*b2*v12)/(s1^4*s2^2*(b1^2/s1^2 + b2^2/s2^2)^2) - (b1*b2^2*v22)/(s1^2*s2^4*(b1^2/s1^2 + b2^2/s2^2)^2)
    diff1[3] = diff1[3] + (2*b1^3*b2*v12)/(s1^5*s2^2*(b1^2/s1^2 + b2^2/s2^2)^2) - (v11*((2*b1^2)/(s1^3*(b1^2/s1^2 + b2^2/s2^2)) - (2*b1^4)/(s1^5*(b1^2/s1^2 + b2^2/s2^2)^2)))/(2*s1^2) - (v11*(b1^2/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - 1))/s1^3 - ((2*(b1^2/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - 1))/s1^3 + ((2*b1^2)/(s1^3*(b1^2/s1^2 + b2^2/s2^2)) - (2*b1^4)/(s1^5*(b1^2/s1^2 + b2^2/s2^2)^2))/s1^2 - (2*b1^2*b2^2)/(s1^3*s2^4*(b1^2/s1^2 + b2^2/s2^2)^2))/(2*((b1^2/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - 1)/s1^2 + (b2^2/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - 1)/s2^2)) + (b1^2*b2^2*v22)/(s1^3*s2^4*(b1^2/s1^2 + b2^2/s2^2)^2) - (2*b1*b2*v12)/(s1^3*s2^2*(b1^2/s1^2 + b2^2/s2^2))
    diff2[2] = diff2[2] + (((2*b2)/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - (2*b2^3)/(s2^4*(b1^2/s1^2 + b2^2/s2^2)^2))/s2^2 - (2*b1^2*b2)/(s1^4*s2^2*(b1^2/s1^2 + b2^2/s2^2)^2))/(2*((b1^2/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - 1)/s1^2 + (b2^2/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - 1)/s2^2)) + (v22*((2*b2)/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - (2*b2^3)/(s2^4*(b1^2/s1^2 + b2^2/s2^2)^2)))/(2*s2^2) + (b1*v12)/(s1^2*s2^2*(b1^2/s1^2 + b2^2/s2^2)) - (b1^2*b2*v11)/(s1^4*s2^2*(b1^2/s1^2 + b2^2/s2^2)^2) - (2*b1*b2^2*v12)/(s1^2*s2^4*(b1^2/s1^2 + b2^2/s2^2)^2)
    diff2[3] = diff2[3] + (2*b1*b2^3*v12)/(s1^2*s2^5*(b1^2/s1^2 + b2^2/s2^2)^2) - (v22*((2*b2^2)/(s2^3*(b1^2/s1^2 + b2^2/s2^2)) - (2*b2^4)/(s2^5*(b1^2/s1^2 + b2^2/s2^2)^2)))/(2*s2^2) - (v22*(b2^2/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - 1))/s2^3 - ((2*(b2^2/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - 1))/s2^3 + ((2*b2^2)/(s2^3*(b1^2/s1^2 + b2^2/s2^2)) - (2*b2^4)/(s2^5*(b1^2/s1^2 + b2^2/s2^2)^2))/s2^2 - (2*b1^2*b2^2)/(s1^4*s2^3*(b1^2/s1^2 + b2^2/s2^2)^2))/(2*((b1^2/(s1^2*(b1^2/s1^2 + b2^2/s2^2)) - 1)/s1^2 + (b2^2/(s2^2*(b1^2/s1^2 + b2^2/s2^2)) - 1)/s2^2)) + (b1^2*b2^2*v11)/(s1^4*s2^3*(b1^2/s1^2 + b2^2/s2^2)^2) - (2*b1*b2*v12)/(s1^2*s2^3*(b1^2/s1^2 + b2^2/s2^2))

    # correct of the exponetional transform
    diff1[3] = diff1[3]*s1
    diff2[3] = diff2[3]*s2

    # add to main gradient
    grad[,k1] = grad[,k1] + Nm[k1,k2]*diff1
    grad[,k2] = grad[,k2] + Nm[k1,k2]*diff2
  }

  #flog.info("lik=%f",-lik)
  -c(grad[1,2:10],grad[2,2:10],grad[3,])
}


#' @export
m2.minirc.estimate <- function(cstats,mstats,method=1,verbose=FALSE) {

  setkey(cstats,j1)
  l1 = 0 #cstats[,wt.mean(m1,N)]
  l2 = 0 #cstats[,wt.mean(m2,N)]

  M1  = acast(mstats,j1~j2,value.var = "m1",fill=0)-l1
  M2  = acast(mstats,j1~j2,value.var = "m2",fill=0)-l2
  V11 = acast(mstats,j1~j2,value.var = "sd1",fill=0)^2
  V22 = acast(mstats,j1~j2,value.var = "sd2",fill=0)^2
  V12 = acast(mstats,j1~j2,value.var = "v12",fill=0)
  Nm  = acast(mstats,j1~j2,value.var = "N",fill=0)
  md  = list(M1=M1,M2=M2,V11=V11,V22=V22,V12=V12,Nm=Nm)
  nf = dim(M1)[[1]]

  if (method==1) {
    theta_1    =  c(rnorm(9),0.5+runif(9), rep(-1,10))
    res = optim(theta_1 ,fn = m2.mini.rc.lik,gr=m2.mini.rc.lik.grad,
                method ="BFGS",control=list(trace=100,maxit=1000,REPORT=1),md=md,type=0,norm=c(1,1))
  }

  # start with B=1
  if (method==2) {
    model0 = list(A=rep(0,nf),B=seq(1,1,l=nf),eps_sd=rep(1,nf))
    theta_0    = c(rnorm(nf-1),rep(-1,nf))
    res = optim(theta_0 ,fn = m2.mini.rc.lik,
                method ="BFGS",control=list(trace=100,maxit=300,REPORT=1),md=md,type=2,model0=model0)

    theta_1    = c(res$par[1:(nf-1)],model0$B[2:nf],res$par[nf:(2*nf-1)])
    res = optim(theta_1 ,fn = m2.mini.rc.lik,gr=m2.mini.rc.lik.grad,
                method ="BFGS",control=list(trace=100,maxit=300,REPORT=1),md=md,type=0)
  }

  A      = c(0,res$par[1:(nf-1)])
  B      = c(1,res$par[nf:(2*nf-2)])
  eps_sd = exp(res$par[(2*nf-1):(3*nf-2)])

  # extract EEm and EEsd
  T1 = (V11 - spread(eps_sd,1,nf)^2)/spread(B,1,nf)^2
  T2 = (V22 - spread(eps_sd,2,nf)^2)/spread(B,2,nf)^2
  EEsd = sqrt(1/( T1^-1 + T2^-1))
  R1 = (M1-spread(A,1,nf))/spread(B,1,nf)
  R2 = (M2-spread(A,2,nf))/spread(B,2,nf)
  EEm = EEsd^2 * ( T1^-1*R1 + T2^-1*R2)

  # finally do the same on stayers, recover Em and Esd
  Em  = (cstats$m1 -l1- A)/B
  Esd = sqrt((cstats$sd1^2 - eps_sd^2)/B^2)

  Esd[is.nan(Esd)]   = 0.001
  EEsd[is.nan(EEsd)] = 0.001
  EEm[is.nan(EEm)]   = 0.001

  model        = list(A1=A,A2=A,B1=B,B2=B,EEm=EEm,EEsd=EEsd,Em=Em,Esd=Esd,Nm=Nm,nf=10,Ns=cstats$N,eps1_sd=eps_sd,eps2_sd=eps_sd)
  stayer_share = sum(cstats$N)/(sum(cstats$N) + sum(mstats$N))
  model$vdec   = m2.mini.vdec(model,1e6,stayer_share,"y1")

  return(model)
}


