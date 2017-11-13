#' Model creation
#'
#' Creates a 2 period mini model. Can be set to be linear or not, with serial correlation or not
#' @export
m2.mini.new <- function(nf,linear=FALSE,serial=F,fixb=F) {
  m = list()
  m$nf = nf

  # generating intercepts and interaction terms
  m$A1 = rnorm(nf)
  m$B1 = exp(rnorm(nf)/2)
  m$A2 = rnorm(nf)
  m$B2 = exp(rnorm(nf)/2)

  # generate the E(alpha|l,l') and sd
  m$EEm  = spread(m$A1,2,nf) +  spread(m$A2,1,nf) + array(rnorm(nf^2),c(nf,nf))
  m$EEsd = exp(array(rnorm(nf^2),c(nf,nf))/2)

  # generate the variance of eps
  m$eps1_sd  = 0.1*exp(rnorm(nf)/2)
  m$eps2_sd  = m$eps1_sd
  m$eps_cor  = runif(1)

  # generate the E(alpha|l) and sd
  m$Em   = sort(rnorm(nf))
  m$Esd  = exp(rnorm(nf)/2)


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

  # order the clusters
  I = order(m$A1 + m$B1*m$Em)
  m$A1 = m$A1[I]
  m$B1 = m$B1[I]
  m$A2 = m$A2[I]
  m$B2 = m$B2[I]
  m$Em = m$Em[I]

  # set normalization
  bnorm = m$B1[1]
  m$B2  = m$B2/bnorm
  m$B1  = m$B1/bnorm
  m$Em  = m$Em*bnorm
  m$EEm = m$EEm*bnorm
  anorm = m$A1[1]
  m$A1  = m$A1-anorm
  m$A2  = m$A2-anorm

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

  L = max(J1)
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
m2.mini.estimate <- function(jdata,sdata,norm=1,model0=c(),method="ns",withx=FALSE) {

  # --------- use LIML on movers to get A1,B1,A2,B2 ---------- #
  Y1=jdata$y1;Y2=jdata$y2;J1=jdata$j1;J2=jdata$j2;
  nf = max(J1)
  if (method=="fixb") {
    rliml = m2.mini.liml.int.fixb(Y1,Y2,J1,J2,norm=norm)
  }  else if (method=="prof") {
    rliml = m2.mini.liml.int.prof(Y1,Y2,J1,J2,norm=norm)
  }  else if (method=="linear") {
    rliml = m2.mini.linear.int(Y1,Y2,J1,J2,norm=norm)
  } else {
    rliml = m2.mini.liml.int2(Y1,Y2,J1,J2,norm=norm)
  }
  B1 = rliml$B1
  B2 = rliml$B2
  N = length(Y1)

  # ---------  use stayers to get E[alpha|l] ------------- #
  EEm = jdata[, (mean(y1)- rliml$A1[j1])/rliml$B1[j1] + (mean(y2)- rliml$A2[j2])/rliml$B2[j2] ,list(j1,j2)]
  EEm = 0.5*acast(EEm,j1~j2,value.var="V1") # what do you think of this 0.5 right here?

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
  setkey(jdata,j1,j2)
  YY1 = c(acast(jdata[,mvar(y1),list(j1,j2)],j2~j1,fill = 0,value.var = "V1"))
  YY2 = c(acast(jdata[,mcov(y1,y2),list(j1,j2)],j2~j1,fill = 0,value.var = "V1")) #jdata[,mcov(y1,y2),list(j1,j2)][,V1]
  YY3 = c(acast(jdata[,mvar(y2),list(j1,j2)],j2~j1,fill = 0,value.var = "V1")) #jdata[,mvar(y2),list(j1,j2)][,V1]
  XX1 = array(0,c(nf^2, nf^2 + 2*nf))
  XX2 = array(0,c(nf^2, nf^2 + 2*nf))
  XX3 = array(0,c(nf^2, nf^2 + 2*nf))
  W = c(acast(jdata[,.N,list(j1,j2)],j2~j1,fill = 0,value.var = "N")) #jdata[,.N,list(j1,j2)][,N]

  for (l1 in 1:nf) for (l2 in 1:nf) {
    ll = l2 + nf*(l1 -1)
    if (W[ll]>0) {
      XX1[ll,ll                 ] = B1[l1]^2
      XX1[ll,nf^2 + l1          ] = 1
      XX2[ll,ll                 ] = B1[l1]*B2[l2]
      XX3[ll,ll                 ] = B2[l2]^2
      XX3[ll,nf^2 + nf + l2     ] = 1
    } else {
      XX1[ll,ll                 ] = 1
      XX3[ll,ll                 ] = 1
      W[ll]=0.1
    }
  }

  Wm = acast(jdata[,.N,list(j1,j2)],j1~j2,value.var="N")
  XX = rbind(XX1,XX2,XX3)

  res     = lm.wfitnn( XX, c(YY1,YY2,YY3), c(W,W,W))$solution
  EEsd    = t(sqrt(pmax(array((res)[1:nf^2],c(nf,nf)),0)))
  eps1_sd = sqrt(pmax((res)[(nf^2+1):(nf^2+nf)],0))
  eps2_sd = sqrt(pmax((res)[(nf^2+nf+1):(nf^2+2*nf)],0))
  if (any(is.na(eps1_sd*eps2_sd))) warning("NAs in Var(nu1|k1,k2)")
  if (any(is.na(EEsd))) warning("NAs in Var(alpha|k1,k2)")

  # ---------- STAYERS: use covariance restrictions  -------------- #
  # function which computes the variance terms
  getEsd <- function(sdata) {
    setkey(sdata,j1)
    mvar <- function(...) var(...)
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
  NNm = model$Nm; NNs = model$Ns
  NNm[!is.finite(NNm)]=0
  NNs[!is.finite(NNs)]=0
  stayer_share = sum(NNs)/(sum(NNs)+sum(NNm))
  model$vdec   = m2.mini.vdec(model,1e6,stayer_share,"y1")

  if (length(model0)>0) {
    rr = addmom(Em,model0$Em,"E(alpha|l1,l2)")
    rr = addmom(model$A1,model0$A1,"A1",rr)
    rr = addmom(model$A2,model0$A2,"A2",rr)
    rr = addmom(B1,model0$B1,"B1",rr)
    rr = addmom(B2,model0$B2,"B2",rr)
    rr = addmom(EEsd,model0$EEsd,"EEsd",rr)
    rr = addmom(eps1_sd,model0$eps1_sd,"eps1_sd",rr)
    rr = addmom(eps2_sd,model0$eps2_sd,"eps2_sd",rr)
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
m2.mini.vdec <- function(model,nsim,stayer_share=1,ydep="y1") {

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
  proj_unc  = lin.proja(sdata.sim,ydep,"k","j1");

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

#' Generate a linear projection decomposition for the model
#' with continuous worker hetergoneity
#'
#' @export
lin.proja <- function(sdata,y_col="y",k_col="k",j_col="j") {
  rr = list()

  sdata2 = copy(data.table(sdata))
  sdata2[,y_imp := get(y_col)]
  sdata2[,k_imp := get(k_col)]
  sdata2[,j     := get(j_col)]

  fit = lm(y_imp ~ k_imp + factor(j),sdata2)
  sdata2$res = residuals(fit)
  pred = predict(fit,type = "terms")
  sdata2$k_hat = pred[,1]
  sdata2$l_hat = pred[,2]

  rr$cc = sdata2[,cov.wt( data.frame(y_imp,k_hat,l_hat,res))$cov]
  rr$rsq1 = summary(fit)$r.squared

  fit2 = lm(y_imp ~ 0+  k_imp:factor(j) + factor(j),sdata2)
  rr$rsq2 = 1-mean(resid(fit2)^2)/var(sdata2$y_imp)

  get.stats <- function(cc) {
    r=list()
    den = cc[2,2] + cc[3,3] + 2 * cc[2,3]
    r$cor_kl = round(cc[2,3]/sqrt(cc[2,2]*cc[3,3]),4)
    r$cov_kl = 2*round(cc[2,3]/den,4)
    r$var_k  = round(cc[2,2]/den,4)
    r$var_l  = round(cc[3,3]/den,4)
    r$rsq    = round((cc[1,1] - cc[4,4])/cc[1,1],4)
    return(r)
  }

  rr$stats = get.stats(rr$cc)
  print(data.frame(rr$stats))
  rr$NNs = sdata[,.N,j1][order(j1)][,N]

  return(rr)
}
