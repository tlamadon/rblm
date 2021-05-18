# Creates, simulates and estimates
# the 2 period interactive model

require(reshape2)
require(limSolve)
require(Matrix)
require(data.table)
require(MASS)
require(quadprog)
catf <- function(...) cat(sprintf(...))

addmom <- function(A1,A2,name,rr=data.frame(),type="mean") {
  rr1 = melt(A1)
  rr2 = melt(A2)
  rr = rbind(rr,data.frame(name=name,val1=rr1$value,val2=rr2$val,type=type))
  return(rr)
}

#' Creates a 4 period mini model. Can be set to be linear or not, with serial correlation or not
#' @export
model.mini4.new <- function(nf,linear=FALSE,serial=F,r1=0.3 + 0.6*runif(1) ,r4=0.3 + 0.6*runif(1),fixb=F,eps_sd_factor=1) {
  m = list()
  m$nf = nf

  # generating intercepts and interaction terms
  m$A1  = rnorm(nf)
  m$B1  = exp(rnorm(nf)/2)
  m$A2  = rnorm(nf)
  m$A2s = rnorm(nf)
  m$B2  = exp(rnorm(nf)/2)
  m$A3  = rnorm(nf)
  m$A3s = rnorm(nf)
  m$B3  = exp(rnorm(nf)/2)
  m$A4  = rnorm(nf)
  m$B4  = exp(rnorm(nf)/2)

  if (fixb) {
    m$B2=m$B1
    m$B3=m$B1
    m$B4=m$B1
  }

  # generat the effect of the future firm
  m$C2 = rnorm(nf)
  m$C3 = rnorm(nf)

  # generat the auto-cor
  m$r1 = r1
  m$r4 = r4

  # generate the E(alpha|l,l') and sd
  m$EEm  = array(rnorm(nf^2),c(nf,nf))
  m$EEsd = exp(array(rnorm(nf^2),c(nf,nf))/2)

  # generate the variance of eps
  m$eps_sd  = exp(rnorm(nf)/2) * eps_sd_factor
  m$nu1_sd  = sqrt(  (1-r1^2) *  m$eps_sd^2 )
  m$nu2_sd  = sqrt(  (1-r4^2) *  m$eps_sd^2 )
  m$eps_cor = runif(1)
  m$eps2m_sd = array(exp(rnorm(nf*nf)/2) * eps_sd_factor,c(nf,nf))
  m$eps2s_sd = exp(rnorm(nf)/2) * eps_sd_factor

  # generate the E(alpha|l) and sd
  m$Em   = sort(rnorm(nf))
  m$Esd  = exp(rnorm(nf)/2)

  # set normalization
  rr = m$B1[1] - m$r1 * m$B2[1]
  m$B1 = m$B1/rr
  m$B2 = m$B2/rr
  m$B3 = m$B3/rr
  m$B4 = m$B4/rr
  m$A4[nf] = m$r4 * m$A3[nf]
  m$C2[1]  = 0
  m$C3[1]  = 0

  if (linear) {
    m$B1[]  = 1
    m$B2[]  = 1
    m$B3[]  = 1
    m$B4[]  = 1
    m$A1[1] = 0
  }

  return(m)
}

lm.wfitc <- function(XX,YY,rw,C1,C0,meq) {

  S = apply(abs(XX),2,max)
  XX2 = XX*spread(1/S,1,dim(XX)[1])
  C12 = C1*spread(1/S,1,dim(C1)[1])

  XXw      = diag(rw) %*% XX2
  Dq       = t(XXw) %*% XX2
  dq       = t(YY %*% XXw)

  # do quadprod
  fit      = solve.QP(Dq,dq,t(C12),C0,meq)

  # rescale
  fit$solution = fit$solution/S

  return(fit)
}

lm.wfitnn <- function(XX,YY,rw) {

  n = dim(XX)[2]
  XX2 = XX
#   S = apply(abs(XX),2,max)
#   XX2 = XX*spread(1/S,1,dim(XX)[1])
#   C12 = C1*spread(1/S,1,dim(C1)[1])

  XXw      = diag(rw) %*% XX2
  Dq       = t(XXw) %*% XX2
  dq       = t(YY %*% XXw)
  C1       = diag(n)
  C0       = rep(0,n)

  #fit      = qprog(Dq,dq,C1,C0)
  #fit$solution = as.numeric(fit$thetahat)

  fit      = solve.QP(Dq,dq,C1,C0,0)
  fit$solution = as.numeric(fit$solution)

  return(fit)
}

#' simulate movers according to the model
#' @export
model.mini4.simulate.movers <- function(model,NNm) {

  J1 = array(0,sum(NNm))
  J2 = array(0,sum(NNm))
  Y1 = array(0,sum(NNm))
  Y2 = array(0,sum(NNm))
  Y3 = array(0,sum(NNm))
  Y4 = array(0,sum(NNm))
  e2 = array(0,sum(NNm))
  e3 = array(0,sum(NNm))
  K  = array(0,sum(NNm))

  A1  = model$A1
  B1  = model$B1
  A2  = model$A2
  B2  = model$B2
  A3  = model$A3
  B3  = model$B3
  A4  = model$A4
  B4  = model$B4
  C2  = model$C2
  C3  = model$C3
  EEm = model$EEm
  EEsd= model$EEsd
  eps2m_sd=model$eps2m_sd
  eps3m_sd=model$eps3m_sd
  nu1_sd=model$nu1_sd
  nu2_sd=model$nu2_sd
  rho32m = model$rho32m
  nf  = model$nf
  r1 = model$r1
  r4 = model$r4

  i =1
  for (l1 in 1:nf) for (l2 in 1:nf) {
    I = i:(i+NNm[l1,l2]-1)
    ni = length(I)
    jj = l1 + nf*(l2 -1)
    J1[I] = l1
    J2[I] = l2

    # draw alpha
    K[I] = EEm[l1,l2] + EEsd[l1,l2]*rnorm(ni)

    # draw epsilon2 and epsilon3 corolated
    eps1 = eps2m_sd[l1,l2] * rnorm(ni)
    if (eps2m_sd[l1,l2]>0) {
      eps2 = eps3m_sd[l1,l2]/eps2m_sd[l1,l2]*rho32m*eps1    +    sqrt(  (1-rho32m^2) *  eps3m_sd[l1,l2]^2 ) * rnorm(ni)
    } else {
      eps2 = sqrt(  (1-rho32m^2) *  eps3m_sd[l1,l2]^2 ) * rnorm(ni)
    }
    # compute Y2 and Y3
    Y2[I]  = A2[l1] + B2[l1]*K[I] + C2[l2] + eps1
    Y3[I]  = A3[l2] + B3[l2]*K[I] + C3[l1] + eps2
    e2[I]  = eps1
    e3[I]  = eps2

    # compute Y1,Y4
    Y1[I]  = r1*Y2[I] + A1[l1] - r1*A2[l1] + (B1[l1] - r1 * B2[l1])*K[I] +  nu1_sd[l1] * rnorm(ni)
    Y4[I]  = r4*Y3[I] + A4[l2] - r4*A3[l2] + (B4[l2] - r4 * B3[l2])*K[I] +  nu2_sd[l2] * rnorm(ni)

    i = i + NNm[l1,l2]
  }

  jdatae = data.table(alpha=K,y1=Y1,y2=Y2,y3=Y3,y4=Y4,j1=J1,j2=J2,e2=e2,e3=e3)
  return(jdatae)
}

test.simulate.movers <- function() {
  jdata[, mean( model$A1[j1] + model$EEm[j1,j2] + model$r1*model$C2[j2] - y1   ),list(j1,j2)]
  jdata[, mean( model$A2[j1] + model$EEm[j1,j2] + model$C2[j2] - y2   ),list(j1,j2)]
  jdata[, mean( model$A3[j2] + model$EEm[j1,j2] + model$C3[j1] - y3   ),list(j1,j2)]
  jdata[, mean( model$A4[j2] + model$EEm[j1,j2] + model$r4*model$C3[j1] - y4   ),list(j1,j2)]



}

# simulate movers according to the model
model.mini4.impute.movers <- function(model,jdatae) {

  A1  = model$A1
  B1  = model$B1
  A2  = model$A2
  B2  = model$B2
  A3  = model$A3
  B3  = model$B3
  A4  = model$A4
  B4  = model$B4
  C2  = model$C2
  C3  = model$C3
  EEm = model$EEm
  EEsd= model$EEsd
  eps2m_sd=model$eps2m_sd
  eps3m_sd=model$eps3m_sd
  nu1_sd=model$nu1_sd
  nu2_sd=model$nu2_sd
  rho32m = model$rho32m
  nf  = model$nf
  r1 = model$r1
  r4 = model$r4

  # ------  impute K, Y1, Y4 on jdata -------
  jdatae.sim = copy(jdatae)
  jdatae.sim[, c('k_imp','y1_imp','y2_imp','y3_imp','y4_imp','e2','e3') := {
    ni = .N
    jj = j1 + nf*(j2-1)
    Ki = EEm[j1,j2] + EEsd[j1,j2]*rnorm(ni)

    # draw epsilon2 and epsilon3 corolated
    eps1 = eps2m_sd[j1,j2] * rnorm(ni)
    if (eps2m_sd[j1,j2]>0) {
      eps2 = eps3m_sd[j1,j2]/eps2m_sd[j1,j2]*rho32m*eps1    +    sqrt(  (1-rho32m^2) *  eps3m_sd[j1,j2]^2 ) * rnorm(ni)
    } else {
      eps2 = sqrt(  (1-rho32m^2) *  eps3m_sd[j1,j2]^2 ) * rnorm(ni)
    }
    Y2  = A2[j1] + B2[j1]*Ki + C2[j2] + eps1
    Y3  = A3[j2] + B3[j2]*Ki + C3[j1] + eps2

    # compute Y1,Y4
    Y1  = r1*Y2 + A1[j1] - r1*A2[j1] + (B1[j1] - r1 * B2[j1])*Ki +  nu1_sd[j1] * rnorm(ni)
    Y4  = r4*Y3 + A4[j2] - r4*A3[j2] + (B4[j2] - r4 * B3[j2])*Ki +  nu2_sd[j2] * rnorm(ni)

    list(Ki,Y1,Y2,Y3,Y4,eps1,eps2)
  },list(j1,j2)]

  return(jdatae.sim)
}

model.mini4.impute.stayers <- function(model,sdatae,ks=c(1,1)) {

  A1  = model$A1
  B1  = model$B1
  A2s  = model$A2s
  A2  = model$A2
  B2  = model$B2
  A3s  = model$A3s
  A3  = model$A3
  B3  = model$B3
  A4  = model$A4
  B4  = model$B4
  C2  = model$C2
  C3  = model$C3
  Em = model$Em
  Esd= model$Esd
  eps2s_sd=model$eps2s_sd
  eps3s_sd=model$eps3s_sd
  nu1_sd=model$nu1_sd
  nu2_sd=model$nu2_sd
  rho32s = model$rho32s
  nf  = model$nf
  r1 = model$r1
  r4 = model$r4

  # ------  impute K, Y1, Y4 on jdata -------
  sdatae.sim = copy(sdatae)
  sdatae.sim[, c('k_imp','y1_imp','y2_imp','y3_imp','y4_imp','e2','e3') := {
    ni = .N
    Ki = Em[j1] + Esd[j1]*rnorm(ni)

    # draw epsilon2 and epsilon3 corolated
    eps1 = eps2s_sd[j1] * rnorms(ni,ks)
    if (eps2s_sd[j1]>0) {
      eps2 = eps3s_sd[j1]/eps2s_sd[j1]*rho32s*eps1    +    sqrt(  (1-rho32s^2) *  eps3s_sd[j1]^2 ) * rnorms(ni,ks)
    } else {
      eps2 = sqrt(  (1-rho32s^2) *  eps3s_sd[j1]^2 ) * rnorms(ni,ks)
    }
    Y2  = A2s[j1] + B2[j1]*Ki + eps1
    Y3  = A3s[j1] + B3[j1]*Ki + eps2

    # compute Y1,Y4
    Y1  = r1*Y2 + A1[j1] - r1*A2[j1] + (B1[j1] - r1 * B2[j1])*Ki +  sqrt(  (1-r1^2) *  nu1_sd[j1]^2 ) * rnorms(ni,ks)
    Y4  = r4*Y3 + A4[j1] - r4*A3[j1] + (B4[j1] - r4 * B3[j1])*Ki +  sqrt(  (1-r4^2) *  nu2_sd[j1]^2 ) * rnorms(ni,ks)

    list(Ki,Y1,Y2,Y3,Y4,eps1,eps2)
  },list(j1)]

  return(sdatae.sim)
}


#' simulate stayers according to the model
#' @export
model.mini4.simulate.stayers <- function(model,NNs,ks=c(1,1)) {

  J1 = array(0,sum(NNs))
  Y1 = array(0,sum(NNs))
  Y2 = array(0,sum(NNs))
  Y3 = array(0,sum(NNs))
  Y4 = array(0,sum(NNs))
  e1 = array(0,sum(NNs))
  K  = array(0,sum(NNs))

  A1   = model$A1
  B1   = model$B1
  A2s  = model$A2s
  A2   = model$A2
  B2   = model$B2
  A3   = model$A3
  A3s  = model$A3s
  B3   = model$B3
  A4   = model$A4
  B4   = model$B4
  r1   = model$r1
  r4   = model$r4
  Em   = model$Em
  Esd  = model$Esd
  eps2s_sd=model$eps2s_sd
  eps3s_sd=model$eps3s_sd
  rho32s=model$rho32s
  nu1_sd=model$nu1_sd
  nu2_sd=model$nu2_sd

  nf  = model$nf

  i =1
  for (l1 in 1:nf) {
    I = i:(i+NNs[l1]-1)
    ni = length(I)
    J1[I] = l1

    # draw alpha
    K[I] = Em[l1] + Esd[l1]*rnorm(ni)

    # draw epsilon2 and epsilon3 corolated
    eps1 = eps2s_sd[l1] * rnorms(ni,ks)
    if (eps2s_sd[l1]>0) {
      eps2 = eps3s_sd[l1]/eps2s_sd[l1]*rho32s*eps1    +    sqrt(  (1-rho32s^2) *  eps3s_sd[l1]^2 ) * rnorms(ni,ks)
    } else {
      eps2 = sqrt(  (1-rho32s^2) *  eps3s_sd[l1]^2 ) * rnorms(ni,ks)
    }
    # compute Y2 and Y3
    Y2[I]  = A2s[l1] + B2[l1]*K[I] + eps1
    Y3[I]  = A3s[l1] + B3[l1]*K[I] + eps2

    # compute Y1,Y4
    Y1[I]  = r1*Y2[I] + A1[l1] - r1*A2[l1] + (B1[l1] - r1 * B2[l1])*K[I] + nu1_sd[l1] * rnorms(ni,ks)
    Y4[I]  = r4*Y3[I] + A4[l1] - r4*A3[l1] + (B4[l1] - r4 * B3[l1])*K[I] + nu2_sd[l1] * rnorms(ni,ks)

    i = i + NNs[l1]
  }

  sdatae = data.table(alpha=K,y1=Y1,y2=Y2,y3=Y3,y4=Y4,j1=J1)
  return(sdatae)
}

test.simulate.stayers <- function() {
  model = model.mini4.new(10,linear = T)
  NNs  = array(100000/model$nf,model$nf)
  sdata = model.mini4.simulate.stayers(model,NNs)

  # checking stayers
  sdata[, mean(  alpha - model$Em[j1] ),list(j1)]
  sdata[, mean( model$A2s[j1] + model$Em[j1]  - y2   ),list(j1)]
  sdata[, with(model,mean( r1*y2 + A1[j1]-r1*A2[j1] + (B1[j1]-r1*B2[j1])*Em[j1]  - y1   )),list(j1)]
  sdata[, with(model,mean( r4*y3 + A4[j1]-r4*A3[j1] + (B4[j1]-r4*B3[j1])*Em[j1]  - y4   )),list(j1)]
}



#' this function uses LIML to estimate a model with a given normalization
#' and we want to do it in a non-stationary way
model.mini2.liml.int <- function(Y1,Y2,J1,J2,norm=1,fixb=F,r1=0,r4=0,alpha=0) {

  L = max(J1)
  N = length(Y1)

  # ----- prepare the different regressors -----
  D1 = array(0,c(N,L))
  I = (1:N) + N*(J1-1)
  D1[I]=1
  D2 = array(0,c(N,L))
  I = (1:N) + N*(J2-1)
  D2[I]=1


  # contrust full matrix
  if (fixb) {
    X1 = D1*spread(Y1,2,L)/(1-r1)-D2*spread(Y2,2,L)/(1-r4)  # construct the matrix for the interaction terms
    xdim = L
  } else {
    X1 = cbind(D1*spread(Y1,2,L),D2*spread(Y2,2,L))  # construct the matrix for the interaction terms
    xdim = 2*L
  }
  X2 = cbind(D1,D2)                                # construct the matrix for the intercepts

  # checking that it fits
  if (FALSE) {
    eps = X1 %*% (1/model$B1) - X2 %*% c(  (model$A1 - model$r1*model$A2)/(model$B1*(1-model$r1)) ,  - (model$A4 - model$r4*model$A3)/(model$B1*(1-model$r4))   )
  }

  # set the normalizations
  Y  = -X1[,norm]
  X1 =  X1[,setdiff(1:xdim,norm)]
  X2 =  X2[,1:(2*L-1)]

  # construct the Z matrix of instruemts (j1,j2)
  Z = array(0,c(N,L^2))
  I = (1:N) + N*(J1-1 + L*(J2-1))
  Z[I]=1

  # make matrices sparse
  Z  = Matrix(Z,sparse=T)
  X2 = Matrix(X2,sparse=T)

  # ------ LIML procedure ----------
  R  = cbind(Y,X1)
  X1 = Matrix(X1,sparse=T)

  Wz = t(R) %*% R  -  ( t(R) %*% Z)  %*% solve( t(Z)  %*% Z  ) %*% ( t(Z) %*% R)
  Wx = t(R) %*% R  -  ( t(R) %*% X2) %*% solve( t(X2) %*% X2 ) %*% ( t(X2) %*% R)

  WW = Wx %*% solve(Wz)
  #WW = Wx %*% ginv(as.matrix(Wz))
  lambdas = eigen(WW)$values
  lambda  = min(lambdas) - alpha/( length(Y) - dim(X1)[2] - dim(X2)[2])

  XX = cbind(X1,X2)
  RR = (1-lambda)*t(XX) %*% XX + lambda * ( t(XX) %*% Z)  %*% solve( t(Z)  %*% Z  ) %*% ( t(Z) %*% XX)
  RY = (1-lambda)*t(XX) %*% Y  + lambda * ( t(XX) %*% Z)  %*% solve( t(Z)  %*% Z  ) %*% ( t(Z) %*% Y)

  # --------- extract the results ----------
  b_liml = as.numeric(solve(RR) %*% RY)
  tau = rep(0,L)
  for (i in 1:(L-1)) { tau[i+ (i>=norm)] = b_liml[i ]}
  tau[norm]=1
  B1 = 1/tau
  if (fixb==FALSE) {
    B2 = - 1/b_liml[L:(2*L-1)]
    A1 = - b_liml[(xdim):(xdim+L-1)] * B1
    A2 = c(b_liml[(xdim+L):(xdim + 2*L-2)],0) * B2
  } else {
    B2 = B1 * (1-r4)/(1-r1)
    A1 = - b_liml[(xdim):(xdim+L-1)]*B1*(1-r1)
    A2 = c(b_liml[(xdim+L):(xdim + 2*L-2)],0) * B2 *(1-r1)
  }

  return(list(A1=A1,A2=A2,B1=B1,B2=B2,b_liml=b_liml))
}

#' estimate 4 period model
#' @export
model.mini4.estimate <- function(jdata,sdata,r1,r4,model0=c(),fixb=F,norm=1,alpha=0) {

  # --------- use LIML on movers to get A1,B1,A2,B2 ----------
  Y1=jdata$y1;Y2=jdata$y2;Y3=jdata$y3;Y4=jdata$y4;J1=jdata$j1;J2=jdata$j2;
  nf = max(J1)
  rliml = model.mini2.liml.int(Y1-r1*Y2,Y4-r4*Y3,J1,J2,norm=norm,fixb = fixb,r1=r1,r4=r4,alpha=alpha)
  AA1 = rliml$A1
  AA2 = rliml$A2
  BB1 = rliml$B1
  BB2 = rliml$B2
  N   = length(Y1)
  Nm  = acast(jdata[,.N,list(j1,j2)],j1~j2,value.var ="N")
  Ns  = sdata[,.N,j1][order(j1)][,N]

  # get E[alpha|l1,l2]
  M1 = acast(jdata[,mean(y1-r1*y2),list(j1,j2)],j1~j2,value.var ="V1")
  M2 = acast(jdata[,mean(y4-r4*y3),list(j1,j2)],j1~j2,value.var ="V1")
  EEm = ( M1 - spread(AA1,2,nf)  )/(2*spread(BB1,2,nf)) + ( M2 - spread(AA2,1,nf)  )/(2*spread(BB2,1,nf))

  # get A,B,C using MEAN RESTRICTIONS
  setkey(jdata,j1,j2)
  YY1 = jdata[,mean(y1),list(j1,j2)][,V1]
  YY2 = jdata[,mean(y2),list(j1,j2)][,V1]
  YY3 = jdata[,mean(y3),list(j1,j2)][,V1]
  YY4 = jdata[,mean(y4),list(j1,j2)][,V1]
  W   = jdata[,.N,list(j1,j2)][,N]

  XX1 = array(0,c(nf^2,10*nf-2))
  XX2 = array(0,c(nf^2,10*nf-2))
  XX3 = array(0,c(nf^2,10*nf-2))
  XX4 = array(0,c(nf^2,10*nf-2))
  L = nf

  for (l1 in 1:nf) for (l2 in 1:nf) {
    ll = l2 + nf*(l1 -1)
    Ls = 0;

    # the A
    XX1[ll,Ls+l1] = 1;Ls= Ls+L
    XX2[ll,Ls+l1] = 1;Ls= Ls+L
    XX3[ll,Ls+l2] = 1;Ls= Ls+L
    XX4[ll,Ls+l2] = 1;Ls= Ls+L

    # the B
    XX1[ll,Ls+l1] = EEm[l1,l2];Ls= Ls+L
    XX2[ll,Ls+l1] = EEm[l1,l2];Ls= Ls+L
    XX3[ll,Ls+l2] = EEm[l1,l2];Ls= Ls+L
    XX4[ll,Ls+l2] = EEm[l1,l2];Ls= Ls+L

    # the C
    if (l2>1) {
      XX1[ll,Ls+l2-1] = r1;
      XX2[ll,Ls+l2-1] = 1 ;
    }
    Ls= Ls+L-1

    if (l1>1) {
      XX3[ll,Ls+l1-1] = 1;
      XX4[ll,Ls+l1-1] = r4 ;
    }
  }

  # construct the constraints
  CM1 = cbind(diag(nf),-r1*diag(nf), array(0,c(nf,8*nf-2)))
  CM2 = cbind(array(0,c(nf,2*nf)), -r4*diag(nf),diag(nf), array(0,c(nf,6*nf-2)) )
  CM3 = cbind(array(0,c(nf,4*nf)), diag(nf),-r1*diag(nf), array(0,c(nf,4*nf-2)))
  CM4 = cbind(array(0,c(nf,6*nf)), -r4*diag(nf),diag(nf), array(0,c(nf,2*nf-2)))

  XX = rbind(XX1,XX2,XX3,XX4)
  YY = c(YY1,YY2,YY3,YY4)
  CM = rbind(CM1,CM2,CM3,CM4)
  C0 = c(AA1,AA2,BB1,BB2)
  meq = 4*nf

  # add 2 contraints, B1=B2 and B3=B4
  if (fixb==T) {
    CM5 =  cbind(array(0,c(nf,4*nf)), diag(nf), -diag(nf), array(0,c(nf,4*nf-2)))
    CM6 =  cbind(array(0,c(nf,6*nf)), diag(nf), -diag(nf), array(0,c(nf,2*nf-2)))
    CM = rbind(CM,CM5,CM6)
    C0 = c(C0,rep(0,2*L))
    meq = meq+2*L
  }

  rs = lm.wfitc(XX,YY,c(W,W,W,W),CM,C0,meq)$solution

  Ls = 0;
  A1 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  A2 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  A3 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  A4 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  B1 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  B2 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  B3 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  B4 = rs[(Ls+1):(Ls+L)];Ls=Ls+L
  C2 = c(0,rs[(Ls+1):(Ls+L-1)]);Ls=Ls+L-1
  C3 = c(0,rs[(Ls+1):(Ls+L-1)]);Ls=Ls+L-1
  model = list(A1=A1,A2=A2,A3=A3,A4=A4,B1=B1,B2=B2,B3=B3,B4=B4,C2=C2,C3=C3,EEm=EEm,r1=r1,r4=r4,nf=nf)

  res_mean_stayer = model.mini4.getstayers(jdata,sdata,model)
  model$Em=res_mean_stayer$Em
  model$A2s=res_mean_stayer$A2s
  model$A3s=res_mean_stayer$A3s

  # get the variances
  res_var = model.mini4.getvar(jdata,sdata,model,r1,r4)
  model$Esd=res_var$Esd
  model$EEsd=res_var$EEsd
  model$nu1_sd=res_var$nu1_sd
  model$nu2_sd=res_var$nu2_sd
  model$fit1 = res_var$fit1
  model$fit2 = res_var$fit2
  model$Evar=res_var$Evar

  model$Ns=Ns
  model$Nm=Nm
  model$fit3 = model.mini4.getvar.stayers(jdata,sdata,model,r1,r4)
  model = model.mini.rho32.stayers(jdata,sdata,model)
  model = model.mini.rho32.movers(jdata,sdata,model)

  if (length(model0)>0) {
    rr = addmom(AA1,model0$A1 - model0$r1*model0$A2,"A1 - r1*A2")
    rr = addmom(AA2,model0$A4 - model0$r4*model0$A3,"A4 - r4*A3",rr)
    rr = addmom(BB1,model0$B1 - model0$r1*model0$B2,"B1 - r1*B2",rr)
    rr = addmom(BB2,model0$B4 - model0$r4*model0$B3,"B4 - r4*B3",rr)
    rr = addmom(EEm,model0$EEm,"E(alpha|l1,l2)",rr)
    rr = addmom(A1,model0$A1,"A1",rr)
    rr = addmom(A2,model0$A2,"A2",rr)
    rr = addmom(A3,model0$A3,"A3",rr)
    rr = addmom(A4,model0$A4,"A4",rr)
    rr = addmom(B1,model0$B1,"B1",rr)
    rr = addmom(B2,model0$B2,"B2",rr)
    rr = addmom(B3,model0$B3,"B3",rr)
    rr = addmom(B4,model0$B4,"B4",rr)
    rr = addmom(C2,model0$C2,"C2",rr)
    rr = addmom(C3,model0$C3,"C3",rr)
    rr = addmom(model$Em,model0$Em,"Em",rr)
    rr = addmom(model$Esd,model0$Esd,"Esd",rr)
    rr = addmom(model$EEsd,model0$EEsd,"EEsd",rr)
    rr = addmom(model$nu1_sd,model0$nu1_sd,"nu1_sd",rr)

    #rr = addmom(model$nu2_sd,model0$nu2_sd,"nu2_sd",rr)
    print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    catf("cor with true model:%f",cor(rr$val1,rr$val2))
  }

  return(model)
}

model.mini4lin.estimate <- function(jdata,sdata,r1,r4,model0=c(),norm=1) {

  # --------- use LIML on movers to get A1,B1,A2,B2 ----------
  Y1=jdata$y1;Y2=jdata$y2;Y3=jdata$y3;Y4=jdata$y4;J1=jdata$j1;J2=jdata$j2;
  nf = max(J1)
  N   = length(Y1)

  # get A,B,C using MEAN RESTRICTIONS
  setkey(jdata,j1,j2)
  YY1 = jdata[,mean(y1),list(j1,j2)][,V1]
  YY2 = jdata[,mean(y2),list(j1,j2)][,V1]
  YY3 = jdata[,mean(y3),list(j1,j2)][,V1]
  YY4 = jdata[,mean(y4),list(j1,j2)][,V1]
  W   = jdata[,.N,list(j1,j2)][,N]
  Nm  = acast(jdata[,.N,list(j1,j2)],j1~j2,value.var ="N")
  Ns  = sdata[,.N,j1][order(j1)][,N]

  XX1 = array(0,c(nf^2,6*nf-2 + nf^2))
  XX2 = array(0,c(nf^2,6*nf-2 + nf^2))
  XX3 = array(0,c(nf^2,6*nf-2 + nf^2))
  XX4 = array(0,c(nf^2,6*nf-2 + nf^2))
  L = nf

  for (l1 in 1:nf) for (l2 in 1:nf) {
    ll = l2 + nf*(l1 -1)
    Ls = 0;

    # the A
    XX1[ll,Ls+l1] = 1;Ls= Ls+L
    XX2[ll,Ls+l1] = 1;Ls= Ls+L
    XX3[ll,Ls+l2] = 1;Ls= Ls+L
    XX4[ll,Ls+l2] = 1;Ls= Ls+L

    # the EEm
    XX1[ll,Ls+ll] = 1;
    XX2[ll,Ls+ll] = 1;
    XX3[ll,Ls+ll] = 1;
    XX4[ll,Ls+ll] = 1;Ls= Ls+L^2

    # the C
    if (l2>1) {
      XX1[ll,Ls+l2-1] = r1;
      XX2[ll,Ls+l2-1] = 1 ;
    }
    Ls= Ls+L-1

    if (l1>1) {
      XX3[ll,Ls+l1-1] = 1;
      XX4[ll,Ls+l1-1] = r4 ;
    }
  }

  XX = rbind(XX1,XX2,XX3,XX4)
  XX = XX[,2:ncol(XX)] # removing the first intercept
  YY = c(YY1,YY2,YY3,YY4)

  rs = coef(lm.wfit(XX,YY,c(W,W,W,W)))

  Ls = 0;
  A1 = c(0,as.numeric(rs[(Ls+1):(Ls+L-1)]));Ls=Ls+L-1
  A2 = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L
  A3 = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L
  A4 = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L
  EEm= t(array( as.numeric(rs[(Ls+1):(Ls+L^2)]),c(L,L)));Ls=Ls+L^2
  C2 = c(0,as.numeric(rs[(Ls+1):(Ls+L-1)]));Ls=Ls+L-1
  C3 = c(0,as.numeric(rs[(Ls+1):(Ls+L-1)]));Ls=Ls+L-1

  model = list(A1=A1,A2=A2,A3=A3,A4=A4,B1=rep(1,nf),B2=rep(1,nf),B3=rep(1,nf),B4=rep(1,nf),C2=C2,C3=C3,EEm=EEm,r1=r1,r4=r4,nf=nf)
  res_mean_stayer = model.mini4.getstayers(jdata,sdata,model)
  model$Em=res_mean_stayer$Em
  model$A2s=res_mean_stayer$A2s
  model$A3s=res_mean_stayer$A3s

  # get the variances
  res_var = model.mini4.getvar(jdata,sdata,model,r1,r4)
  model$Esd=res_var$Esd
  model$EEsd=res_var$EEsd
  model$nu1_sd=res_var$nu1_sd
  model$nu2_sd=res_var$nu2_sd
  model$fit1 = res_var$fit1
  model$fit2 = res_var$fit2
  model$Evar=res_var$Evar

  # get the variances from the full stayers variance structure
  model$fit3 = model.mini4.getvar.stayers(jdata,sdata,model,r1,r4)

  model$Ns=Ns
  model$Nm=Nm

  model = model.mini.rho32.stayers(jdata,sdata,model)
  model = model.mini.rho32.movers(jdata,sdata,model)

  if (length(model0)>0) {
    rr = addmom(A1,model0$A1,"A1")
    rr = addmom(A2,model0$A2,"A2",rr)
    rr = addmom(A3,model0$A3,"A3",rr)
    rr = addmom(A4,model0$A4,"A4",rr)
    rr = addmom(C2,model0$C2,"C2",rr)
    rr = addmom(C3,model0$C3,"C3",rr)
    rr = addmom(EEm,model0$EEm,"EEm",rr)
    rr = addmom(model$Em,model0$Em,"Em",rr)
    rr = addmom(model$A2s,model0$A2s,"A2s",rr)
    rr = addmom(model$A3s,model0$A3s,"A3s",rr)
    rr = addmom(model$Esd,model0$Esd,"Esd",rr)
    rr = addmom(model$EEsd,model0$EEsd,"EEsd",rr)
    rr = addmom(model$nu1_sd,model0$nu1_sd,"nu1_sd",rr)
    rr = addmom(model$nu2_sd,model0$nu2_sd,"nu2_sd",rr)
    print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    catf("cor with true model:%f",cor(rr$val1,rr$val2))
  }

  return(model)
}

# extract Var(alpha|l1,l2) and Var(epsion|l)
model.mini4.getvar <- function(jdata,sdata,model,r1,r4) {

  # USE MOVERS TO GET Var(epsilon|l1)
  setkey(jdata,j1,j2)
  YY1 = jdata[, var(y1-r1*y2)           ,list(j1,j2)][,V1]
  YY2 = jdata[, cov(y1-r1*y2, y2)       ,list(j1,j2)][,V1]
  YY3 = jdata[, cov(y1-r1*y2, y3)       ,list(j1,j2)][,V1]
  YY4 = jdata[, cov(y1-r1*y2, y4-r4*y3) ,list(j1,j2)][,V1]
  YY5 = jdata[, cov(y2      , y4-r4*y3) ,list(j1,j2)][,V1]
  YY6 = jdata[, cov(y3      , y4-r4*y3) ,list(j1,j2)][,V1]
  YY7 = jdata[, var(y4-r4*y3)           ,list(j1,j2)][,V1]
  W   = jdata[,.N,list(j1,j2)][,N]

  nf = model$nf
  XX1 = array(0,c(nf^2,nf^2 + 2*nf))
  XX2 = array(0,c(nf^2,nf^2 + 2*nf))
  XX3 = array(0,c(nf^2,nf^2 + 2*nf))
  XX4 = array(0,c(nf^2,nf^2 + 2*nf))
  XX5 = array(0,c(nf^2,nf^2 + 2*nf))
  XX6 = array(0,c(nf^2,nf^2 + 2*nf))
  XX7 = array(0,c(nf^2,nf^2 + 2*nf))
  L = nf

  B1 = model$B1
  B2 = model$B2
  B3 = model$B3
  B4 = model$B4

  BB1 = B1 - model$r1*B2
  BB2 = B4 - model$r4*B3

  for (l1 in 1:nf) for (l2 in 1:nf) {
    ll = l2 + nf*(l1 -1)
    Ls = 0;

    # the Var(alpha|l1,l2)
    XX1[ll,ll] = BB1[l1]^2
    XX2[ll,ll] = BB1[l1]*B2[l1]
    XX3[ll,ll] = BB1[l1]*B3[l2]
    XX4[ll,ll] = BB1[l1]*BB2[l2]
    XX5[ll,ll] = BB2[l2]*B2[l1]
    XX6[ll,ll] = BB2[l2]*B3[l2]
    XX7[ll,ll] = BB2[l2]^2
    Ls= Ls+L^2

    # the var(nu|l)
    XX1[ll,Ls+l1] = 1;Ls= Ls+L
    XX7[ll,Ls+l2] = 1
  }

  XX = rbind(XX1,XX2,XX3,XX4,XX5,XX6,XX7)
  YY =     c(YY1,YY2,YY3,YY4,YY5,YY6,YY7)
  CM = diag(nf^2 + 2*nf)
  C0 = rep(0,nf^2 + 2*nf)

  fit1 = lm.wfitnn(XX,YY,c(W,W,W,W,W,W,W))
  #fit1 = lm.wfitc(XX,YY,c(W,W,W,W,W,W,W),CM,C0,0)
  res = fit1$solution
  obj1 = t( YY - XX %*% res ) %*% diag(c(W,W,W,W,W,W,W)) %*% ( YY - XX %*% res )

  EEsd    = t(sqrt(pmax(array((res)[1:nf^2],c(nf,nf)),0)))
  nu1_sd = sqrt(pmax((res)[(nf^2+1):(nf^2+nf)],0))
  nu2_sd = sqrt(pmax((res)[(nf^2+nf+1):(nf^2+2*nf)],0))

  # -------------  USE STAYERS TO GET Var(alpha|l1) ----------------
  setkey(sdata,j1)
  YY1 = sdata[, var(y1-r1*y2)           ,list(j1)][,V1] - nu1_sd^2
  YY2 = sdata[, cov(y1-r1*y2,y4-r4*y3)  ,list(j1)][,V1]
  YY3 = sdata[, var(y4-r4*y3)           ,list(j1)][,V1] - nu2_sd^2
  W   = sdata[,.N,list(j1)][,N]

  XX1 = array(0,c(nf,nf))
  XX2 = array(0,c(nf,nf))
  XX3 = array(0,c(nf,nf))
  L = nf

  for (l1 in 1:nf)  {
    # the Var(alpha|l1,l2)
    XX1[l1,l1] = BB1[l1]^2
    XX2[l1,l1] = BB1[l1]*BB2[l1]
    XX3[l1,l1] = BB2[l1]^2
  }

  XX = rbind(XX1,XX2,XX3)
  YY =     c(YY1,YY2,YY3)
  CM = diag(nf)
  C0 = rep(0,nf)
  #browser()

  fit2 = lm.wfitnn(XX,YY,c(W,W,W))
  res = fit2$solution

  obj2 = t( YY - XX %*% res ) %*% diag(c(W,W,W)) %*% ( YY - XX %*% res )
  Esd = sqrt(pmax(res,0))

  return(list(Esd=Esd,EEsd=EEsd,nu1_sd=nu1_sd,nu2_sd=nu2_sd,fit1=obj1,fit2=obj2,Evar=res))
}

test.model.mini4.getvar <- function() {
  model = model.mini4.new(10,r1=0.4,r4=0.7,eps_sd_factor = 1)
  NNm   = array(20000/(model$nf^2),c(model$nf,model$nf))
  NNs  = array(100000/model$nf,model$nf)
  jdata = model.mini4.simulate.movers(model,NNm)
  sdata = model.mini4.simulate.stayers(model,NNs)

  model1 = model.mini4.getvar(jdata,sdata,model,r1=model$r1,r4=model$r4)
  plot(model$nu1_sd,model1$nu1_sd)
  plot(model$nu2_sd,model1$nu2_sd)
}




model.mini.rho32.movers <- function(jdatae,sdata,model) {
  # get the variances of the epsilons
  eps2_sd = sqrt(pmax(acast(jdatae[, var(y2), list(j1,j2)],j1~j2,value.var = "V1") - model$EEsd^2,0))
  eps3_sd = sqrt(pmax(acast(jdatae[, var(y3), list(j1,j2)],j1~j2,value.var = "V1") - model$EEsd^2,0))

  # get the correlation
  rtmp = jdatae[,list(y = cov(y2,y3) - model$EEsd[j1,j2]^2, x = eps2_sd[j1,j2]*eps3_sd[j1,j2],.N),list(j1,j2)]
  fit = lm(y~0+x,rtmp,weights = rtmp$N)

  model$rho32m   = coef(fit)[[1]]
  model$eps2m_sd = eps2_sd
  model$eps3m_sd = eps3_sd

  return(model)
}

model.mini.rho32.stayers <- function(jdatae,sdata,model) {
  # get the variances of the epsilons
  eps2s_sd = sqrt(pmax(sdata[, var(y2), j1][,V1] - model$Esd^2,0))
  eps3s_sd = sqrt(pmax(sdata[, var(y3), j1][,V1] - model$Esd^2,0))

  # get the correlation
  rtmp = sdata[,list(y = cov(y2,y3) - model$Esd[j1]^2, x = eps2s_sd[j1]*eps3s_sd[j1],.N),j1]
  fit = lm(y~0+x,rtmp,weights = rtmp$N)

  model$rho32s   = coef(fit)[[1]]
  model$eps3s_sd = eps3s_sd
  model$eps2s_sd = eps2s_sd

  return(model)
}


# extract Var(alpha|l1,l2) and Var(epsion|l)
model.mini4.getvar.stayers <- function(jdata,sdata,model,r1,r4) {

  B1 = model$B1
  B2 = model$B2
  B3 = model$B3
  B4 = model$B4

  BB1 = B1 - r1*B2
  BB2 = B4 - r4*B3

  Esd = model$Esd
  nu1_sd = model$nu1_sd
  nu2_sd = model$nu2_sd

  # USE MOVERS TO GET Var(epsilon|l1)
  setkey(sdata,j1)
  R1  = sdata[, var(y1-r1*y2)            - BB1[j1]^2       * Esd[j1]^2 - nu1_sd[j1]^2 ,list(j1)][,V1]
  R2  = sdata[, cov(y1-r1*y2, y2)        - BB1[j1]*B2[j1]  * Esd[j1]^2                ,list(j1)][,V1]
  R3  = sdata[, cov(y1-r1*y2, y3)        - BB1[j1]*B3[j1]  * Esd[j1]^2                ,list(j1)][,V1]
  R4  = sdata[, cov(y1-r1*y2, y4-r4*y3)  - BB1[j1]*BB2[j1] * Esd[j1]^2                ,list(j1)][,V1]
  R5  = sdata[, cov(y2      , y4-r4*y3)  - BB2[j1]*B2[j1]  * Esd[j1]^2                ,list(j1)][,V1]
  R6  = sdata[, cov(y3      , y4-r4*y3)  - BB2[j1]*B3[j1]  * Esd[j1]^2                ,list(j1)][,V1]
  R7  = sdata[, var(y4-r4*y3)            - BB2[j1]^2       * Esd[j1]^2 - nu2_sd[j1]^2 ,list(j1)][,V1]
  W   = sdata[,.N,list(j1)][,N]

  obj1 = sum( c(R1,R2,R3,R4,R5,R6,R7)^2 * c(W,W,W,W,W,W,W) )

  return(obj1)
}



#' we extract E(alpha|l) and A2s and A3s jointly
model.mini4.getstayers <- function(jdata,sdata,model) {

  B1 = model$B1
  B2 = model$B2
  B3 = model$B3
  B4 = model$B4
  A1 = model$A1
  A2 = model$A2
  A3 = model$A3
  A4 = model$A4

  BB1 = B1 - model$r1*B2
  BB2 = B4 - model$r4*B3
  AA1 = A1 - model$r1*A2
  AA2 = A4 - model$r4*A3
  r1  = model$r1
  r4  = model$r4
  nf=model$nf

  setkey(sdata,j1)
  YY1 = sdata[, mean(y1 -r1*y2 -AA1[j1])     ,list(j1)][,V1]
  YY2 = sdata[, mean(y2)                     ,list(j1)][,V1]
  YY3 = sdata[, mean(y3)                     ,list(j1)][,V1]
  YY4 = sdata[, mean(y4-r4*y3 - AA2[j1])     ,list(j1)][,V1]
  W   = sdata[,.N,list(j1)][,N]

  XX1 = array(0,c(nf,3*nf))
  XX2 = array(0,c(nf,3*nf))
  XX3 = array(0,c(nf,3*nf))
  XX4 = array(0,c(nf,3*nf))
  L = nf

  for (l1 in 1:nf)  {
    XX1[l1,l1]   = BB1[l1]
    XX2[l1,l1]   = B2[l1]
    XX3[l1,l1]   = B3[l1]
    XX4[l1,l1]   = BB2[l1]

    # A2s and A3s
    XX2[l1,L+l1]   = 1
    XX3[l1,L+L+l1] = 1
  }

  XX = rbind(XX1,XX2,XX3,XX4)
  YY =     c(YY1,YY2,YY3,YY4)

  fit2 = lm.wfit(XX,YY,c(W,W,W,W))
  rs = coef(fit2)

  Ls = 0;
  Em  = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L
  A2s = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L
  A3s = as.numeric(rs[(Ls+1):(Ls+L)]);Ls=Ls+L

  return(list(Em=Em,A2s=A2s,A3s=A3s))
}


model.mini4.vardec <- function(model,verbose=F) {
  # compute the common b
  B0 = wtd.mean(model$B2*model$Esd^2,model$Ns)/wtd.mean(model$Esd^2,model$Ns)
  At = model$A2s + (model$B2 - B0)*model$Em

  # compute the variances
  v_psi    = wtd.var(At,model$Ns)
  v_alpha  = B0^2 * (  wtd.var(model$Em,model$Ns)  +  wtd.mean(model$Esd^2,model$Ns)   )
  v_2cov   = 2*cov.wt(cbind( At , B0 *model$Em),as.numeric(model$Ns))$cov[1,2]
  v_eps    = wtd.mean(model$eps2s_sd^2,model$Ns)
  v_cor    = v_2cov/(2 * sqrt(v_alpha*v_psi))

  if (verbose) {
    catf("v_psi= %4.4f   v_alpha= %4.4f   2cov= %4.4f   v_eps= %4.4f \n",v_psi,v_alpha,v_2cov,v_eps)
    catf("  psi= %4.4f   v_alpha= %4.4f   2cov= %4.4f   cor  = %4.4f \n",100*v_psi/(v_psi+v_alpha+v_2cov),100*v_alpha/(v_psi+v_alpha+v_2cov),100*v_2cov/(v_psi+v_alpha+v_2cov),v_cor)
    catf("V(E(alpha|l))= %f \n", wtd.var(model$Em,model$Ns))
  }

  return(list(v_psi=v_psi,v_alpha=v_alpha,v_2cov=v_2cov,v_eps=v_eps))
}

model.mini4.plot <- function(m) {

  dd = data.frame()
  L = m$nf
  dd = rbind(dd,data.frame(l=1:L,y=m$A1,name="a",t=1,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$A2,name="a",t=2,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$A3,name="a",t=3,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$A4,name="a",t=4,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$A2s,name="as",t=2,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$A3s,name="as",t=3,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$C2,name="xi",t=2,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$C3,name="xi",t=3,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$B1,name="B",t=1,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$B2,name="B",t=2,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$B3,name="B",t=3,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$B4,name="B",t=4,rho1=m$r1,rho2=m$r4),
             data.frame(l=1:L,y=m$Em - m$Em[1],name="alpha_l",t=2,rho1=m$r1,rho2=m$r4))

  ggplot(dd,aes(x=factor(l),y=y,group=t,color=factor(t))) + geom_line() + geom_point() +
    theme_bw() + facet_wrap(~name,scale="free")
}

model.mini4.plotw <- function(model) {

  rr = data.table(l=1:length(model$A1),Em = model$Em,Esd = model$Esd,
                  N = model$Ns,A1=model$A1,B1=model$B1,A2s=model$A2s,B2=model$B2)

  alpha_m  = rr[, wtd.mean(Em,N)]
  alpha_sd = sqrt(rr[, wtd.mean(Esd^2,N) + wtd.var(Em,N) ])

  rr2 = rr[, list( y1= (qnorm((1:10)/11)*alpha_sd + alpha_m)* B1+ A1 ,y2= (qnorm((1:10)/11)*alpha_sd + alpha_m)* B2+ A2s, k=1:10),l ]
  rr2 = melt(rr2,id.vars = c('l','k'))
  ggplot(rr2,aes(x=l,y=value,color=factor(k))) + geom_line() + geom_point() + theme_bw() + facet_wrap(~variable)

}


# extract Var(alpha|l1,l2) and Var(epsion|l)
model.mini4.getvar.stayers.unc <- function(sdata,r1,r4,weights=rep(1,7)) {

  # USE MOVERS TO GET Var(epsilon|l1)
  setkey(sdata,j1)
  YY1 = sdata[, var(y1-r1*y2)           ,j1][,V1]
  YY2 = sdata[, cov(y1-r1*y2, y2)       ,j1][,V1]
  YY3 = sdata[, cov(y1-r1*y2, y3)       ,j1][,V1]
  YY4 = sdata[, cov(y1-r1*y2, y4-r4*y3) ,j1][,V1]
  YY5 = sdata[, cov(y2      , y4-r4*y3) ,j1][,V1]
  YY6 = sdata[, cov(y3      , y4-r4*y3) ,j1][,V1]
  YY7 = sdata[, var(y4-r4*y3)           ,j1][,V1]
  W   = sdata[,.N,j1][,N]

  nf  = max(sdata$j1)
  XX1 = array(0,c(nf,3*nf))
  XX2 = array(0,c(nf,3*nf))
  XX3 = array(0,c(nf,3*nf))
  XX4 = array(0,c(nf,3*nf))
  XX5 = array(0,c(nf,3*nf))
  XX6 = array(0,c(nf,3*nf))
  XX7 = array(0,c(nf,3*nf))
  L = nf

  for (l1 in 1:nf) {
    ll = l1
    Ls = 0;

    # the Var(alpha|l1,l2)
    XX1[ll,ll] = (1-r1)^2
    XX2[ll,ll] = (1-r1)
    XX3[ll,ll] = (1-r1)
    XX4[ll,ll] = (1-r1)*(1-r4)
    XX5[ll,ll] = (1-r4)
    XX6[ll,ll] = (1-r4)
    XX7[ll,ll] = (1-r4)^2
    Ls= Ls+L

    # the var(nu|l)
    XX1[ll,Ls+l1] = 1;Ls= Ls+L
    XX7[ll,Ls+l1] = 1
  }

  XX = rbind(XX1,XX2,XX3,XX4,XX5,XX6,XX7)
  YY =     c(YY1,YY2,YY3,YY4,YY5,YY6,YY7)

  fit1 = lm.wfitnn(XX,YY,c(W*weights[1],W*weights[2],W*weights[3],W*weights[4],W*weights[5],W*weights[6],W*weights[7]))
  res  = fit1$solution
  obj1 = t( YY - XX %*% res ) %*% diag(c(W,W,W,W,W,W,W)) %*% ( YY - XX %*% res )

  EEsd    = t(sqrt(pmax(array((res)[1:nf^2],c(nf,nf)),0)))
  nu1_sd = sqrt(pmax((res)[(nf^2+1):(nf^2+nf)],0))
  nu2_sd = sqrt(pmax((res)[(nf^2+nf+1):(nf^2+2*nf)],0))

  return(list(EEsd=EEsd,nu1_sd=nu1_sd,nu2_sd=nu2_sd,fit1=obj1,res=colMeans(rdim(YY - XX %*% res,10,7)^2)))
}

# extract Var(alpha|l1,l2) and Var(epsion|l)
model.mini4.getvar.stayers.unc2 <- function(sdata,r1,r4,rt,weights=rep(1,7),model0=c()) {

  # USE MOVERS TO GET Var(epsilon|l1)
  setkey(sdata,j1)
  YY1 = sdata[, var(y1)                 ,j1][,V1]
  YY2 = sdata[, var(y2)                 ,j1][,V1]
  YY3 = sdata[, var(y3)                 ,j1][,V1]
  YY4 = sdata[, var(y4)                 ,j1][,V1]
  YY5 = sdata[, cov(y1,y2)              ,j1][,V1]
  YY6 = sdata[, cov(y1,y3)              ,j1][,V1]
  YY7 = sdata[, cov(y1,y4)              ,j1][,V1]
  YY8 = sdata[, cov(y2,y3)              ,j1][,V1]
  YY9 = sdata[, cov(y2,y4)              ,j1][,V1]
  YY0 = sdata[, cov(y3,y4)              ,j1][,V1]
  W   = sdata[,.N,j1][,N]

  nf  = max(sdata$j1)
  XX1 = array(0,c(nf,5*nf))
  XX2 = array(0,c(nf,5*nf))
  XX3 = array(0,c(nf,5*nf))
  XX4 = array(0,c(nf,5*nf))
  XX5 = array(0,c(nf,5*nf))
  XX6 = array(0,c(nf,5*nf))
  XX7 = array(0,c(nf,5*nf))
  XX8 = array(0,c(nf,5*nf))
  XX9 = array(0,c(nf,5*nf))
  XX0 = array(0,c(nf,5*nf))
  L = nf

  for (l1 in 1:nf) {
    ll = l1
    Ls = 0;

    # the Var(alpha|l1,l2)
    XX1[ll,ll] = 1
    XX2[ll,ll] = 1
    XX3[ll,ll] = 1
    XX4[ll,ll] = 1
    XX5[ll,ll] = 1
    XX6[ll,ll] = 1
    XX7[ll,ll] = 1
    XX8[ll,ll] = 1
    XX9[ll,ll] = 1
    XX0[ll,ll] = 1
    Ls= Ls+L

    # Var(nu_1|l)
    XX1[ll,Ls+l1] = 1;
    Ls= Ls+L

    # Var(nu_2|l)
    XX1[ll,Ls+l1] = r1^2;
    XX2[ll,Ls+l1] = 1;
    XX3[ll,Ls+l1] = rt^2
    XX4[ll,Ls+l1] = r1^2 * rt^2
    XX5[ll,Ls+l1] = r1
    XX6[ll,Ls+l1] = r1*rt
    XX7[ll,Ls+l1] = r1*rt*r4
    XX8[ll,Ls+l1] = rt
    XX9[ll,Ls+l1] = rt*r4
    XX0[ll,Ls+l1] = rt^2*r4
    Ls= Ls+L

    # Var(nu_3|l)
    XX3[ll,Ls+l1] = 1
    XX4[ll,Ls+l1] = r4^2
    XX0[ll,Ls+l1] = r4
    Ls= Ls+L

    # Var(nu_4|l)
    XX4[ll,Ls+l1] = 1
  }

  XX = rbind(XX1,XX2,XX3,XX4,XX5,XX6,XX7,XX8,XX9,XX0)
  YY =     c(YY1,YY2,YY3,YY4,YY5,YY6,YY7,YY8,YY9,YY0)

  fit1 = lm.wfitnn(XX,YY,c(W,W,W,W,W,W,W,W,W,W))
  res  = fit1$solution
  obj1 = t( YY - XX %*% res ) %*% diag(c(W,W,W,W,W,W,W,W,W,W)) %*% ( YY - XX %*% res )

  BEsd   = sqrt(pmax(res[(0*nf+1):(1*nf)],0))
  nu1_sd = sqrt(pmax(res[(1*nf+1):(2*nf)],0))
  nu2_sd = sqrt(pmax(res[(2*nf+1):(3*nf)],0))
  nu3_sd = sqrt(pmax(res[(3*nf+1):(4*nf)],0))
  nu4_sd = sqrt(pmax(res[(4*nf+1):(5*nf)],0))

  if(length(model0)>0) {
      rr = addmom(nu1_sd,model0$nu1_sd,"nu1")
      rr = addmom(nu4_sd,model0$nu2_sd,"nu4",rr)
      rr = addmom(nu2_sd,model0$eps2s_sd,"eps2",rr)
      rr = addmom(nu3_sd,sqrt(1-model0$rho32s^2)*model0$eps3s_sd,"eps3",rr)
      rr = addmom(BEsd,model0$Esd,"Esd",rr)

      print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))

  }


  return(list(BEsd=BEsd,nu1_sd=nu1_sd,nu2_sd=nu2_sd,nu4_sd=nu4_sd,nu3_sd=nu3_sd,fit1=obj1,res=colMeans(rdim(YY - XX %*% res,10,10)^2)))
}



# extract Var(alpha|l1,l2) and Var(epsion|l)
model.mini4.getvar.movers.unc <- function(jdata,r1,r4) {

  # USE MOVERS TO GET Var(epsilon|l1)
  setkey(jdata,j1,j2)
  YY1 = jdata[, var(y1-r1*y2)           ,list(j1,j2)][,V1]
  YY2 = jdata[, cov(y1-r1*y2, y2)       ,list(j1,j2)][,V1]
  YY3 = jdata[, cov(y1-r1*y2, y3)       ,list(j1,j2)][,V1]
  YY4 = jdata[, cov(y1-r1*y2, y4-r4*y3) ,list(j1,j2)][,V1]
  YY5 = jdata[, cov(y2      , y4-r4*y3) ,list(j1,j2)][,V1]
  YY6 = jdata[, cov(y3      , y4-r4*y3) ,list(j1,j2)][,V1]
  YY7 = jdata[, var(y4-r4*y3)           ,list(j1,j2)][,V1]
  W   = jdata[,.N,list(j1,j2)][,N]

  nf  = max(sdata$j1)
  XX1 = array(0,c(nf^2,nf^2+2*nf))
  XX2 = array(0,c(nf^2,nf^2+2*nf))
  XX3 = array(0,c(nf^2,nf^2+2*nf))
  XX4 = array(0,c(nf^2,nf^2+2*nf))
  XX5 = array(0,c(nf^2,nf^2+2*nf))
  XX6 = array(0,c(nf^2,nf^2+2*nf))
  XX7 = array(0,c(nf^2,nf^2+2*nf))
  L = nf

  for (l1 in 1:nf) for (l2 in 1:nf) {
    ll = l2 + nf*(l1 -1)
    Ls = 0;

    # the Var(alpha|l1,l2)
    XX1[ll,ll] = (1-r1)^2
    XX2[ll,ll] = (1-r1)
    XX3[ll,ll] = (1-r1)
    XX4[ll,ll] = (1-r1)*(1-r4)
    XX5[ll,ll] = (1-r4)
    XX6[ll,ll] = (1-r4)
    XX7[ll,ll] = (1-r4)^2
    Ls= Ls+L^2

    # the var(nu|l)
    XX1[ll,Ls+l1] = 1;Ls= Ls+L
    XX7[ll,Ls+l1] = 1
  }

  XX = rbind(XX1,XX2,XX3,XX4,XX5,XX6,XX7)
  YY =     c(YY1,YY2,YY3,YY4,YY5,YY6,YY7)

  fit1 = lm.wfitnn(XX,YY,c(W,W,W,W,W,W,W))
  res  = fit1$solution
  obj1 = t( YY - XX %*% res ) %*% diag(c(W,W,W,W,W,W,W)) %*% ( YY - XX %*% res )

  EEsd    = t(sqrt(pmax(array((res)[1:nf^2],c(nf,nf)),0)))
  nu1_sd = sqrt(pmax((res)[(nf^2+1):(nf^2+nf)],0))
  nu2_sd = sqrt(pmax((res)[(nf^2+nf+1):(nf^2+2*nf)],0))

  return(list(EEsd=EEsd,nu1_sd=nu1_sd,nu2_sd=nu2_sd,fit1=obj1))
}


model.mini4.test <- function() {

  # Estimate full model
  model  = model.mini4.new(10,r1=0.4,r4=0.8,eps_sd_factor = 1)
  NNm    = array(50000/(model$nf^2),c(model$nf,model$nf))
  NNs    = array(100000/model$nf,model$nf)
  jdata  = model.mini4.simulate.movers(model,NNm)
  sdata  = model.mini4.simulate.stayers(model,NNs)
  model1 = model.mini4.estimate(jdata,sdata,model$r1,model$r4,model0 = model);
  catf("r1=%f r4=%f r32=%f\n",model$r1,model$r4,model$eps_cor)
  tmp=model.mini4.vardec(model1,verbose=T);
  model$Ns = NNs; model$Nm = NNm
  tmp=model.mini4.vardec(model,verbose=T);

  model.mini4.plot(model1)

  # Estimate full model with fixed Bs
  model = model.mini4.new(10,fixb=T,r1=0.3,r4=0.6)
  NNm   = array(40000/(model$nf^2),c(model$nf,model$nf))
  NNs  = array(100000/model$nf,model$nf)
  jdata = model.mini4.simulate.movers(model,NNm)
  sdata = model.mini4.simulate.stayers(model,NNs)
  model1=model.mini4.estimate(jdata,sdata,model$r1,model$r4,model0 = model,fixb=T);
  catf("r1=%f r4=%f r32=%f\n",model$r1,model$r4,model$eps_cor)

  # Estimate linear model
  model = model.mini4.new(10,linear = T)
  NNm   = array(20000/(model$nf^2),c(model$nf,model$nf))
  NNs  = array(1000000/model$nf,model$nf)
  jdata = model.mini4.simulate.movers(model,NNm)
  sdata = model.mini4.simulate.stayers(model,NNs)
  model1 = model.mini4lin.estimate(jdata,sdata,model$r1,model$r4,model0=model)
  #model.mini4.vardec(model1)

  # do a grid on r1 and check fit
  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  for (i in 1:nrow(dd)) {
    model1 = model.mini4lin.estimate(jdata,sdata,dd$r1[i],dd$r4[i])
    dd$fit1[i] = model1$fit1
    dd$fit2[i] = model1$fit2
    dd$fit3[i] = model1$fit3
    vv = model.mini4.vardec(model1)
    dd$v_psi[i] = vv$v_psi
    dd$v_alpha[i] = vv$v_alpha
    dd$v_2cov[i] = vv$v_2cov
    dd$v_eps[i] = vv$v_eps
    catf("[%i/%i] r1=%f r4=%f fit3=%f \n",i,nrow(dd),dd$r1[i],dd$r4[i],model1$fit3)
  }
  ggplot(dd,aes(x=rho,y=fit3)) + geom_point() + geom_line() +
    geom_vline(xintercept=model$r1,color="red",linetype=2) + theme_bw()

  wireframe(fit3 ~ r1 + r4,dd)

  # testing

  model1 = model.mini4.getvar(jdata,sdata,model,r1=model$r1,r4=model$r4)
  plot(model$nu1_sd,model1$nu1_sd)


  model.mini4.getvar.stayers


  # load the data
  load("../figures/src/mini4p-linear-10.dat",verbose=T)
  model.mini4.vardec(model1)
  model.mini4.vardec(model2)


  Dmat = diag(c(1059.46526020664,2197.73648768006,5232.6060981711,5986.39143708911,4619.8257387377,2956.56357907282,5872.3510132438,5735.03707424825,2926.84811684315,1496.80395221145))
  Amat = diag(c(1.000000e+00,3.296800e-01,5.575260e-02,1.301464e-01,8.098885e-03,8.432904e-05,3.076671e-03,1.026615e-01,1.911045e-01,2.350272e+00))
  dvec = c(14.012253278123,-5.74691605253042,-63.4752811846924,-41.8857494513356,15.5578043209264,-20.73392047329,-70.9155514168175,-61.0727647947296,8.54506221361423,57.6462542228647)
  bvec = rep(0,10)
  solve.QP(Dmat,as.numeric(dvec),Amat,bvec,meq = 0,factorized = F)

  LowRankQP(Dmat,dvec,array(0,c(1,10)),c(0),uvec=rep(10000,10))

  require(LowRankQP)

  load("../quadprog_text.dat",verbose=T)
  solve.QP(Dmat,dvec,Amat,bvec)
  qprog(Dmat,dvec,Amat,bvec)

  # ========= REESTIMATING THE LINEAR MODEL ================
  load("../figures/src/mini4p-linear-10.dat",verbose=T)
  NNm   = model1_atm$Nm
  NNs   = model1_atm$Ns
  jdata = model.mini4.simulate.movers(model1_atm,NNm)
  sdata = model.mini4.simulate.stayers(model1_atm,NNs)
  model1 = model.mini4lin.estimate(jdata,sdata,model1_atm$r1,model1_atm$r4,model0=model1_atm)

  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  for (i in 1:nrow(dd)) {
    model1 = model.mini4lin.estimate(jdata,sdata,dd$r1[i],dd$r4[i])
    dd$fit1[i] = model1$fit1
    dd$fit2[i] = model1$fit2
    dd$fit3[i] = model1$fit3
  }

  dd$fit1b = pmin(dd$fit1,median(dd$fit1))
  wireframe(-fit1b~r1+r4,dd)
  dd$fit3b = pmin(dd$fit3,median(dd$fit3))
  wireframe(-fit3b~r1+r4,dd)

  g1=ggplot(dd,aes(x=r1,group=r4,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit1)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,group=r1,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="Red") + geom_vline(xintercept=dd$r4[which.min(dd$fit1)],linetype=2,color="blue")
  g3=ggplot(dd,aes(x=r1,group=r4,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit1)],linetype=2,color="blue")
  g4=ggplot(dd,aes(x=r4,group=r1,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="Red") + geom_vline(xintercept=dd$r4[which.min(dd$fit1)],linetype=2,color="blue")
  blm:::multiplot(g1,g2,g3,g4,cols=2)


  dd[which.min(dd$fit1),]

  # ========= REESTIMATING THE FIXB MODEL ================
  load("../figures/src/mini4p-fixb-10.dat",verbose=T)
  model0 = model1_atm
  NNm    = 20*model0$Nm
  NNs    = 20*model0$Ns
  jdata  = model.mini4.simulate.movers(model0,NNm)
  sdata  = model.mini4.simulate.stayers(model0,NNs)
  model1 = model.mini4.estimate(jdata,sdata,model0$r1,model0$r4,model0=model0,fixb=T,norm=1,alpha=0)

  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  pb <- txtProgressBar(min = 0, max = nrow(dd), style = 3)
  for (i in 1:nrow(dd)) {
    model1 = model.mini4.estimate(jdata,sdata,dd$r1[i],dd$r4[i],fixb=T)
    dd$fit1[i] = model1$fit1
    dd$fit2[i] = model1$fit2
    dd$fit3[i] = model1$fit3
    setTxtProgressBar(pb, i)
  }
  close(pb)

  dd$fit1b = pmin(dd$fit1,median(dd$fit1))
  wireframe(-fit1b~r1+r4,dd)
  dd$fit3b = pmin(dd$fit3,median(dd$fit3))
  wireframe(-fit3b~r1+r4,dd)

  g1=ggplot(dd,aes(x=r1,group=r4,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit1)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,group=r1,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="Red") + geom_vline(xintercept=dd$r4[which.min(dd$fit1)],linetype=2,color="blue")
  g3=ggplot(dd,aes(x=r1,group=r4,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit3)],linetype=2,color="blue")
  g4=ggplot(dd,aes(x=r4,group=r1,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="Red") + geom_vline(xintercept=dd$r4[which.min(dd$fit3)],linetype=2,color="blue")
  blm:::multiplot(g1,g2,g3,g4,cols=2)

  # ========= REESTIMATING THE FIXB MODEL ================
  load("../figures/src/mini4p-allb-10.dat",verbose=T)
  model0 = model1_atm
  NNm    = 20*model0$Nm
  NNs    = model0$Ns
  jdata  = model.mini4.simulate.movers(model0,NNm)
  sdata  = model.mini4.simulate.stayers(model0,NNs)
  model1 = model.mini4.estimate(jdata,sdata,model0$r1,model0$r4,model0=model0,fixb=F,alpha=1)

  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  pb <- txtProgressBar(min = 0, max = nrow(dd), style = 3)
  for (i in 1:nrow(dd)) {
    model1 = model.mini4.estimate(jdata,sdata,dd$r1[i],dd$r4[i],fixb=F,alpha=1)
    dd$fit1[i] = model1$fit1
    dd$fit2[i] = model1$fit2
    dd$fit3[i] = model1$fit3
    setTxtProgressBar(pb, i)
  }
  close(pb)

  dd$fit1b = pmin(dd$fit1,median(dd$fit1))
  wireframe(-fit1b~r1+r4,dd)
  dd$fit3b = pmin(dd$fit3,median(dd$fit3))
  wireframe(-fit3b~r1+r4,dd)

  g1=ggplot(dd,aes(x=r1,group=r4,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit1)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,group=r1,y=fit1)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="red") + geom_vline(xintercept=dd$r4[which.min(dd$fit1)],linetype=2,color="blue")
  g3=ggplot(dd,aes(x=r1,group=r4,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit3)],linetype=2,color="blue")
  g4=ggplot(dd,aes(x=r4,group=r1,y=fit3b)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="red") + geom_vline(xintercept=dd$r4[which.min(dd$fit3)],linetype=2,color="blue")
  blm:::multiplot(g1,g2,g3,g4,cols=2)

  dd[which.min(dd$fit),]

  # -------------- STAYERS ---------------------
  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  pb <- txtProgressBar(min = 0, max = nrow(dd), style = 3)
  dd[,paste("res",1:7,sep='')]=0
  for (i in 1:nrow(dd)) {
    model1 = model.mini4.getvar.stayers.unc(sdata,dd$r1[i],dd$r4[i])
    dd$fit[i] = model1$fit1
    dd[i,paste("res",1:7,sep='')]=model1$res
    setTxtProgressBar(pb, i)
  }
  close(pb)
  g1=ggplot(dd,aes(x=r1,group=r4,y=fit)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,group=r1,y=fit)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="red") + geom_vline(xintercept=dd$r4[which.min(dd$fit)],linetype=2,color="blue")
  blm:::multiplot(g1,g2,cols=2)
  dd[which.min(dd$fit),]

  # -------------- MOVERS ---------------------
  dd = expand.grid(r1=seq(0,0.99,l=25),r4=seq(0,0.99,l=25),fit=Inf)
  pb <- txtProgressBar(min = 0, max = nrow(dd), style = 3)
  for (i in 1:nrow(dd)) {
    model1 = model.mini4.getvar.movers.unc(jdata,dd$r1[i],dd$r4[i])
    dd$fit[i] = model1$fit1
    setTxtProgressBar(pb, i)
  }
  close(pb)
  g1=ggplot(dd,aes(x=r1,group=r4,y=fit)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,group=r1,y=fit)) + geom_line() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="red") + geom_vline(xintercept=dd$r4[which.min(dd$fit)],linetype=2,color="blue")
  blm:::multiplot(g1,g2,cols=2)




  # ------- SIMULATE FROM MIXTURE MODEL AND TRY TO EXTRACT rhos -----------
  load("../figures/src/em-endo-levelmodel-6x10.dat",verbose=T)
  model0 = res_het



}

# use ks=c(1,1) for normal and c(0.5,0.1) for high kurtosis
rnorms <- function(n,ks) {
  s2= ks[1]
  lambda = ks[2]
  s1 = sqrt(1- (1-lambda)*s2^2)/sqrt(lambda)
  s1 = sqrt(1- (1-lambda)*s2^2)/sqrt(lambda)
  co = runif(n)<lambda
  X  = rnorm(n)*( co*s1 + (1-co)*s2)
  return(X)
}



test.kurtosis <- function() {
  # simulate from a function with given kurosis and skewness and mean 0 and variance 1

  # 2 point normal mixture lamba, m1, m2, s1, s2
  # but mean 0 and variance 1 so we have constraints!
  # take (m1,lambda,s1,s2)

  sim.mix <- function(lambda,s2,n=1000) {
    s1 = sqrt(1- (1-lambda)*s2^2)/sqrt(lambda)
    co = runif(n)<lambda
    X  = rnorm(n)*( co*s1 + (1-co)*s2)
    return(c(var(X),skewness(X),kurtosis(X)))
  }

  rnorms <- function(n,lambda,s2) {
    s1 = sqrt(1- (1-lambda)*s2^2)/sqrt(lambda)
    s1 = sqrt(1- (1-lambda)*s2^2)/sqrt(lambda)
    co = runif(n)<lambda
    X  = rnorm(n)*( co*s1 + (1-co)*s2)
    return(X)
  }

  ff  <- function(theta) sum( (sim.mix(theta[1],theta[2],theta[3],theta[4]) - c(1,2,10))^2)
  ff2 <- function(theta) sim.mix(theta[1],theta[2],theta[3],theta[4])
  res = optim(c(0,0.1,1,10),ff)
  ff2(res$par)


}

test.getrhos.unc2 <- function() {
  load("../figures/src/mini4p-linear-10.dat",verbose=T)
  model0 = model1_atm
  NNm   = model0$Nm
  NNs   = 15*model0$Ns
  jdata = model.mini4.simulate.movers(model0,NNm)

  model0$r1=0.6
  model0$r4=0.7
  model0$rho32s=0.4
  sdata = model.mini4.simulate.stayers(model0,NNs)

  # run at true rhos
  res = model.mini4.getvar.stayers.unc2(sdata,model0$r1,model0$r4,model0$rho32s,model0=model0)

  nr = 15
  dd = expand.grid(r1=seq(0,0.99,l=nr),r4=seq(0,0.99,l=nr),rt=seq(0,0.99,l=nr),fit=Inf)
  pb <- txtProgressBar(min = 0, max = nrow(dd), style = 3)
  dd[,paste("res",1:10,sep='')]=0
  for (i in 1:nrow(dd)) {
    model1 = model.mini4.getvar.stayers.unc2(sdata,dd$r1[i],dd$r4[i],dd$rt[i])
    dd$fit[i] = model1$fit1
    dd[i,paste("res",1:10,sep='')]=model1$res
    setTxtProgressBar(pb, i)
  }
  close(pb)
  g1=ggplot(dd,aes(x=r1,y=fit)) + geom_point() + theme_bw() + geom_vline(xintercept=model0$r1,linetype=2,color="red") + geom_vline(xintercept=dd$r1[which.min(dd$fit)],linetype=2,color="blue")
  g2=ggplot(dd,aes(x=r4,y=fit)) + geom_point() + theme_bw() + geom_vline(xintercept=model0$r4,linetype=2,color="red") + geom_vline(xintercept=dd$r4[which.min(dd$fit)],linetype=2,color="blue")
  g3=ggplot(dd,aes(x=rt,y=fit)) + geom_point() + theme_bw() + geom_vline(xintercept=model0$rho32s,linetype=2,color="red") + geom_vline(xintercept=dd$rt[which.min(dd$fit)],linetype=2,color="blue")
  multiplot(g1,g2,g3,cols=3)
  dd[which.min(dd$fit),]

  # testing strong kurtosis effects
  load("../figures/src/mini4p-fixb-10.dat",verbose=T)
  model0 = model1_atm
  NNm   = model0$Nm
  NNs   = 5*model0$Ns
  jdata = model.mini4.simulate.movers(model0,NNm)

  model0$r1=0.4
  model0$r4=0.6
  model0$rho32s=0.3
  sdata = model.mini4.simulate.stayers(model0,NNs,ks = c(0.5,0.1))
  sdata = model.mini4.simulate.stayers(model0,NNs,ks = c(1,1))

  # estimate mixture of normal
  modelr = em.endo.level.new(6,10)
  ctrl   = em.control(nplot=10,check_lik=F,fixb=T,est_rho=T)
  res    = em.endo.level.rhoext.stayers.unc(jdata,sdata,modelr,ctrl)


}
