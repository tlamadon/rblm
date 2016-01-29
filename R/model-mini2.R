# Creates, simulates and estimates
# the 2 period interactive model

require(reshape2)
require(limSolve)
require(Matrix)
require(data.table)
require(MASS)

addmom <- function(A1,A2,name,rr=data.frame(),type="mean") {
  rr1 = melt(A1)
  rr2 = melt(A2)
  rr = rbind(rr,data.frame(name=name,val1=rr1$value,val2=rr2$val,type=type))
  return(rr)
}

#' Creates a 2 period mini model. Can be set to be linear or not, with serial correlation or not
model.mini2.new <- function(nf,linear=FALSE,serial=F) {
  m = list()
  m$nf = nf

  # generating intercepts and interaction terms 
  m$A1 = rnorm(nf)
  m$B1 = exp(rnorm(nf)/2) 
  m$A2 = rnorm(nf)
  m$B2 = exp(rnorm(nf)/2) 

  # generate the E(alpha|l,l') and sd
  m$EEm  = array(rnorm(nf^2),c(nf,nf))
  m$EEsd = exp(array(rnorm(nf^2),c(nf,nf))/2)

  # generate the variance of eps
  m$eps_sd  = exp(rnorm(nf)/2) 
  m$eps_cor = runif(1)

  # generate the E(alpha|l) and sd
  m$Em   = sort(rnorm(nf))
  m$Esd  = exp(rnorm(nf)/2)

  # set normalization
  m$B2 = m$B2/m$B1[1]
  m$B1 = m$B1/m$B1[1]
  m$A2[m$nf] = 0
  
  if (linear) {
    m$B1[] = 1
    m$B2[] = 1
  }

  if (serial==F) {
    m$eps_cor = 0
  }

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

# simulate movers according to the model
model.mini2.simulate.movers <- function(model,NNm) {

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
  eps_sd=model$eps_sd
  eps_cor=model$eps_cor
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
    eps1 = eps_sd[l1] * rnorm(ni)
    eps2 = eps_sd[l2]/eps_sd[l1]*eps_cor*eps1    +    sqrt(  (1-eps_cor^2) *  eps_sd[l2]^2 ) * rnorm(ni)

    Y1[I]  = A1[l1] + B1[l1]*K[I] + eps1 
    Y2[I]  = A2[l2] + B2[l2]*K[I] + eps2 
    e1[I]  = eps1
    e2[I]  = eps2
        
    i = i + NNm[l1,l2]
  }
  
  jdatae = data.table(alpha=K,y1=Y1,y2=Y2,j1=J1,j2=J2,e1=e1,e2=e2)
  return(jdatae)   
}

# simulate movers according to the model
model.mini2.simulate.stayers <- function(model,NNs) {

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
  eps_sd=model$eps_sd

  nf  = model$nf
  
  i =1 
  for (l1 in 1:nf) {
    I = i:(i+NNs[l1]-1)
    ni = length(I)
    J1[I] = l1
    
    # draw alpha
    K[I] = Em[l1] + Esd[l1]*rnorm(ni)
    
    # draw Y2, Y3
    e1 = rnorm(ni) * eps_sd[l1]
    e2 = rnorm(ni) * eps_sd[l1]
    Y1[I]  = A1[l1] + B1[l1]*K[I] + e1
    Y2[I]  = A2[l1] + B2[l1]*K[I] + e2
        
    i = i + NNs[l1]
  }
  
  sdatae = data.table(alpha=K,y1=Y1,y2=Y2,j1=J1)
  return(sdatae)  
}


#' this function uses LIML to estimate a model with a given normalization
#' and we want to do it in a non-stationary way
model.mini2.liml.int <- function(Y1,Y2,J1,J2,norm=1) {
    
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
  X1 = cbind(D1*spread(Y1,2,L),D2*spread(Y2,2,L))  # construct the matrix for the interaction terms
  X2 = cbind(D1,D2)                                # construct the matrix for the intercepts
  
  # set the normalizations
  Y  = -X1[,norm]
  X1 =  X1[,setdiff(1:(2*L),norm)]
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
  lambda  = min(lambdas)
  
  XX = cBind(X1,X2)
  RR = (1-lambda)*t(XX) %*% XX + lambda * ( t(XX) %*% Z)  %*% solve( t(Z)  %*% Z  ) %*% ( t(Z) %*% XX)
  RY = (1-lambda)*t(XX) %*% Y  + lambda * ( t(XX) %*% Z)  %*% solve( t(Z)  %*% Z  ) %*% ( t(Z) %*% Y)
  
  # --------- extract the results ----------
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


model.mini2.estimate <- function(jdata,sdata,norm=1,model0=c()) {

  # --------- use LIML on movers to get A1,B1,A2,B2 ----------
  Y1=jdata$y1;Y2=jdata$y2;J1=jdata$j1;J2=jdata$j2;
  nf = max(J1)
  rliml = model.mini2.liml.int(Y1,Y2,J1,J2,norm=norm)
  B1 = rliml$B1
  B2 = rliml$B2
  N = length(Y1)
  
  # ---------  use stayers to get E[alpha|l] -------------
  Em = (sdata[,mean(y1),j1][order(j1),V1] - rliml$A1)/rliml$B1

  # ---------- MOVERS: use covariance restrictions  --------------
  # we start by computing Var(Y1), Var(Y2) and Cov(Y1,Y2)
  setkey(jdata,j1,j2)
  YY1 = jdata[,var(y1),list(j1,j2)][,V1]
  YY2 = jdata[,cov(y1,y2),list(j1,j2)][,V1]
  YY3 = jdata[,var(y2),list(j1,j2)][,V1]
  XX1 = array(0,c(nf^2, nf^2 + 2*nf))
  XX2 = array(0,c(nf^2, nf^2 + 2*nf))
  XX3 = array(0,c(nf^2, nf^2 + 2*nf))
  
  for (l1 in 1:nf) for (l2 in 1:nf) {
    ll = l2 + nf*(l1 -1)
    XX1[ll,ll                 ] = B1[l1]^2
    XX1[ll,nf^2 + l1          ] = 1
    XX2[ll,ll                 ] = B1[l1]*B2[l2]
    XX3[ll,ll                 ] = B2[l2]^2
    XX3[ll,nf^2 + nf + l2     ] = 1
  }
  W = jdata[,.N,list(j1,j2)][,N]
  Wm = acast(jdata[,.N,list(j1,j2)],j1~j2)
  
  XX = rbind(XX1,XX2,XX3)
  C1 = diag(nf^2 + 2*nf)
  C0 = rep(0,nf^2 + 2*nf)
  
  res     = lm.wfitc( XX, c(YY1,YY2,YY3), c(W,W,W), C1,C0,0)$solution
  EEsd    = t(sqrt(pmax(array((res)[1:nf^2],c(nf,nf)),0)))
  eps1_sd = sqrt(pmax((res)[(nf^2+1):(nf^2+nf)],0))
  eps2_sd = sqrt(pmax((res)[(nf^2+nf+1):(nf^2+2*nf)],0))
  
  # ---------- STAYERS: use covariance restrictions  --------------
  setkey(sdata,j1)
  YY1 = sdata[,var(y1),list(j1)][,V1] - eps1_sd^2
  YY2 = sdata[,var(y2),list(j1)][,V1] - eps2_sd^2
  XX1 = diag(B1^2)
  XX2 = diag(B2^2)
  W = sdata[,.N,list(j1)][,N]
  
  XX = rbind(XX1,XX2)
  C1 = diag(nf)
  C0 = rep(0,nf)
  res2 = lm.wfitc( XX, c(YY1,YY2), c(W,W), C1,C0,0)$solution
  Esd = sqrt(pmax(res2,0))
  
  model = list(A1=rliml$A1,B1=B1,A2=rliml$A2,B2=B2,EEsd=EEsd,Em=Em,eps1_sd=eps1_sd,Esd = Esd,Ns=W,Nm=Wm)
  
  if (length(model0)>0) {
    rr = addmom(Em,model0$Em,"E(alpha|l1,l2)")
    rr = addmom(model$A1,model0$A1,"A1",rr)
    rr = addmom(model$A2,model0$A2,"A2",rr)
    rr = addmom(B1,model0$B1,"B1",rr)
    rr = addmom(B2,model0$B2,"B2",rr)
    rr = addmom(EEsd,model0$EEsd,"EEsd",rr)
    rr = addmom(eps1_sd,model0$eps_sd,"eps_sd",rr)
    rr = addmom(Esd,model0$Esd,"Esd",rr)
    print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    catf("cor with true model:%f",cor(rr$val1,rr$val2))
  }
  
  return(model)
}

model.show <- function(model) {
  
  # compute the cavariance terms
  
  v1 = wtd.var(model$A1, model$Ns)
  v2 = wtd.mean(model$eps1_sd^2, model$Ns)
  v3 = wtd.var(model$B1 * model$Em, model$Ns ) + wtd.mean(model$B1^2 * model$Esd^2, model$Ns  ) 
  v4 = wtd.var(rbind(model$A1,model$B1*model$Em),model$Ns )
  catf(" V(alpha)=%f \t V(psi)=%f \t 2*cov")
  
  
}



model.mini2.test <- function() {
  
  model = model.mini2.new(10,serial = F)
  NNm   = array(200000/(model$nf^2),c(model$nf,model$nf))
  NNs  = array(300000/model$nf,model$nf)
  jdata = model.mini2.simulate.movers(model,NNm)
  sdata = model.mini2.simulate.stayers(model,NNs)
  

  
  
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

# 
extract.variances <- function(jdata,sdata,b_l,rr=list(),na.rm=F) {

  # compute the variance of espilon conditional on l1, l2
  alpha_var_l1l2 = acast(jdata[, list(V1 = cov(y1,y2) /( b_l[j1] *b_l[j2])) , list(j1,j2)],j1~j2)
  
  # compute the variance of espilon conditional on l1, l2  
  espilon_var_l1l2 = jdata[, list( var(y1) - cov(y1,y2) * b_l[j1]/b_l[j2] ,.N) , list(j1,j2)]
  if (na.rm==TRUE) espilon_var_l1l2[ is.na(V1), V1:=0 ];
  
  # compute the variance of espilon conditional on l1, l2
  epsilon_var_l  = espilon_var_l1l2[ , weighted.mean(V1,N) , j1][, V1]
  
  # compute the variance of alpha conditional on l
  alpha_var_l = (sdata[,var(y),j][,V1] - epsilon_var_l)/b_l^2
  
  rr$espilon_var_l1l2 = acast(espilon_var_l1l2,j1~j2,value.var = "V1")
  rr$epsilon_var_l    = epsilon_var_l
  rr$alpha_var_l      = alpha_var_l
  rr$alpha_var_l1l2  = alpha_var_l1l2
  
  return(rr)
}

