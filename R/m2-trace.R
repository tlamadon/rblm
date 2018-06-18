#' some work on trace formula

#' Extracts the largest connected set from data using f1,f2
#' movers
#' @export
get.largest.conset.fid <- function(jdata) {
  # combine th 2 directions
  jdata2 = rBind(jdata[,list(f1,f2)],jdata[,list(f1=f2,f2=f1)])
  AD = pmin(acast(jdata2[,.N,list(f1,f2)],value.var = "N",f1~f2,drop=FALSE,fill=0),1)
  # compute connected sets
  cs = conComp(AD,2)
  # extract largest
  cs = data.table(f1=names(cs),set=as.numeric(cs))
  cs.size = cs[,.N,set][order(-N)]
  return(cs[  set == cs.size$set[1]  , f1])
}


#' create a model for testing trace estimation
#' @export
m2.trace.new <- function(nf = 200 , nm = 10000, eps_sd = 1.5) {
  p = list(dsize = c(exp(3.199292), 1+exp(-1.247662)) , nf = nf , nm = nm, eps_sd = eps_sd  )
  S = rpareto(p$nf,p$dsize[1],p$dsize[2])
  psi = rnorm(p$nf)

  model     = copy(p)
  model$S   = S
  model$psi = psi

  return(model)
}

#' simulates for trace estimation
#' @export
m2.trace.simulate.old <- function(model) {
  JJ = array(0,c(model$nm,model$nf))
  AA = array(0,c(model$nf,model$nf))
  F1 = rep(0,model$nm)
  F2 = rep(0,model$nm)
  for (i in 1:model$nm) {
    ii = sample.int(model$nf,2,prob = model$S)
    JJ[i,ii[1]] = 1
    JJ[i,ii[2]] = -1
    AA[ii[1],ii[2]]=1
    AA[ii[2],ii[1]]=1
    F1[i] = ii[1]
    F2[i] = ii[2]
  }
  colnames(AA) = 1:model$nf
  rownames(AA) = 1:model$nf

  D = JJ %*% model$psi + rnorm(model$nm)*model$eps_sd

  return(list(JJ=JJ,AA=AA,D=D,F1=F1,F2=F2))
}

#' simulates for trace estimation
#' @export
m2.trace.simulate <- function(model) {
  F1   = rep(0,model$nm)
  F2   = rep(0,model$nm)
  psi1 = rep(0,model$nm)
  psi2 = rep(0,model$nm)
  for (i in 1:model$nm) {
    ii = sample.int(model$nf,2,prob = model$S)
    F1[i] = ii[1]
    F2[i] = ii[2]
    psi1[i] = model$psi[ii[1]]
    psi2[i] = model$psi[ii[2]]
  }

  jdata = data.table(f1=paste(F1),f2=paste(F2),psi1=psi1,psi2=psi2)
  jdata[, y1 := psi1 + rnorm(.N)*model$eps_sd/sqrt(2)]
  jdata[, y2 := psi2 + rnorm(.N)*model$eps_sd/sqrt(2)]

  sdata = data.table(f1=paste(1:model$nf),psi1=model$psi,size=model$S)
  sdata = sdata[,list(y1 = psi1 + rnorm(ceiling(size))*model$eps_sd/sqrt(2),
                      y2 = psi1 + rnorm(ceiling(size))*model$eps_sd/sqrt(2)), list(f1,size,psi1) ]
  sdata$x=1

  return(list(sdata=sdata,jdata=jdata))
}



#' gets the connected set, then
#' @export
m2.trace.estimate <- function(sim, model0=NA,hetero=FALSE) {

  # EXTRACT CONNECTED SET
  jdata = sim$jdata
  f1s   = get.largest.conset.fid(jdata)
  flog.info("connected set %i/%i",length(f1s),length(unique(sim$sdata$f1)))

  jdata = jdata[f1 %in% f1s][f2 %in% f1s]

  # index firms with integers
  fids = data.table(f1=f1s,nfid=1:length(f1s))
  setkey(fids,f1)
  setkey(jdata,f1)
  jdata[,f1i := fids[jdata,nfid]]
  setkey(jdata,f2)
  jdata[,f2i := fids[jdata,nfid]]

  # extract size in stayers, and mean
  fsize = rbind(sim$sdata[,list(y1,f1)],sim$jdata[,list(y1,f1)])[,list(N=.N,mw=mean(y1)),f1]
  setkey(fids,f1)
  setkey(fsize,f1)
  fids[,size := fsize[fids,N]]
  fids[,mw   := fsize[fids,mw]]
  fids[is.na(size),size:=0] # pads missing size with 0

  # CONSTRUCT SPARSE DESIGN MATRIX + NORMALIZATION
  nf = pmax(max(jdata$f1i),max(jdata$f2i))
  dd = jdata[,list(m1 = y2 - y1,f1i,f2i)]
  dd = rbind(dd,data.frame(m1=0,f1i=1,f2i=1))
  N = nrow(dd)
  dd[,v1:=-1][,v2:=1]
  dd[,c1:=f1i][,c2:=f2i]
  dd[N,v1:=0]
  JJ = sparseMatrix2(1:N,dd$c1,dd$v1,c(N,nf)) + sparseMatrix2(1:N,dd$c2,dd$v2,c(N,nf))
  S  = fids[order(nfid),size]

  # COMPUTE INVERSE OF DESIGN MATRIX
  M = SparseM::t(JJ) %*% JJ
  Minv = SparseM::solve(M,nnzlmax=1e8,tmpmax=1e8,nsubmax=1e8)

  # compute firms FE
  psi = as.numeric(SparseM::as.matrix(Minv %*% ( SparseM::t(JJ) %*% dd$m1)))
  E   = dd$m1 - JJ%*%psi
  E   = E[1:(N-1)]

  if (!any(is.na(model0))) {
    psi0 = model0$psi[as.integer(fids[order(nfid),f1])]
    psi0 = psi0-psi0[1]
    flog.info("corr= %f", cor(psi0,psi))
  }

  # extract homoskedastic error
  ff = nrow(sim$jdata)/(nrow(sim$jdata)-nf)
  var_e = var(E)*nrow(sim$jdata)/(nrow(sim$jdata)-nf)


  # COMPUTE THE TRACE FORMULA - USING THE WEIGHTING OF STAYERS
  tr_correction = var_e*( sum( SparseM::diag(Minv)*S )/sum(S) - sum( S* (Minv %*% rep(1/nf,nf)) )/(sum(S))  )

  # heteroskedastic variance using BLM groups
  if (hetero==TRUE) {
    jdata[, psi1 := psi[f1i]]
    jdata[, psi2 := psi[f2i]]

    R = array(0,c(nf,nf))
    for (i1 in 1:10) for (i2 in 1:10) {
      I = jdata[,list(V1=1:.N,j1,j2)][j1==i1 & j2==i2,V1] # select the rows associated with given (j1,j2)if
      if (length(I)==0) next;
      if (length(I)<5) {
        R = R + var_e*SparseM::t(JJ[I,]) %*% JJ[I,]
      } else {
        n1 = length(I)
        n2 = length(jdata[j1==i1 & j2==i2,unique(c(f1,f2))])
        #flog.info("%i %i",n1,n2)
        R = R + ff*jdata[j1==i1 & j2==i2,var( y2-psi2 - y1 + psi1)]*SparseM::t(JJ[I,]) %*% JJ[I,]
      }
    }

    M2 = Minv %*% R %*% Minv
    tr_correction =  sum( SparseM::diag(M2)*S )/sum(S) - sum( S* (M2 %*% rep(1/nf,nf)) )/(sum(S))
  }

  # COMPUTE THE VARIANCE OF PSI
  fids[, psi := psi[nfid]]
  var_psi_hat = fids[,wt.var(psi,size)]
  tot = sim$sdata[,var(y1)]
  btw = sim$sdata[,list(mean(y1),.N),f1][,wt.var(V1,N)]

  if (!any(is.na(model0))) {
    fids[, psi0 := psi0[nfid]]
    flog.info("var_true=%f  var_akm=%f var2=%f trace=%f ", fids[,wt.var(psi0,S)], fids[,wt.var(psi,S)], fids[,wt.var(psi,S)]- tr_correction,tr_correction)
  } else {
    flog.info("tot=%f btwf=%f var_akm=%f var2=%f trace=%f ",tot,btw, fids[,wt.var(psi,S)], fids[,wt.var(psi,S)]- tr_correction,tr_correction)
  }

  res = list(fids=fids,eps_sd = sqrt(var_e), trace_correction = tr_correction, var_psi= var_psi_hat,tot=tot,btw=btw )
}



#' Ridge AKM
#' @export
m2.firmfe.pen  <- function(sim, model, lambda=1,holdout=0.1) {

  # index firms with integers
  jdata = sim$jdata
  f1s   = unique(c(sim$jdata[,unique(f1,f2)],unique(sim$sdata[,unique(f1)] )))
  fids  = data.table(f1=f1s,nfid=1:length(f1s))
  setkey(fids,f1)
  setkey(jdata,f1)
  jdata[,f1i := fids[jdata,nfid]]
  setkey(jdata,f2)
  jdata[,f2i := fids[jdata,nfid]]

  # extract size in stayers
  fsize = rbind(sim$sdata[,list(y1,f1,j1)],sim$jdata[,list(y1,f1,j1)])[,.N,list(f1,j1)]
  setkey(fids,f1)
  setkey(fsize,f1)
  fids[,size := fsize[fids,N]]
  fids[is.na(size),size:=0] # pads missing size with 0
  fids[,j1 := fsize[fids,j1]]

  # create the matrix
  nf = pmax(max(jdata$f1i),max(jdata$f2i))
  dd = jdata[,list(m1 = y2 - y1,f1i,f2i)]
  N = nrow(dd)
  dd[,v1:=-1][,v2:=1]
  dd[,c1:=f1i][,c2:=f2i]
  JJ = sparseMatrix2(1:N,dd$c1,dd$v1,c(N,nf)) + sparseMatrix2(1:N,dd$c2,dd$v2,c(N,nf))

  # get the vector of sizes
  S  = fids[order(nfid),size]

  # create the penalty, just append the values
  R  = as.numeric(model$A1)[fids[order(nfid),j1]]
  RJ = as(length(R),"matrix.diag.csr")

  # we create a hold out (by setting a random set of obsevartions to 0)
  I = sample.int(N,ceiling(holdout*N),replace = FALSE)

  XX = rbind(JJ,RJ)
  YY = c(dd$m1,R)
  WW1 = rep(1,nrow(dd))
  WW1[I] = 0
  WW = c( WW1, lambda * rep(1,length(R)) )
  fit = SparseM::slm.wfit(XX,YY,WW,nnzlmax=1e8,tmpmax=1e8,nsubmax=1e8)

  # compute the prediction error on holdout
  MSE = mean( ( dd$m1[I]- JJ[I,] %*% fit$coefficients)^2)

  fids[, psi := fit$coefficients[nfid]]
  flog.info(" Ridge AKM lambda=%f var(psi)=%f mse=%f",lambda,fids[,wt.var(psi,size)]/sim$sdata[,var(y1)],MSE)
  list(fids=fids,mse=MSE)
}

check.data <- function(sim) {
  sf1 = sim$sdata[,unique(f1)]
  jf1 = sim$jdata[,unique(f1)]
  jf2 = sim$jdata[,unique(f2)]

  # compute firms in sdata but not in jdata
  flog.info(" %i firms in sdata but not in jdata",   length(setdiff( sf1, c(jf1,jf2) )  ))
  flog.info(" %i firms in jdata.f1 but not in sdata",length(setdiff( jf1, sf1 )))
  flog.info(" %i firms in jdata.f2 but not in sdata",length(setdiff( jf2, sf1 )))
  flog.info(" %i firms in p2 but not in p1",length(setdiff( jf2, c(sf1,jf1 ))))
}



test.m2.trace <- function() {
  require(rmutil)

  model   = m2.trace.new(nf = 1000,nm = 3000,eps_sd=1.5)
  sim     = m2.trace.simulate(model)

  # cluster
  ms    = grouping.getMeasures(sim,"ecdf",Nw=20,y_var = "y1")
  grps  = grouping.classify.once(ms,k = 5,nstart = 1000,iter.max = 200,step=100)
  sim   = grouping.append(sim,grps$best_cluster)

  res     = m2.trace.estimate(sim,model0=model,hetero=T)
  res     = m2.trace.estimate(sim,model0=model,hetero=F)

  # estimate BLM
  res     = m2.mini.estimate(sim$jdata,sim$sdata,method = "linear.ss")

  # estimate RIDGE AKM
  rr = m2.firmfe.pen(sim,res,lambda=10)

  # check if things match at all
  ff = sim$sdata[,list(psi=psi1[1],N=.N),list(f1,j1)]
  ff[, psi_g := res$A1[j1],j1]
  setkey(ff,f1)
  setkey(rr,f1)
  ff[, psi_hat := rr[ff,psi]]
  ff[,cov(psi,psi_hat)]
  ff[,cov(psi,psi_g)]

  ff[,wt.var(psi,N)]
  ff[,wt.var(psi_g,N)]
  ff[,wt.var(psi_hat,N)]

  II = sim$jdata[,list(f1=c(f1,f2))][,.N,f1][N %in% c(2,3,4),f1]

  dd = data.frame()
  for (lambda in 10^(seq(-6,6,l=10))) {
    akm.ridge = m2.firmfe.pen(sim,res,lambda=lambda,holdout = 0.4)
    rr = akm.ridge$fids
    ff = sim$sdata[,list(psi=psi1[1],N=.N),list(f1,j1)]
    ff[, psi_g := res$A1[j1],j1]
    setkey(ff,f1)
    setkey(rr,f1)
    ff[, psi_hat := rr[ff,psi]]

    flog.info("var0=%f lambda=%f var_psi=%f mse0=%f",ff[,wt.var(psi,N)],lambda,ff[,wt.var(psi_hat,N)],ff[II,wt.mean((psi-psi_hat)^2,N)])
    dd = rbind(dd,data.frame(lambda=log(lambda)/log(10), var0=ff[,wt.var(psi,N)],var_ridge=ff[,wt.var(psi_hat,N)],mse=akm.ridge$mse))
  }

}

