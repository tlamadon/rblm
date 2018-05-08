# ------------    CONSTRAINTS ---------------

# imposese that a_k1_l1 - a_k2_l1 = a_k1_l2 - a_k2_l2 for all k1,k2,l1,l2
cons.lin_add <- function(nk,nf) {
  LL = array(0,c(nf-1,nf))
  for (l in 1:(nf-1)) { LL[l,l]=1; LL[l,l+1]=-1}
  KK = array(0,c(nk-1,nk))
  for (k in 1:(nk-1)) { KK[k,k]=1; KK[k,k+1]=-1}
  C1     = kronecker(LL,KK)
  H1      = rep(0,dim(C1)[1])
  meq     = dim(C1)[1]
  return(list(C=C1,H=H1,meq=meq,nk=nk,nf=nf))
}

cons.lin <- function(nk,nf,gap=0) {
  KK = array(0,c(nk-1,nk))
  for (k in 1:(nk-1)) { KK[k,k]=1; KK[k,k+1]=-1}
  
  LL = array(0,c(nf-1,nf))
  for (l in 1:(nf-1)) { LL[l,l]=1; LL[l,l+1]=-1}
  C1    = kronecker(LL,KK)
  H1    = rep(0,dim(C1)[1])
  meq   = dim(C1)[1]
  
  return(list(C=C1,H=H1,meq=meq,nk=nk,nf=nf))
}


cons.akm <- function(nk,nf,gap=0) {
  KK = array(0,c(nk-1,nk)) 
  LL = array(0,c(nf-1,nf))
  for (k in 1:(nk-1)) { KK[k,k]=1; KK[k,k+1]=-1}
  for (l in 1:(nf-1)) { LL[l,l]=1; LL[l,l+1]=-1}
  
  C1      = kronecker(LL,KK)
  H1      = rep(0,dim(C1)[1])
  meq     = dim(C1)[1]
  
  return(list(C=C1,H=H1,meq=meq,nk=nk,nf=nf))
}


cons.akmmono <- function(nk,nf,gap=0) {
  KK = array(0,c(nk-1,nk))
  for (k in 1:(nk-1)) { KK[k,k]=1; KK[k,k+1]=-1}
  C1     = -kronecker(diag(nf),KK)
  H1     = rep(gap,nf*(nk-1))
  meq    = 0
  
  LL = array(0,c(nf-1,nf))
  for (l in 1:(nf-1)) { LL[l,l]=1; LL[l,l+1]=-1}
  C1b     = kronecker(LL,KK)
  C1      = rbind(C1b,C1)
  H1      = c(rep(0,dim(C1b)[1]),H1)
  meq     = dim(C1b)[1]
  
  return(list(C=C1,H=H1,meq=meq,nk=nk,nf=nf))
}

cons.mono_k <- function(nk,nf,gap=0) {
  KK = array(0,c(nk-1,nk))
  for (k in 1:(nk-1)) { KK[k,k]=1; KK[k,k+1]=-1}
  C1     = -kronecker(diag(nf),KK)
  H1     = rep(gap,nf*(nk-1))
  meq    = 0
  return(list(C=C1,H=H1,meq=meq,nk=nk,nf=nf))
}

cons.fixb <- function(nk,nf,nt=4) {
  CC = cons.mono_k(nk,nf)$C
  MM = array( 0,c(nt-1,nt))
  for (i in 1:(nt-1)) {
    MM[i,i]=1
    MM[i,i+1]=-1
  }
  CC  = kronecker(MM,CC)
  H1  = rep(0, (nt-1)*(nk-1)*nf )
  meq = (nt-1)*(nk-1)*nf 
  return(list(C=CC,H=H1,meq=meq,nk=nk,nf=nf))
}

cons.biggerthan <-function(nk,nf,gap=0) {
  CC = diag(nk*nf)
  H1 = rep(gap,nk*nf)
  meq=0
  return(list(C=CC,H=H1,meq=meq,nk=nk,nf=nf))
}

cons.lin_para <- function(nk,nf) {
  LL = array(0,c(nf-1,nf))
  for (l in 1:(nf-1)) { LL[l,l]=1; LL[l,l+1]=-1}
  C1     = -kronecker(LL,diag(nk))
  H1     = rep(0,(nf-1)*nk)
  meq    = (nf-1)*nk
  return(list(C=C1,H=H1,meq=meq,nk=nk,nf=nf))
}

cons.none <- function(nk,nf) {
  meq=0
  C1 = matrix(0,1,nk*nf)
  H1 = c(0)
  return(list(C=C1,H=H1,meq=meq,nk=nk,nf=nf))
}

cons.bind <- function(c1,c2) {
  
  if (c1$meq==0) {
    I11 = c(); I12 = 1:length(c1$H);
  } else if (c1$meq==length(c1$H)) {
    I11 = 1:length(c1$H); I12 = c();
  } else {
    I11 = 1:c1$meq; I12 = (c1$meq+1):length(c1$H)
  }
  if (c2$meq==0) {
    I21 = c(); I22 = 1:length(c2$H);
  } else if (c2$meq==length(c2$H)) {
    I21 = 1:length(c2$H); I22 = c();
  } else {
    I21 = 1:c2$meq; I22 = (c2$meq+1):length(c2$H)
  }
  
  c1$C = rBind(c1$C[I11,],c2$C[I21,],c1$C[I12,],c2$C[I22,])
  c1$H = c(c1$H[I11],c2$H[I21],c1$H[I12],c2$H[I22])
  c1$meq = c1$meq + c2$meq
  return(c1)
}

# add right padding
cons.pad <- function(c1,l,r) {
  c1$C = cBind(matrix(0,dim(c1$C)[1],l),c1$C,matrix(0,dim(c1$C)[1],r))
  return(c1)
}

cons.get <-function(name,val,nk,nf) {
  if (name=="none") {
    return(cons.none(nk,nf))
  } else if (name=="mono_k") {
    return(cons.mono_k(nk,nf,val))
  } else if (name=="para") {
    return(cons.lin_para(nk,nf))
  } else if (name=="akm") {
    return(cons.akm(nk,nf,val))
  } else if (name=="akmmono") {
    return(cons.akmmono(nk,nf,val))
  } else if (name=="lin") {
    return(cons.lin(nk,nf))
  } else {
    error("unkown constraint")
  }
}

cons.sum <- function(nk,nf) {
  CC = kronecker(diag(nf),t(rep(1,nk)))
  C1     = kronecker(diag(nf),t(rep(1,nk)))
  H1     = rep(0,nf)
  meq    = nf
  return(list(C=C1,H=H1,meq=meq,nk=nk,nf=nf))
}
