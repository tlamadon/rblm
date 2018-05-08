catf <- function(...) cat(sprintf(...))

#' this is a utility function to generate
#' multidimensional arrays - like the spread function in fortran
#' @export
spread <- function (A, loc, dims) {
  if (!(is.array(A))) {
    A = array(A, dim = c(length(A)))
  }
  adims = dim(A)
  l = length(loc)
  if (max(loc) > length(dim(A)) + l) {
    stop("incorrect dimensions in spread")
  }
  sdim = c(dim(A), dims)
  edim = c()
  oi = 1
  ni = length(dim(A)) + 1
  for (i in c(1:(length(dim(A)) + l))) {
    if (i %in% loc) {
      edim = c(edim, ni)
      ni = ni + 1
    }
    else {
      edim = c(edim, oi)
      oi = oi + 1
    }
  }
  return(aperm(array(A, dim = sdim), edim))
}

hist2 <- function(Y1,Y2,wsup) {
  n = length(wsup)
  H = array(0,c(n,n))
  for (i in 1:n) for (j in 1:n) {
    H[i,j] = sum(  (Y1 < wsup[i]) & (Y2 < wsup[j]) )
  }
  H[n,n] = length(Y1)
  H = H/H[n,n]
  return(H)
}

hist1 <- function(Y1,wsup) {
  n = length(wsup)
  H = array(0,c(n))
  for (i in 1:n) {
    H[i] = sum(  (Y1 < wsup[i]) )
  }
  H[n] = length(Y1)
  H = H/H[n]
  return(H)
}

# allow to use 2 different supports
hist2d <- function(Y1,Y2,wsup) {
  n = length(wsup)
  H = array(0,c(n,n))
  for (i in 1:n) for (j in 1:n) {
    H[i,j] = sum(  (Y1 < wsup[i]) & (Y1 >= wsup[i-1]) & (Y2 < wsup[j]) & (Y1 >= wsup[i-1])  )
  }
  H[n,n] = length(Y1)
  H = H/H[n,n]
  return(H)
}

# smoothed histogram
hist2s <- function(Y1,Y2,wsup,h) {
  n = length(wsup)
  H = array(0,c(n,n))
  for (i in 1:n) for (j in 1:n) {
    H[i,j] = sum(  pnorm( (wsup[i] - Y1 )/h ) *  pnorm( (wsup[j] - Y2 )/h ) )
  }
  H = H/H[n,n]
  return(H)
}

#' @export
rdim <- function(A,...) {
  dd <- list(...);
  if (length(dd)==1) {
    dim(A)<-dd[[1]]
  } else {
    dim(A) <- dd
  }
  return(A)
}

tic.new <- function() {
  t = Sys.time()
  tt = list(all.start=t,last.time=t,loop=0,timers=list())

  tic.toc <- function(name="") {
    t = Sys.time()
    if (name=="") {
      return(tt)
    }

    if (name %in% names(tt$timers)) {
      tm = tt$timers[[name]]
      tm$count = tm$count+1
      tm$total = tm$total + t - tt$last.time
      tt$timers[[name]] = tm
    } else {
      tm = list(count=1,total=t - tt$last.time)
      tt$timers[[name]] = tm
    }
    tt$last.time=t;
    tt <<- tt
  }

  return(tic.toc)
}

#' order cluster by increasing wage
#' @export
cluster.order <- function(sim) {
  sim$sdata = sim$sdata[!is.na(j1)]
  I = sim$sdata[,mean(y1),j1][order(j1)][,rank(V1)]
  sim$sdata[,j1:=I[j1]][,j2:=I[j2]]
  sim$jdata[,j1:=I[j1]][,j2:=I[j2]]
  return(sim)
}

mvar <- function(x) {
  if (length(x)<=1) return(0);
  return(var(x))
}
mcov <- function(x,y) {
  if (length(x)<=1) return(0);
  return(cov(x,y))
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

# fits a weighted ols with non-negative constraints
lm.wfitnn <- function(XX,YY,rw,floor = 0) {

  n = dim(XX)[2]
  XX2 = XX
  #   S = apply(abs(XX),2,max)
  #   XX2 = XX*spread(1/S,1,dim(XX)[1])
  #   C12 = C1*spread(1/S,1,dim(C1)[1])

  XXw      = diag(rw) %*% XX2
  Dq       = t(XXw) %*% XX2
  dq       = t(YY %*% XXw)
  C1       = diag(n)
  C0       = rep(floor,n)

  #fit      = qprog(Dq,dq,C1,C0)
  #fit$solution = as.numeric(fit$thetahat)
  fit      = solve.QP(Dq,dq,C1,C0)

  return(fit)
}

# fits a linear problem with weights under constraints
slm.wfitc <- function(XX,YY,rw,CS,scaling=0) {
  nk = CS$nk
  nf = CS$nf
  YY = as.numeric(YY)
  XX = as.matrix.csr(XX)
  # to make sure the problem is positive semi definite, we add
  # the equality constraints to the XX matrix! nice, no?

  if (CS$meq>0) {
    XXb = rbind(XX,  as.matrix.csr(CS$C[1:CS$meq,]))
    YYb = c(YY,CS$H[1:CS$meq])
    rwb  = c(rw,rep(1,CS$meq))
  } else {
    XXb = XX
    YYb = YY
    rwb = rw
  }

  t2 = as(dim(XXb)[1],"matrix.diag.csr")
  t2@ra = rwb
  XXw = t2 %*% XXb
  Dq       = SparseM:::as.matrix(SparseM:::t(XXw) %*% XXb)
  dq       = SparseM:::t(YYb %*% XXw)

  # scaling
  #
  if (scaling>0) {
    sc <- norm(Dq,"2")^scaling
  } else {
    sc=1
  }

  # do quadprod
  tryCatch({
    fit      = solve.QP(Dq/sc,dq/sc,t(CS$C)/sc,CS$H/sc,CS$meq)
  }, error = function(err) {
    browser()
  })

  return(fit)
}

#' Computes graph connectedness among the movers
#' within each type and returns the smalless value
#' @export
model.connectiveness <- function(model,all=FALSE) {
  EV = rep(0,model$nk)
  pk1 = rdim(model$pk1,model$nf,model$nf,model$nk)
  dd_post = data.table(melt(pk1,c('j1','j2','k')))
  pp = model$NNm/sum(model$NNm)
  dd_post <- dd_post[, pr_j1j2 := pp[j1,j2],list(j1,j2)  ]
  dd_post <- dd_post[, pr_j1j2k := pr_j1j2*value]

  for (kk in 1:model$nk) {
    # compute adjency matrix
    A1 = acast(dd_post[k==kk, list(pr=pr_j1j2k/sum(pr_j1j2k),j2,j1)],j1~j2,value.var = "pr")
    A2 = acast(dd_post[k==kk, list(pr=pr_j1j2k/sum(pr_j1j2k),j2,j1)],j2~j1,value.var = "pr")
    # construct Laplacian
    A = 0.5*A1 + 0.5*A2
    D = diag( rowSums(A)^(-0.5) )
    L = diag(model$nf) - D%*%A%*%D
    EV[kk] = sort(eigen(L)$values)[2]
    #print(eigen(L)$values)
  }
  if (all==TRUE) return(EV);
  return(min(abs(EV)))
}

#' plots the wages of a model
#' @export
m2.mixt.wplot <- function(Wm) {
  dd = melt(Wm,c('l','k'))
  ggplot(dd,aes(x=factor(l),color=factor(k),group=factor(k),y=value)) + geom_line() + theme_bw()
}

#' @export
mplot <- function(M) {
  mm = melt(M,c('i','j'))
  #mm$i = factor(mm$i)
  #mm$j = factor(mm$j)
  mm = mm[mm$value>0,]
  ggplot(mm,aes(x=j,y=i,fill=value)) + geom_tile() + theme_bw() + scale_y_reverse()
}

#' plots the proportions of a model
#' @export
m2.mixt.pplot <- function(pk0) {
  dd = melt(pk0,c('l','k'))
  ggplot(dd,aes(x=factor(l),y=value,fill=factor(k))) + geom_bar(position="stack",stat = "identity") + theme_bw()
}
