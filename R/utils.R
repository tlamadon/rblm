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

kmeansW.repeat <- function(x, centers, weight = rep(1, nrow(x)), iter.max = 10, nstart = 1,step=20) {
  
  cat(sprintf("running weigthed kmeans step=%i total=%i\n nobs=%i nmeasures=%i \n", step, nstart, dim(x)[1], dim(x)[2] ))
  
  best = Inf
  best_clus = NA
  
  for (i in 1:round(nstart/step)) {
    clus = kmeansW(x,centers,weight,iter.max,nstart=step)
    
    tot = sum(clus$withinss)
    if (tot < best) {
      best_clus = clus
      best = tot
      cat(sprintf("[%2i%%] tot=%f best=%f <<<< \n", round(100*i*step/nstart), tot,  best))
    } else {
      cat(sprintf("[%2i%%] tot=%f best=%f \n", round(100*i*step/nstart), tot,  best))
    }
        
  }
   
  return(clus)
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





