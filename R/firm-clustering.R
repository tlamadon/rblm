#' extract information
#' @export
grouping.infos <- function(cdata) {

  cdata[, age:= min(aret) - birthyear]

  rr = cdata[(fid!="") & (aret==min(aret)),{
    r = list()
    r$ni = .N
    r$nj = length(unique(fid))

    r$educ1 = .SD[educ==1,.N]
    r$educ2 = .SD[educ==2,.N]
    r$educ3 = .SD[educ==3,.N]

    r$worker_0_30   = .SD[  (age <= 30) ,.N]
    r$worker_30_50  = .SD[  (age > 30) & (age <=50) ,.N]
    r$worker_50_100 = .SD[  (age >50) ,.N]

    r$female = .SD[female=="female",.N]
    r$cohort    = .SD[, mean(birthyear)]

    r$mw     = mean(lw)
    r$wsd    = var(lw)
    r$bwfsd  = .SD[,rep(mean(lw),.N),fid][,var(V1)]
    r$bwcsd  = 0

    r$bwfq50 = .SD[, rep(quantile(lw,0.5),.N), fid][,var(V1)]
    r$bwcq50 = 0
    r$bwfq10 = .SD[, rep(quantile(lw,0.1),.N), fid][,var(V1)]
    r$bwcq10 = 0
    r$bwfq90 = .SD[, rep(quantile(lw,0.9),.N), fid][,var(V1)]
    r$bwcq90 = 0

    r$bwfvar = .SD[, rep(var(lw),.N), fid][,var(V1,na.rm=T)]
    r$bwcvar = 0

    r$pm      = .SD[,list(prod=prod[1]),fid][,mean(log(prod),na.rm=T)]
    r$psd     = .SD[,list(prod=prod[1]),fid][,var(log(prod),na.rm=T)]

    r$sizem    = .SD[,.N,fid][,mean(N)]
    r$sizem2   = .SD[,rep(.N,.N),fid][,mean(V1)]
    r$sizemed  = .SD[,.N,fid][,median(N)+0.0]
    r$sizemed2 = .SD[,rep(.N,.N),fid][,median(V1)+0.0]

    r
  }, clus]

  rr2 = cdata[(fid!="") & (aret==min(aret)),{
    r = list()
    r$ni = .N
    r$nj = length(unique(fid))

    r$educ1 = .SD[educ==1,.N]
    r$educ2 = .SD[educ==2,.N]
    r$educ3 = .SD[educ==3,.N]

    r$worker_0_30   = .SD[  (age <= 30) ,.N]
    r$worker_30_50  = .SD[  (age > 30) & (age <=50) ,.N]
    r$worker_50_100 = .SD[  (age >50) ,.N]

    r$female = .SD[female=="female",.N]
    r$cohort    = .SD[, mean(birthyear)]

    r$mw      = mean(lw)
    r$wsd     = var(lw)
    r$bwfsd   = .SD[,rep(mean(lw),.N),fid][, var(V1)]
    r$bwcsd   = .SD[,rep(mean(lw),.N),clus][,var(V1)]

    r$bwfq50 = .SD[, rep(quantile(lw,0.5),.N), fid][,var(V1)]
    r$bwcq50 = .SD[, rep(quantile(lw,0.5),.N), clus][,var(V1)]
    r$bwfq10 = .SD[, rep(quantile(lw,0.1),.N), fid][,var(V1)]
    r$bwcq10 = .SD[, rep(quantile(lw,0.1),.N), clus][,var(V1)]
    r$bwfq90 = .SD[, rep(quantile(lw,0.9),.N), fid][,var(V1)]
    r$bwcq90 = .SD[, rep(quantile(lw,0.9),.N), clus][,var(V1)]

    r$bwfvar = .SD[, rep(var(lw),.N), fid][,var(V1,na.rm=T)]
    r$bwcvar = .SD[, rep(var(lw),.N), clus][,var(V1,na.rm=T)]

    r$pm      = .SD[,list(prod=prod[1]),fid][,mean(log(prod),na.rm=T)]
    r$psd     = .SD[,list(prod=prod[1]),fid][,sd(log(prod),na.rm=T)]

    r$sizem    = .SD[,.N,fid][,mean(N)]
    r$sizem2   = .SD[,rep(.N,.N),fid][,mean(V1)]
    r$sizemed  = .SD[,.N,fid][,median(N)+0.0]
    r$sizemed2 = .SD[,rep(.N,.N),fid][,median(V1)+0.0]

    r
  }]

  rr2$clus=0
  rr = rbind(rr,rr2)
  setkey(rr,clus)

  return(rr)
}

#' internal function that runs Kmean with multiple starting values
kmeansW.repeat <- function(x, centers, weight = rep(1, nrow(x)), iter.max = 100, nstart = 1,step=20) {

  flog.info("running weigthed kmeans step=%i total=%i",step, nstart)
  flog.info("nobs=%i nmeasures=%i", dim(x)[1], dim(x)[2] )

  best = Inf
  best_clus = NA
  allres = data.frame()

  for (i in 1:round(nstart/step)) {

    clus = kmeansW(x,centers,weight,iter.max,nstart=step)
    tot = sum(clus$withinss)
    allres = rbind(allres,data.frame(i=i,tot=tot,nstrat=step))
    if (tot < best) {
      best_clus = clus
      best = tot
      flog.info("[%2i%%] tot=%f best=%f <<<<", round(100*i*step/nstart), tot,  best)
    } else {
      flog.info("[%2i%%] tot=%f best=%f", round(100*i*step/nstart), tot,  best)
    }
  }
  return(list(best=best_clus,all=allres))
}

#' Extract the measurement matrix to be given to the
#' classification algorithm
#' @export
grouping.getMeasures <- function(sim,y_var="y1",measure="ecdf",Nw=40,model=NA,exclude_movers=FALSE,user_fn=NA) {

  sdata = sim$sdata
  if (exclude_movers) {
    wids=sim$jdata[,wid]
    sdata = sdata[!(wid %in% wids)]
  }

  # extract dep variable
  sdata[,ydep := get(y_var)]

  N = sdata[,length(unique(f1))]
  qs = seq(1/(Nw+1),Nw/(Nw+1),l=Nw)
  flog.info("processing %i firms",N)

  # fix quantiles for wages
  wv    = as.numeric(sdata[,quantile(ydep,probs = qs)])
  t.cum = function(Y,W) {
    nj = length(W)
    R = rep(0,nj)
    for (j in 1:nj) {R[j] = sum(Y<=W[j])}
    return(R)
  }
  t.discrepancy = function(Y,W,E) {
    nj = length(W)
    R = rep(0,nj)
    for (j in 1:nj) {R[j] = sum(  ( (Y<=W[j]) - E[j])^2) }
    return(R)
  }

  rd=0
  # from here we have to be careful with the ordering of the firms
  # in the weighting matrix and in the measurement matrix

  flog.info("computing measures...")
  # for each firm get the empirical cdf
  if (str_detect(measure,"user_")) {
    rs = sdata[,user_fn(ydep),f1]
    M      = acast(rs,f1~m,value.var = "value")/sqrt(Nw)
  } else if (measure == "quant") {
    rs     = sdata[,list(cdf = quantile(ydep,probs=qs), s = 1:Nw,N=.N),f1]
    M      = acast(rs,f1~s,value.var = "cdf")/sqrt(Nw)
  } else if (measure == "ecdf") {
    # for each firm get the quantiles
    rs     = sdata[,list(cdf = as.numeric(t.cum(ydep,wv)/.N), s = 1:Nw),f1]
    M      = acast(rs,f1~s,value.var = "cdf")/sqrt(Nw)
    #rd     = sdata[, W[f1]*sum(t.discrepancy(ydep,wv,sqrt(Nw)*M[f1,]))/(.N^2),f1][,sum(V1)]/(N*Nw)
  } else if (measure == "ecdf_true") {
    # for each firm get the quantiles
    rs     = sdata[,list(cdf = pnorm(wv,mean= am[1] + psi[1], sd= sqrt( (model$eta_sd_m0 + model$eta_sd_m1*am[1])^2 + model$eps_sd^2)), s = 1:Nw),f1]
    M      = acast(rs,f1~s,value.var = "cdf")/sqrt(Nw)
    rd     = 0
  } else if (measure == "meanvar") {
    # use means instead
    rs_mean = sdata[,list(mean(ydep),var(ydep)),f1]
    M      = rs_mean[,list(V1,V2)]
  } else if (measure == "meanvar_true") {
    # use means instead
    rs_mean = sdata[,list( model$v_m0 + model$v_m1*am[1] + psi[1] ,(model$eta_sd_m0 + model$eta_sd_m1*am[1]))^2 + model$eps_sd^2 ,f1]
    M       = rs_mean[,list(V1,V2)]
  } else {
    stop("unkown measure")
  }

  flog.info("computing weights...")
  # extract weights, in the same order as the M matrix
  fids = rownames(M)
  W   = rep(0,length(fids))
  sdata.size = sdata[,.N,f1]
  sdata.size[,f1s:=paste(f1)]
  setkey(sdata.size,f1s)
  W = sdata.size[fids,N]

  # for (i in 1:length(fids)) {
  #   W[i] = sdata.size[f1==fids[i],N]
  # }

  # rescale the weights
  # W      = W/sqrt(sum(W^2))

  tsize = sum(W^2)/sum(W)
  if (is.na(rd)) rd=0;

  sdata[,ydep:=NULL]

  return(list(M=M,W=W,discrepency=rd,measure=measure,Nw=Nw,tsize=tsize,N=N,fids=rownames(M),rs=rs,wage.quant=wv))
}

#' clusters firms based on their cross-sectional wage distributions
#'
#' @param sdata cross sectional data, needs a column j (firm id) and w (log wage)
#' @param Nw number of points to use for wage distributionsdsd
#' @param ksupp vector of different number of groups to try
#' @param nstart (default:1000) total number of starting values
#' @param iter.max (default:100) max nunmber of step for each repetition
#' @param M  you can pass the matrix measurements, requires also weights W (pass on the truth for instance)
#' @param measures specify the type of measures to use (mean and var, quantiles, etc...)
#' @export
grouping.classify <- function(measures,ksupp = ceiling( (1:(60)^(1/1.3))^1.3),nstart=1000,iter.max=200,stop=FALSE,verbose=1,cval=1) {

  N = measures$N
  flog.info("clustering T=%f, Nw=%i , measure=%s",measures$tsize,measures$Nw,measures$measure);
  ksupp = setdiff(ksupp,1)

  rk = data.frame()
  rk_all = list()
  # we evaluate all values
  if (stop==FALSE) {
    for (k in ksupp) {
      if (k>nrow(measures$M)) next;  # go to next element if not enough row
      if (k<2) next;
      kres = kmeansW(measures$M,centers = k,weight=measures$W,nstart = nstart,iter.max = iter.max)
      rk=rbind(rk,data.frame(k=k,Q=sum(kres$withinss/N)))
      flog.info("k=%i WSS=%f nstart=%i nfrims=%i",k,sum(kres$withinss),nstart,N)
      rk_all[[k]] = kres$cluster
      if ((stop==TRUE) & (sum(kres$withinss/N) < cval*measures$discrepency)) break;
    }

    # search more efficiently
  } else {
    k_max = pmin(max(ksupp),nrow(measures$M)-1)
    k_min = pmax(min(ksupp),2)

    #evaluate k_max
    kres = kmeansW(measures$M,centers = k_max,weight=measures$W,nstart = nstart,iter.max = iter.max)
    flog.info("k=%i WSS=%f nstart=%i nfrims=%i",k,sum(kres$withinss),nstart,N)
    rk=rbind(rk,data.frame(k=k_max,Q=sum(kres$withinss/N)))
    kres_max = kres; rk_all[[k_max]] = kres$cluster
    #evaluate K_min
    kres = kmeansW(measures$M,centers = k_min,weight=measures$W,nstart = nstart,iter.max = iter.max)
    flog.info("k=%i WSS=%f nstart=%i nfrims=%i",k,sum(kres$withinss),nstart,N)
    rk=rbind(rk,data.frame(k=k_min,Q=sum(kres$withinss/N)))
    kres_min = kres; rk_all[[k_min]] = kres$cluster

    if ((sum(kres_min$withinss/N) > cval*measures$discrepency) &
        (sum(kres_max$withinss/N) < cval*measures$discrepency)) {
      for (i in k_min:k_max) {
        k_mid = floor( (k_max + k_min)/2)

        if ( k_mid == k_min) {
          break;
        } else {
          kres = kmeansW(measures$M,centers = k_mid,weight=measures$W,nstart = nstart,iter.max = iter.max)
          flog.info("k=%i WSS=%f nstart=%i nfrims=%i",k,sum(kres$withinss),nstart,N)
          rk=rbind(rk,data.frame(k=k_mid,Q=sum(kres$withinss/N)))
          rk_all[[k_mid]] = kres$cluster

          # check which interval to pick
          if (sum(kres$withinss/N) < cval*measures$discrepency) {
            k_max = k_mid
            kres_max = kres
          } else{
            k_min = k_mid
            kres_min = kres
          }
        }
      } # endfor
    }
  }

  rk$cval    = cval
  rk$tsize   = measures$tsize
  rk$vh      = measures$discrepency
  rk$Nw      = ncol(measures$M)
  rk$measure = measures$measure

  # extract first to below discrepency
  if (any(rk$Q<cval*measures$discrepency)) {
    best_k = min(subset(rk,Q<cval*measures$discrepency)$k)
  } else {
    best_k = max(ksupp)
  }
  return(list(summary=rk,all=rk_all,best_k=best_k,best_cluster=rk_all[[best_k]],discrepency=measures$discrepency,measures=measures))
}

#' clusters firms based on their cross-sectional wage distributions
#'
#' @param sdata cross sectional data, needs a column j (firm id) and w (log wage)
#' @param Nw number of points to use for wage distributionsdsd
#' @param k  number of groups
#' @param nstart (default:1000) total number of starting values
#' @param iter.max (default:100) max nunmber of step for each repetition
#' @param measures object created using grouping.getMeasures
#' @param step step size in the repeating
#' @export
grouping.classify.once <- function(measures,k=10,nstart=1000,iter.max=200,step=20) {

  N = measures$N
  flog.info("clustering T=%f, Nw=%i , measure=%s",measures$tsize,measures$Nw,measures$measure);
  kres = kmeansW.repeat(measures$M,centers = k,weight=measures$W,nstart = nstart,iter.max = iter.max,step=step)$best
  rk = data.frame(k=k,Q=sum(kres$withinss/N))
  flog.info("k=%i WSS=%f nstart=%i nfrims=%i",k,sum(kres$withinss),nstart,N)

  rk$cval    = 0
  rk$tsize   = measures$tsize
  rk$vh      = measures$discrepency
  rk$Nw      = ncol(measures$M)
  rk$measure = measures$measure

  return(list(summary=rk,cluster=kres$cluster,best_cluster=kres$cluster,kres=kres,measures=measures))
}

#' Compute the objective function of the clustering
grouping.computeobj <- function(grps) {

  # we collect the info
  rs = grps$measures$rs

  # we attach the weights
  W        = grps$measures$W
  names(W) = grps$measures$fids
  W = as.list(W)
  rs[, weight:= W[[f1]],f1]

  # attach the groups
  G        = grps$cluster
  names(G) = grps$measures$fids
  G = as.list(G)
  rs[, j1 := G[[f1]],f1]

  # compute centroids
  rs[, center := wt.mean(cdf,weight), list(j1,s)]

  # compute weighted sum of squares
  rs[, .N*wt.mean((cdf - center  )^2,weight)]
}


#' Append result of a grouping to a data-set
#'
#' @param sim a data set with jdata and sdata
#' @param clus a clustering result from the classify function
#' @param drop if true firms in sim, not in cluster are dropped
#'
#' @export
grouping.append <- function(sim,clus,drop=F,sort=T){

  # remove cluster if any
  sim$jdata[, j1:=0]
  sim$jdata[, j2:=0]
  sim$sdata[, j1:=0]
  sim$sdata[, j2:=0]

  # create a temporary table to link
  dd = data.table(f = names(clus), j = as.integer(clus))

  # linking
  setkey(dd,f)
  setkey(sim$sdata,f1)
  sim$sdata[,j1 := dd[sim$sdata,j]]
  setkey(sim$jdata,f1)
  sim$jdata[,j1 := dd[sim$jdata,j]]
  setkey(sim$jdata,f2)
  sim$jdata[,j2 := dd[sim$jdata,j]]
  sim$sdata[,j2:=j1]

  # dropping if asked
  if (drop==T) {
    nvals = sim$sdata[,sum(is.na(j1*j2))] + sim$jdata[,sum(is.na(j1*j2))]
    if (nvals>0) {
      flog.warn(" %i observations are missing clusters",nvals)
    }
    sim$jdata <- sim$jdata[j1!=0][j2!=0][!is.na(j1)][!is.na(j2)]
    sim$sdata <- sim$sdata[j1!=0][!is.na(j1)][!is.na(j2)]
  }

  # sort if asked
  if (sort==T) {
    sim <- cluster.order(sim)
  }

  return(sim)
}
