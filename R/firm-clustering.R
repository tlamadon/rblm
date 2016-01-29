#' clusters firms based on their cross-sectional wage distributions
#'
#' @param adata cross sectional data
#' @param ncluster number of clusters to form
#' @param nw number of points to use for wage distributions
#' @param plot default is false, plots the results
#' @export
cluster.firms.data <- function(adata,ncluster=7,nw=40,nstart=100,plot=FALSE,merge=TRUE,max.iter=300,step=20) {
  
  # construct a matrix with the ecdf of wages 
  # together with number of observations
  cat("creating distributions for each firm\n")
  wsup = quantile(adata$lw, (0:nw)/(nw),na.rm=TRUE) 
  dd = adata[fid!="",{ 
    d= melt(hist1(lw,wsup),c('y1')); 
    #d = data.frame()
    #d = rbind(d,data.frame(y1=nw+2,value=mean(lw)))
    #d = rbind(d,data.frame(y1=nw+3,value=sd(lw)))
    d$N = .N; d
  }, fid]  
  dd = cast(dd,fid+N~y1)
  ll = length(names(dd))
  AA = data.matrix(dd[,3:ll])
  WW = dd$N
  fnames = dd$fid
  colnames(AA) <- NULL
  
  cat("starting kmean ... \n")
  clus = kmeansW.repeat(AA,ncluster,WW,iter.max=max.iter,nstart=nstart,step=step) 
    
  # create a matrix of centroid for each rwo in AA
  CC = clus$centers[clus$cluster, ]
  DD = rowMeans( (AA - CC)^2 )
  
  # get belongings
  cluster = data.table(fid=fnames,clus = clus$cluster,dist = DD)
  
  if (merge==FALSE) {
    return(cluster)
  }
  
  setkey(cluster,fid)
  setkey(adata,fid)
  adata[,clus := cluster[adata][,clus]]
  
  # order clusters by mean wage
  cat("done, ordering clusters by mean wage \n")
  I = adata[is.finite(clus),mean(lw),list(clus)][,list(i=rank(V1),clus=clus)]
  setkey(I,clus); I = I[,i]  
  adata[,clus := I[clus]]
  
  if (plot) {
    adata.sum = adata[,list(mw=mean(lw,na.rm=T),my=mean(lprod,na.rm=T),
                            size=log(.N),vw=var(lw,na.rm=T),ge=sum(female)/.N,hed=sum(educ==3)/.N),clus]
    adata.sum=melt(adata.sum,id=c("clus","mw"))
    print(ggplot(adata.sum,aes(x=mw,y=value,color=factor(clus))) + geom_point() + theme_bw() + facet_wrap(~variable,scale="free"))
  }
  
  return(adata)
}

#' clusters firms based on their cross-sectional wage distributions
#'
#' @param adata cross sectional data
#' @param ncluster number of clusters to form
#' @param nw number of points to use for wage distributions
#' @param plot default is false, plots the results
#' @export
cluster.firms.data.meanvar <- function(sdata,ncluster=7,nstart=100,plot=FALSE,merge=TRUE,max.iter=300,step=20) {
  
  # construct a matrix with the ecdf of wages 
  # together with number of observations
  cat("creating moments for each firm\n")
  dd = sdata[fid!="", list(m1 = mean(y1) , m2 = sd(y1), N = .N),fid][!is.na(m1*m2)]
  AA = data.matrix(dd[,list(m1,m2)])
  WW = dd$N
  fnames = dd$fid

  cat("starting kmean ... \n")
  clus = blm:::kmeansW.repeat(AA,ncluster,WW,iter.max=max.iter,nstart=nstart,step=step) 
  
  # create a matrix of centroid for each rwo in AA
  CC = clus$centers[clus$cluster, ]
  DD = rowMeans( (AA - CC)^2 )
  
  # get order by mean wage
  I = rank(clus$centers[,1])
  
  # get belongings
  cluster = data.table(fid=fnames,clus = I[clus$cluster],dist = DD)
  return(cluster)
}

cluster.append.jdata <- function(jdata,clus) {
  setkey(clus,fid)
  setkey(jdata,f1)
  jdata[,j1 := clus[jdata,clus]]
  setkey(jdata,f2)
  jdata[,j2 := clus[jdata,clus]]
  return(jdata)
}

cluster.append.sdata <- function(sdata,clus) {
  setkey(clus,fid)
  setkey(sdata,fid)
  sdata[,j := clus[sdata,clus]]
  return(sdata)
}



#' objective function for clustering based on movers. This function clusters
#' first based on the joint distribution of wages in current firm and wage
#' in the other firm (either coming or moving to). Then it clusters using 
#' the distribution of only wages in the current firm. It compares the 
#' two clusters.
#' @export
cluster.firms.data.movers <- function(jdata,ncluster=7,nw=40,nstart=100,plot=FALSE) {
  
  #jdata2 = jdata[, list( f=c(f1,f2), yh=c(y1,y2), yo=c(y2,y1) ), wid]    
  jdata2 = jdata[, list( f=f1, yh=y1, yo=y2)]    
  wsup  = quantile(c(jdata2$yh,jdata$yo), (0:nw)/(nw),na.rm=TRUE)
  
  # select firms with enough movers
  fids = jdata2[,.N,f][N>=8,f]
  jdata2 = jdata2[f %in% fids]
  
  # PART I : using joint y in firm and y in any other firm
  # ------------------------------------------------------
  # construct the joint density of (yh,yo) for each firm  
  setkey(jdata2,f)
  dd = jdata2[, { d = melt( hist2(yh,yo,wsup)  ,c('y1','y2')); d$N=.N; d} , f]  
  dd = cast(dd,f + N ~y1+y2)
  ll = length(names(dd))
  AA = data.matrix(dd[,3:ll])
  WW = dd$N
  fnames = dd$f
  colnames(AA) <- NULL
  
  clus = kmeansW(AA,ncluster,WW,iter.max=200,nstart=nstart) 
  
  # compute distance from centroid!
  
  
  
  
  dd = data.table(fid=fnames,clus = as.numeric(clus$cluster))
  
  return(dd)    
}



#' clusters firms based on their cross-sectional wage distributions
#'
#' @param adata cross sectional data
#' @param ncluster number of clusters to form
#' @param nw number of points to use for wage distributions
#' @param plot default is false, plots the results
#' @export
cluster.firms.data.all <- function(adata,ncluster=7,nw=40,nstart=100,plot=FALSE,check=FALSE) {
  
  wsup = quantile(adata$lw, (0:nw)/(nw),na.rm=TRUE) 
  
  # I want to select the movers only, then I want to be able to 
  # subset the distribution to use only firms with enough movers
  # - distribution of wages in cross-section
  # - distribtuion of wages for movers
  # - number of workers
  # - number of movers
  setkey(adata,wid)
  m.ids = adata[,list( V1 = all(fid[1]==fid) , fid=fid[1], .N),list(wid,aret)][V1==TRUE][N==4]
  adata2 = adata[wid %in% unique(m.ids[,wid])]
  m.ids = m.ids[,list(all(fid[1]==fid)),wid][V1==FALSE,wid]
  adata2$movers = adata2$wid %in% m.ids
  
  mms = adata2[, {
    # get the crossection wage distribution
    hh  = hist(lw,wsup,plot=FALSE)  
    hh2 = .SD[movers==TRUE, hist(lw,wsup,plot=FALSE) ]
    
  }, fid]
  
  
  
  # select firms with at least 10 movers
  # ------------------------------------
  
  
  
  # part I: cluster using full crossection  
  # --------------------------------------
  fids = adata[,.N,]
  adata2 = adata[]
  
  
  # construct a matrix with the ecdf of wages 
  wsup = quantile(adata$lw, (0:nw)/(nw),na.rm=TRUE) 
  dd = adata[fid!="",hist(lw,wsup,plot=FALSE),fid]
  dd = cast(dd,fid~breaks,value="counts")
  AA = data.matrix(dd[,2:42])
  fnames = dd$fid
  colnames(AA) <- NULL
  
  # cluster
  AA2 = t(apply(AA,1,function(x) {x=cumsum(as.numeric(x)); x=x/x[length(x)]; return(x)}))
  clus = kmeans(AA2,ncluster,iter.max=100,nstart=nstart) 
  
  # part II: cluster using full crossection and movers 
  # --------------------------------------------------
  
  # select 
  
  # construct a matrix with the ecdf of wages of stayers
  year1 = adata[,min(aret)]  
  year2 = adata[,min(aret)]  
  movers.ids = adata[,list(all(fid[1]==fid)), wid][V1==FALSE,wid] 
  jdata = 
  
  stayers.ids = unique()
  dd = adata[]
  wsup = quantile(adata$lw, (0:nw)/(nw),na.rm=TRUE) 
  dd = adata[fid!="",hist(lw,wsup,plot=FALSE),fid]
  dd = cast(dd,fid~breaks,value="counts")
  AA = data.matrix(dd[,2:42])
  fnames = dd$fid
  colnames(AA) <- NULL
  
  # cluster
  AA2 = t(apply(AA,1,function(x) {x=cumsum(as.numeric(x)); x=x/x[length(x)]; return(x)}))
  clus = kmeans(AA2,ncluster,iter.max=100,nstart=nstart) 
  
  
  
  
  # append to data
  cluster = data.table(fid=fnames,clus = clus$cluster)
  setkey(cluster,fid)
  setkey(adata,fid)
  adata[,clus := cluster[adata][,clus]]
  
  # reorder clusters by mean wage
  I = adata[is.finite(clus),mean(lw),list(clus)][,list(i=rank(V1),clus=clus)]
  setkey(I,clus); I = I[,i]  
  adata[,clus := I[clus]]
  
  if (plot) {
    adata.sum = adata[,list(mw=mean(lw,na.rm=T),my=mean(lprod,na.rm=T),
                            size=log(.N),vw=var(lw,na.rm=T),ge=sum(female)/.N,hed=sum(educ==3)/.N),clus]
    adata.sum=melt(adata.sum,id=c("clus","mw"))
    print(ggplot(adata.sum,aes(x=mw,y=value,color=factor(clus))) + geom_point() + theme_bw() + facet_wrap(~variable,scale="free"))
  }
  
  if (check)
  
  return(adata)
}



#' clusters firms based on their cross-sectional wage distributions
#'
#' @param adata cross sectional data
#' @param ncluster number of clusters to form
#' @param nw number of points to use for wage distributions
#' @param plot default is false, plots the results
#' @export
cluster.firms.cross <- function(sdata,ncluster=7,nw=40,nstart=100,plot=FALSE) {
  
  # construct a matrix with the ecdf of wages 
  wsup = quantile(sdata$lw, (0:nw)/(nw),na.rm=TRUE) 
  dd = sdata[,hist(lw,wsup,plot=FALSE),f]
  dd = cast(dd,f~breaks,value="counts")
  AA = data.matrix(dd[,2:42])
  fnames = dd$f
  colnames(AA) <- NULL
  
  AA = t(apply(AA,1,function(x) {x=cumsum(as.numeric(x)); x=x/x[length(x)]; return(x)}))
  clus = kmeans(AA,ncluster,iter.max=100,nstart=nstart) 
  
  cluster = data.table(f=fnames,clus = clus$cluster)
  setkey(cluster,f)
  setkey(sdata,f)
  sdata[,clus := cluster[sdata][,clus]]
  
  # order clusters by mean wage
  I = sdata[is.finite(clus),mean(y),list(clus)][,list(i=rank(V1),clus=clus)]
  setkey(I,clus); I = I[,i]  
  sdata[,clus := I[clus]]
  
  print(ggplot(sdata[,list(size=.N),list(j,clus)],aes(x=j,y=clus,size=(size))) + geom_point() + theme_bw() + geom_abline(linetype=2))
  
  return(sdata)
}




#' objective function for clustering based on movers. This function clusters
#' first based on the joint distribution of wages in current firm and wage
#' in the other firm (either coming or moving to). Then it clusters using 
#' the distribution of only wages in the current firm. It compares the 
#' two clusters.
#' @export
cluster.firms.movers <- function(jdata,ncluster=7,nw=40,nstart=100,plot=FALSE) {
  
  jdata2 = jdata[, list( f=c(f1,f2), yh=c(y1,y2), yo=c(y2,y1) ), wid]    
  jdata2 = jdata[, list( f=f1, yh=y1, yo=y2)]    
  wsup  = quantile(c(jdata2$yh,jdata$yo), (0:nw)/(nw),na.rm=TRUE)
  
  # select firms with enough movers
  fids = jdata2[,.N,f][N>=8,f]
  jdata2 = jdata2[f %in% fids]
  
  # PART I : using joint y in firm and y in any other firm
  # ------------------------------------------------------
  # construct the joint density of (yh,yo) for each firm  
  setkey(jdata2,f)
  dd = jdata2[, { d = melt( hist2(yh,yo,wsup)  ,c('y1','y2')); d$N=.N; d} , f]  
  dd = cast(dd,f + N ~y1+y2)
  ll = length(names(dd))
  AA1 = data.matrix(dd[,3:ll])
  WW1 = dd$N
  fnames1 = dd$f
  colnames(AA1) <- NULL
  
  # PART II : using only y in current firm
  # --------------------------------------
  # construct the joint density of (yh,yo) for each firm  
  dd = jdata2[,{ d= melt(hist1(yh,wsup),c('y1')); d$N = .N; d} , f]  
  dd = cast(dd,f+N~y1)
  ll = length(names(dd))
  AA2 = data.matrix(dd[,3:ll])
  WW2 = dd$N
  fnames2 = dd$f
  colnames(AA2) <- NULL
  
  clus1 = kmeansW(AA1,7,WW1,iter.max=200,nstart=10) 
  clus2 = kmeansW(AA2,7,WW2,iter.max=200,nstart=10) 
  
  gg = data.table(c1 = clus$cluster, c2 = clus2$cluster)
  gg.sort = gg[,quantile(c2,0.5),c1]
  gg.sort = gg[,mean(c2),c1]
  setkey(gg.sort,c1)
  new_label = gg.sort[,rank(V1)]
  gg[,c1_bis := new_label[c1]]
  
  gg2 = gg[,.N,list(c1_bis,c2)]
  
  ggplot(gg2,aes(x=c1_bis,y=c2,size=N)) + geom_point()+theme_bw()
 
}

#' plots the distribution of wage differences
#' conditional on mobility
kline.plot <- function(jdata) {
  
  # KLINE PLOT
  # compute wage growth conditional on j1 and j2
  setkey(jdata,j1,j2)
  wsup = quantile(jdata$y1,seq(0.1,0.9,l=nw))
  jdata$y1m = jdata$y1 - mean(jdata$y1)
  jdata$y2m = jdata$y2 - mean(jdata$y2)
  
  dd  = jdata[,{ d= melt(hist1(y1m-y2m,wsup),c('y1')); d$N = .N;d$x=wsup; d} , list(j1,j2)]  
  dd2 = jdata[,{ d= melt(hist1(y2m-y1m,wsup),c('y1')); d$N = .N;d$x=wsup; d} , list(j1,j2)]  
  dd2$j1 = dd$j2; dd2$j2 = dd$j1
  ggplot(dd,aes(x=x,y=value)) + geom_line() + facet_grid(j1~j2) + geom_line(data=dd2,color='red') + theme_bw()  
}

#' extract information
cluster.infos <- function(cdata) {
    
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




