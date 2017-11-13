

#' Prepare the data for BLM from an employer-employee matched data
#'
#' This relies on the data with the following variables:
#'  - wid worker id
#'  - fid firm id
#'  - y is the log-wage
#'  - t is time
#'
#' @export
jdata.prepare <- function(data) {

  # extract the movers
  mids = data[,length(unique(fid)),wid][V1>0,wid]

  # construct the event data (using wages around the move)
  # we keep all moves, meaning all consecutive spells that are at different firms
  jdata = sdata[wid %in% mids]
  setkey(jdata,wid,t)
  jdata[, iobs := 1:.N, wid]
  setkey(jdata,wid,iobs)
  jdata[, fid.p1 := jdata[J(wid,iobs-1),fid]]
  jdata = jdata[fid!=fid.p1]

  return(jdata)
}

#' plots the conditional quantile distribution for each transitions (l,l')
#' @export
plot.trquant <- function(res,jdata,jdata2=NULL) {

  nf = dim(res$wm)[2]
  # for each j1, j2, compute conditional quantiles of y2 | y1
  rr = data.frame()
  for (l1 in 1:nf) for (l2 in 1:nf) for (tau in ((1:4)*0.2)) {
    fit = rq(y2 ~ bs(y1, df = 8), tau = tau,  data = jdata[(j1==l1) & (j2==l2)])
    y1  = quantile(jdata[(j1==l1) & (j2==l2), y1], seq(0.01,0.99,l=100),na.rm=T)
    y2_hat = predict(fit,newdata = data.frame(y1 =y1))
    rr = rbind(rr,data.frame(tau=tau,y1=y1,y2=y2_hat,j1=l1,j2=l2))
  }

  if (is.null(jdata2)) {
    ggplot(rr,aes(x=y1,y=y2,group=factor(tau))) + geom_line() + facet_grid(j1~j2,scales="free") + theme_bw()
  } else {

    rr2 = data.frame()
    for (l1 in 1:nf) for (l2 in 1:nf) for (tau in ((1:4)*0.2)) {
      fit = rq(y2 ~ bs(y1, df = 8), tau = tau,  data = jdata2[(j1==l1) & (j2==l2)])
      y1  = quantile(jdata2[(j1==l1) & (j2==l2), y1], seq(0.01,0.99,l=100),na.rm=T)
      y2_hat = predict(fit,newdata = data.frame(y1 =y1))
      rr2 = rbind(rr2,data.frame(tau=tau,y1=y1,y2=y2_hat,j1=l1,j2=l2))
    }

    ggplot(rr,aes(x=y1,y=y2,group=factor(tau))) + geom_line(color="blue") + geom_line(data=rr2,color="red",linetype=2) + facet_grid(j1~j2,scales="free") + theme_bw()

  }
  #facet_wrap(j1~j2,scales="free")
}

#' we want to look at the effect of
#' movers on value added.
plot.vaeffect <- function() {

  # I want to compute the wage bill of workers moving in
  # > conditional on j2, sum(y1)
  # > conditional on j1, sum(y2)
  dd1 = jdata[,list(W1=sum(exp(y1)),P2=q2[1]*s2[1],j=j2[1],s2=s2[1]),list(f=f2)]
  dd2 = jdata[,list(W2=sum(exp(y2)),P1=q1[1]*s1[1],j=j1[1],s1=s1[1]),list(f=f1)]
  setkey(dd1,f)
  setkey(dd2,f)
  dd = dd1[dd2]
  rr = dd[is.finite(W1), {r = coef(lm(P2-P1 ~ W1 + W2 + I(s2-s1), .SD) ); list(r,names(r))   },j]

  ggplot(rr,aes(x=j,y=V1)) + geom_line() + facet_wrap(~V2,scales="free") + theme_bw() + geom_hline(yintercept=0,linetype=2)

  # remove cluster mean
  jdata2 = jdata
  jdata2[, y1 := y1 - mean(y1), list(j1,j2)]
  jdata2[, y2 := y2 - mean(y2), list(j1,j2)]
  dd1 = jdata2[,list(W1=sum(exp(y1)),P2=q2[1]*s2[1],j=j2[1],s2=s2[1]),list(f=f2)]
  dd2 = jdata2[,list(W2=sum(exp(y2)),P1=q1[1]*s1[1],j=j1[1],s1=s1[1]),list(f=f1)]
  setkey(dd1,f)
  setkey(dd2,f)
  dd = dd1[dd2]
  rr = dd[is.finite(W1), {r = coef(lm(P2-P1 ~ W1 + W2 + I(s2-s1), .SD) ); list(r,names(r))   },j]

  # check the return to education in each cluster

  rr = adata[, {r = coef(lm(lw ~ educ, .SD) ); list(r,names(r))   },j]


  rr = adata[is.finite(clus), {r = coef(lm(l0 ~ 1 + age + I(age^2), .SD) ); list(r,names(r))   },list(clus,female,educ)]

  rr = adata[is.finite(clus),list(value=quantile(lw,seq(0.1,0.9,0.1)), q=seq(0.1,0.9,0.1)),clus]

}



#' computes simple statistics
get.sample.stats <- function(adata,cdata,jdata) {

  tmp <- function(ddd,name) {
    rr = ddd[fid!="",{
    r = list()
    r$nn = .N
    r$ni = length(unique(wid))
    r$nj = length(unique(fid))

    r$t1 = min(aret)
    r$t2 = max(aret)

    r$mw     = mean(lw)
    r$wsd    = var(lw)

    r$sizem    = .SD[,.N,fid][,mean(N)]
    r$sizem2   = .SD[,rep(.N,.N),fid][,mean(V1)]
    r$sizemed  = .SD[,.N,fid][,median(N)+0.0]
    r$sizemed2 = .SD[,rep(.N,.N),fid][,median(V1)+0.0]

    # number of firms with more than 10 workers
    df = .SD[,.N,fid]
    r$firms_bt_5  = df[N>=5,.N]
    r$firms_bt_10 = df[N>=10,.N]
    r$firms_bt_50 = df[N>=50,.N]

    # number of movers
    m.ids = .SD[,list(all(fid[1]==fid),fid),wid][V1==FALSE]
    r$workers_movers = m.ids[,length(unique(wid))]

    # number of firms with 0,at least 1, at least 5, at least 10 movers
    d.firms = m.ids[,.N,fid]
    r$firms_movers_bt_0 = d.firms[N>=0,.N]
    r$firms_movers_bt_3 = d.firms[N>=3,.N]
    r$firms_movers_bt_5 = d.firms[N>=5,.N]
    r$firms_movers_bt_10 = d.firms[N>=10,.N]
    r$firms_movers_bt_20 = d.firms[N>=20,.N]

    r
  }]

  rr$name=name
  return(rr)
  }

  rr1 = tmp(adata,"all")
  rr2 = tmp(cdata,"3years")
  jdata2 = rbind(  jdata[,list(aret=2002,lw=y1,fid=f1 ,wid=wid)  ] , jdata[,list(aret=2004,lw=y2,fid=f2,wid)  ] )
  rr3 = tmp(jdata2,"3years_movers")

  return(rbind(rr1,rr2,rr3))
}

#' provides statistics
dstats <- function(dd) {
  # checks for columns
  ns = names(dd)
  catf <- function(...) cat(sprintf(...))
  coms <- function(x) prettyNum(x,big.mark=",",scientific=F)

  if ( !( "aret" %in% ns ) & ( "t" %in% ns)) {
    dd[,aret:=t]
  }

  # basic stats
  catf("N=%s ni=%s nj=%s year1=%i year2=%i \n", coms(dd[,.N]), coms(dd[,length(unique(wid))]),coms(dd[,length(unique(fid))]),dd[,min(aret)],dd[,max(aret)])

  if ( !( "lw" %in% ns ) & ( "y" %in% ns)) {
    dd[,lw:=y]
  }

  # compute wage states
  catf("wage  mean=%f sd=%f \n", dd[,mean(lw)], dd[,sd(lw)])

}

rawdata.stats <- function(adata) {

  catf <- function(...) cat(sprintf(...))
  coms <- function(x) prettyNum(x,big.mark=",",scientific=F)

  rr = list()

  # some numbers
  rr$wid_n = rdata[,length(unique(wid))]
  rr$fid_n = rdata[,length(unique(wid))]
  rr$unemployement = rdata[, list(p=mean(employed),.N), list(aret,quarter,educ)]
  rr$employed_with_nafirm = rdata[employed==1, mean(fid=="")]

  # compute transition rates
  rr$transition =  rdata[,list(.N),list(from,educ)][, p := N/sum(N),list(educ)]
  TR = acast(rr$transition, educ ~ from , value.var = 'p')
  colnames(TR) <- c('ee','j2j','u2e', 'u2u' , 'e2u' , 'na')
  rr$TR = TR

  # === look at employment stats
  wid.info = rdata[, list(fcount = length(unique(fid)), ms = sum(monthsworked), employed=sum(employed), nobs = .N) , list(wid,aret)]
  rr$worker_year_obs = wid.info[,.N,nobs]

  wid.info = wid.info[nobs==4]
  rr$employment = wid.info[,.N,list(onefirmonly  = fcount==1, fullyearemployed = ms==12,aret)]

  # select only fully employed in any year
  wids = union(  wid.info[(aret==2002) & (fcount==1) & (ms==12) ,wid] ,  wid.info[(aret==2004) & (fcount==1) & (ms==12) ,wid] , wid.info[(aret==2002) & (fcount==1) & (ms==12) ,wid])
  rdata2 = rdata[wid %in% wids]


  # select only fully employed in 2002 and 2003
  wids = intersect(  wid.info[(aret==2002) & (fcount==1) & (ms==12) ,wid] ,  wid.info[(aret==2004) & (fcount==1) & (ms==12) ,wid] )
  rdata2 = rdata[wid %in% wids]

  rr$transition2 =  rdata2[,list(.N),list(from,educ)][, p := N/sum(N),list(educ)]
  TR = acast(rr$transition2, educ ~ from , value.var = 'p')
  colnames(TR) <- c('ee','j2j','u2e', 'u2u' , 'e2u' , 'na')
  rr$TR2 = TR

  rr$transition2003 =  rdata2[aret==2003,list(.N),list(from,educ)][, p := N/sum(N),list(educ)]
  TR = acast(rr$transition2, educ ~ from , value.var = 'p')
  colnames(TR) <- c('ee','j2j','u2e', 'u2u' , 'e2u' , 'na')
  rr$TR3 = TR


  # aggregate

  # remerge info
  setkey(rdata,wid)
  setkey(wid.info,wid)

  TR = acast(rr$transition, aret+ quarter ~ from , value.var = 'p')
  columns(TR) <- c('ee','j2j','u2e', 'u2u' , 'e2u' , 'na')

  # select 2002 to 2004
  cdata = adata[aret %in% c(2002,2003,2004)]

  catf("number of workers: %i\n", cdata[,length(unique(wid))])
  catf("number of firms: %i\n", cdata[,length(unique(fid))])
  catf("number of workers in at least 2 firms: %i\n", cdata[fcount>1,length(unique(wid))])
  catf("number of workers unemployed for at least 1 quarter: %i\n", cdata[ucount>0,length(unique(wid))])


  # select males and fully employed

  catf("number of workers: %i\n", cdata[,length(unique(wid))])
  catf("number of firms: %i\n", cdata[,length(unique(fid))])
}


