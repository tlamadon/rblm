#' loads the data for given gender 
#' and education. Reformat variables, 
#' and merge in the firm data.
loadData <- function(gender='male',educ=3,cached=TRUE,path="../../qtrdata2012") {
  
  gender = ifelse(gender=='male',0,1)

  if (cached) {
    load(paste(path,'/selectedf',gender,'educ',educ,'.dat',sep=''))
    return(data)
  } 
  
  cat(paste(Sys.time(),"\n"))
  cat('data needs to be reloaded \n') 
  cat(sprintf('loading file %s\n',paste(path,'/selectedf',gender,'educ',educ,'.dta',sep=''))) 
  data  <- data.table(read.dta(paste(path,'/selectedf',gender,'educ',educ,'.dta',sep='')))
  cat('data loaded\n')  
  
  # some renaming
  # --------------------
  data[ , fid := peorgnrdnr2008025]
  data[ , wid := dnr2008025]
  data[ , lw  := logrealmanadslon]
  
  # Merging firm data
  # -----------------
  cat('merging firm information\n')  
  data2  <- data.table(read.dta(paste(path,'/selectedfirms9708.dta',sep='')))
  cat(sprintf("total firms: %i \n",data2[,length(unique(peorgnrdnr2008025))]))
  data2  <- subset(data2,valueadded>0 & ant_anst>0) 
  cat(sprintf("total firms after: %i \n",data2[,length(unique(peorgnrdnr2008025))]))
  data2[,prod  := valueadded/ant_anst]
  data2[,lprod := log(prod)]
  data2 = data2[is.finite(lprod)]
  cat(sprintf("total firms after removing negative VA: %i \n",data2[,length(unique(peorgnrdnr2008025))]))
  fit <- lm(lprod~factor(aret)*industry,data2)
  data2[ , lprodr := residuals(fit)]
  # removing time effects from prod
  data2[ , fid:= peorgnrdnr2008025]
  setkey(data2,fid,aret)
  
  # merging industry and prod data to worker data
  setkey(data,fid,aret)
  data2=data2[,list(fid,aret,prod,lprodr,valueadded,industry,ant_anst)] 
  setkey(data2,fid,aret)
  data = data2[data]
  
  # preparing dummy to show if transition within year
  # =================================================
  # compute dummy that states if worker is changing job within the
  # year
  cat('preparing some worker dummies\n')  
  setkey(data,wid,aret)
  dd.w = data[,list(jchange = length(unique(fid))>1),list(wid,aret)]
  setkey(dd.w,wid,aret)
  data = dd.w[data]
  
  data = data[,list(wid,fid,aret,quarter,birthyear,age,educ,female,from,jobmobility,jchange,kapink,monthsworked,employed,logrealmanadslon,industry,prod,lprodr,valueadded,ant_anst)]

  rm(data2)
  cat('done, saving the file for later use\n')    
  save(data,file=paste(path,'/selectedf',gender,'educ',educ,'.dat',sep=''))
  
  cat(paste(Sys.time(),"\n"))
  
  return(data)
  #save(data,file='qtrdata2012/worker_process.dat')
}


#' this function 
loadDataRaw <- function(year_anchor,resid=FALSE,path="L:/Tibo/qtrdata",cache=FALSE) {
  adata = data.table()
  cat("creating data file\n")
  for (ed in 1:3) for (gd in c('male')) {
    data = loadData(gd,ed,cache,path=path)                  
    data$educ   = ed
    data$female = gd
    data = data[ aret %in% c(year_anchor-2,year_anchor,year_anchor-1,year_anchor+1,year_anchor+2)]   

    # merging 
    adata = rbind(adata,data)
    cat(sprintf("done with %s-%i\n",gd,ed))
    flush.console()
  }

  cat("done creating data file, saving for future use\n")
  flush.console()
  return(adata)
}


# loads the data for males, removes time effect, selects only full employment years
loadDataGapYearsWithSpells <- function(year_anchor,resid=FALSE,path="L:\\Tibo\\qtrdata") {
  
  filename = sprintf("%s\\all-spells-anchor%i.dat",path,year_anchor)
  adata = data.table()
  
  cat(sprintf("preparing file %s \n",filename))
  flush.console()
  
  if (file.exists(filename)) {
    cat(sprintf("loading %s from disk \n",filename))
    load(file=filename)
  } else {
    cat("creating data file\n")
    for (ed in 1:3) for (gd in c('male')) {
      data = loadData(gd,ed,FALSE)
            
      # find type of transitions for each worker between year1 and year5
      dtrans = data[,  list( ucount = sum(fid==""), fcount = length(unique(fid))   ), wid  ]
      
      # make the data yearly, keep only full time employements
      setkey(data,wid,aret)
      data = data[, {
        rs = list(birthyear=birthyear[1],fid = fid[1],prod = prod[1], ant_anst = ant_anst[1], industry = industry[1],valueadded = valueadded[1])
        rs$ms     = sum(monthsworked)
        rs$fcount = length(unique(fid))
        rs$lw     = mean(logrealmanadslon,na.rm=T)
        rs
      }, list(wid,aret)]

      # append info about next and past income and month worked
      setkey(data,wid,aret)
      data[,lw.l1:=data[J(wid,aret-1),lw]]
      data[,ms.l1:=data[J(wid,aret-1),ms]]
      data[,lw.f1:=data[J(wid,aret+1),lw]]
      data[,ms.f1:=data[J(wid,aret+1),ms]]
      # EM: We select people with 12 months of employment, one firm id and firm id not blank.
      data = data[ms==12][fcount==1][fid!=""]
      
      # extract relevant spells
      # create lag fid
      setkey(data,wid,aret)
      data[, fid.l1 := data[J(wid,aret-1),fid][,fid]]
      data[, fid.l1 := ifelse(is.na(fid.l1),"na",fid.l1)]
      data[, spell  := cumsum( fid != fid.l1), wid]
            
      # merge in the move info
      setkey(dtrans,wid)
      setkey(data,wid)
      data[,fcount:=NULL]
      data[,fid.l1:=NULL]
      data = dtrans[data]
      data$educ   = ed
      data$female = gd

      # remove time effect
      data[, lw0:=lw ]
      rr = data[aret==year_anchor ,mean(lw),birthyear]
      data[, lw:= lw - mean(lw), list(birthyear,aret)]
      setkey(rr,birthyear)
      setkey(data,birthyear)
      data = rr[data]
      data[,lw := lw + V1][,V1:=NULL]
      data = data[!is.na(lw)]

      # merge the data
      adata = rbind(adata,data)
      cat(sprintf("done with %s-%i\n",gd,ed))
      flush.console()
    }
    cat("done creating data file, saving for future use\n")
    flush.console()
    save(adata,file=filename)
  }
  
  cat("done!\n") 
  return(adata)
}




data.summary <- function(adata) {
  
  r = list()
  
  r$year1 = adata[,min(aret)]
  r$year2 = adata[,max(aret)]
  
  # total number of workers
  r$total_worker = length(unique(adata$wid))
  r$total_firm   = length(unique(adata$fid))
  
  # number of firms with more than 10 workers
  df = adata[,.N,fid]
  r$firms_bt_10 = df[N>=10,.N]
  r$firms_bt_50 = df[N>=50,.N]
  
  # number of workers present in 2 periods
  setkey(adata,wid)
  r$workers_in2periods = adata[,list( length(unique(aret)) >=2 ), wid][V1==TRUE,.N]
  
  # number of movers
  m.ids = adata[,list( V1 = all(fid[1]==fid) , fid=fid[1], .N),list(wid,aret)][V1==TRUE][N==4]
  r$workers_full_periods = m.ids[,length(unique(wid))]
  dcross = adata[wid %in% m.ids[,wid]][aret==r$year1][,list(lw=mean(lw),fid=fid[1]),wid]
  m.ids = m.ids[,list(all(fid[1]==fid),fid),wid][V1==FALSE]
  r$workers_movers = m.ids[,length(unique(wid))]
  
  # number of firms with 0,at least 1, at least 5, at least 10 movers
  d.firms = m.ids[,.N,fid]
  r$firms_movers_bt_0 = d.firms[N>=0,.N]
  r$firms_movers_bt_3 = d.firms[N>=3,.N]
  r$firms_movers_bt_5 = d.firms[N>=5,.N]
  r$firms_movers_bt_10 = d.firms[N>=10,.N]
  r$firms_movers_bt_20 = d.firms[N>=20,.N]
  
  # wages
  r$wage_var = dcross[fid!=""][is.finite(lw),sd(lw)]
  r$wage_var_within = dcross[fid!=""][is.finite(lw)][,lw - mean(lw),fid][,sd(V1)]
  r$wage_var_between = dcross[fid!=""][is.finite(lw)][,mean(lw),fid][,sd(V1)]
  
  # for (n in names(r)) { cat(sprintf("%20s: %f \n", n , r[[n]]))}
  return(r)
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


#' extract spell data from the full data set.
#' We want to select a very short sample
get.shortPanelFromSpells <- function(adata,year1=2002,year2=2004) {

  # append past and future wages
  # append start and end date of the spell
  dd2 = adata[female=="male"]
  
  dd2[, t1:= min(aret), list(wid,spell)]
  dd2[, t2:= max(aret), list(wid,spell)]
  dd2[, lwb:= mean(lw,na.rm=T),list(wid,spell)]
  setkey(dd2,wid,spell,aret)
  dd2[, lw.l1  := dd2[J(wid,spell,aret-1),list(lw)][,lw]]
  dd2[, lw.f1  := dd2[J(wid,spell,aret+1),list(lw)][,lw]]
  setkey(dd2,wid,aret)
  dd2[, lw.l1uc := dd2[J(wid,aret-1),list(lw)][,lw]]
  dd2[, lw.f1uc := dd2[J(wid,aret+1),list(lw)][,lw]]
  dd2[, lw.l2uc := dd2[J(wid,aret-2),list(lw)][,lw]]
  dd2[, lw.f2uc := dd2[J(wid,aret+2),list(lw)][,lw]]
  dd2 = dd2[ (aret>=year1) & (aret <= year2)]
    
  # for each spell, we compute the mean earnings
  # setkey(adata,wid,spell)
  # cdata = adata[,list(lw=mean(lw),fid=fid[1],prod=prod[1],birthyear=birthyear[1],female=female[1],educ=educ[1],
  #                    t1=min(aret),t2=max(aret),len=.N,ucount=ucount[1]),list(wid,spell)]
  # compute demographics
  
  dd2[ , ageg := cut(birthyear,c(0,1957,1968,2010))]
  dd2[ , xg := interaction(educ,female,ageg)]
  dd2[ , x := as.integer(xg)]
  
  return(dd2)
}

#' we extract movers wages in year1 and year2
get.MoversShortPanel <- function(cdata,continuingFirms=T) {
  
  year1 = min(cdata$aret)
  year2 = max(cdata$aret)
  
  # keep firms present in the 2 years
  if (continuingFirms) {
    fids = unique(cdata[, (year1 %in% aret) & (year2 %in% aret), fid][V1==TRUE,fid])
    dd2 = cdata[fid %in% fids]
  } else {
    dd2 = cdata
  }
  
  # select workers with 2 spells
  dd2   = dd2[(aret==year1) | (aret==year2)]
  wids   = dd2[, .N , wid][N>=2,wid]
  jdata  = dd2[wid %in% wids]
  
  # select best 2 spells  
  setkey(jdata,wid,aret)
  jdata = jdata[, list(j1 = clus[1],j2=clus[2],y1=lw[1],y2=lw[2],yb1=lwb[1],yb2=lwb[2],
                       t11=t1[1],t12=t2[1],t21=t1[2],t22=t2[2],
                       y1.l1=lw.l1[1],y2.f1=lw.f1[2],spell1 = spell[1],spell2=spell[2],
                       f1=fid[1],f2=fid[2],x=x[1],ucount=ucount[1],s1=ant_anst[1],s2=ant_anst[2]),wid]
  
  jdata = jdata[f1!=f2]
    
  return(jdata)
}

#' we extract movers wages in year1 and year2
get.MoversShortPanelEnd <- function(cdata) {
  
  year1 = min(cdata$aret)
  year2 = max(cdata$aret)
  
  # keep firms present in the 2 years
  fids = unique(cdata[, (year1 %in% aret) & (year2 %in% aret), fid][V1==TRUE,fid])
  dd2 = cdata[fid %in% fids]
  
  # select workers with 2 spells
  dd2   = dd2[(aret==year1) | (aret==year2)]
  wids   = dd2[, .N , wid][N>=2,wid]
  jdata  = dd2[wid %in% wids]
  
  # select best 2 spells  
  setkey(jdata,wid,aret)
  jdata = jdata[, list(j1 = clus[1],j2=clus[2],y2=lw[1],y3=lw[2],
                       t11=t1[1],t12=t2[1],t21=t1[2],t22=t2[2],
                       y1=lw.l1uc[1],y4=lw.f1uc[2],
                       y0=lw.l2uc[1],y5=lw.f2uc[2],
                       y25 = lw.f1[1],
                       spell1 = spell[1],spell2=spell[2],
                       f1=fid[1],f2=fid[2],x=x[1],ucount=ucount[1],s1=ant_anst[1],s2=ant_anst[2]),wid]
  
  jdata = jdata[f1!=f2][!is.na(y1*y2*y3*y4)]
  
  return(jdata)
}

#' we extract movers wages in year1 and year2
get.StayersShortPanelEnd2 <- function(cdata) {
  
  year1 = min(cdata$aret)
  year2 = max(cdata$aret)
  
  # keep firms present in the 2 years
  fids = unique(cdata[, (year1 %in% aret) & (year2 %in% aret), fid][V1==TRUE,fid])
  dd2  = cdata[fid %in% fids]
  
  # select workers only in 1 firm
  dd2   = dd2[(aret==year1) | (aret==year2)]
  wids   = dd2[, .N , wid][N>=2,wid]
  sdatae2  = dd2[wid %in% wids]
  
  # select best 2 spells  
  setkey(sdatae2,wid,aret)
  sdatae2 = sdatae2[, list(j1 = clus[1],y2=lw[1],y3=lw[2],
                       t11=t1[1],t12=t2[1],t21=t1[2],t22=t2[2],
                       y1=lw.l1uc[1],y4=lw.f1uc[2],
                       y25 = lw.f1[1],
                       spell1 = spell[1],spell2=spell[2],
                       y0=lw.l2uc[1],y5=lw.f2uc[2],
                       f1=fid[1],f2=fid[2],x=x[1],ucount=ucount[1],s1=ant_anst[1],s2=ant_anst[2]),wid]
  
  sdatae2 = sdatae2[f1==f2][!is.na(y1*y2*y3*y4)]
  
  return(sdatae2)
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


