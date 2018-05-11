# This is an em for the exognous mobility case.
# It will estimate a non-stationary model, and will be able to impose monotonicity

# ------------- Initiliazing functions ---------------------


#' create a random model for EM with
#' endogenous mobility with multinomial pr
#' @export
m2.mixt.new <-function(nk,nf,fixb=F,stationary=F) {

  model = list()
  # model for Y1|Y2,l,k for movers and stayes
  model$A1    = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))
  model$S1    = 0.3*array(1+0.5*runif(nf*nk),c(nf,nk))
  # model for Y4|Y3,l,k for movers and stayes
  model$A2    = array(0.9*(1 + 0.5*rnorm(nf*nk)),c(nf,nk))
  model$S2    = 0.3*array(1+0.5*runif(nf*nk),c(nf,nk))
  # model for p(K | l ,l') for movers
  model$pk1    = rdirichlet(nf*nf,rep(1,nk))
  # model for p(K | l ,l') for stayers
  model$pk0    = rdirichlet(nf,rep(1,nk))
  dim(model$pk0) = c(1,dim(model$pk0))

  model$nk    = nk
  model$nf    = nf

  for (l in 1:nf) {
    model$A1[l,] = sort(model$A1[l,])
    model$A2[l,] = sort(model$A2[l,])
  }

  if (fixb) {
    model$A2 = spread(rowMeans(model$A2),2,nk) + model$A1 - spread(rowMeans(model$A1),2,nk)
  }

  if (stationary) {
    model$A2 = model$A1
  }

  return(model)
}





# ------------- Simulating functions ---------------------

#' Using the model, simulates a dataset of movers
#' @export
m2.mixt.simulate.movers <- function(model,NNm=NA) {

  J1 = array(0,sum(NNm))
  J2 = array(0,sum(NNm))
  Y1 = array(0,sum(NNm))
  Y2 = array(0,sum(NNm))
  K  = array(0,sum(NNm))

  A1  = model$A1
  A2  = model$A2
  S1  = model$S1
  S2  = model$S2
  pk1  = model$pk1
  nk  = model$nk
  nf  = model$nf

  i =1
  for (l1 in 1:nf) for (l2 in 1:nf) {
    I = i:(i+NNm[l1,l2]-1)
    ni = length(I)
    jj = l1 + nf*(l2 -1)
    J1[I] = l1
    J2[I] = l2

    # draw k
    Ki = sample.int(nk,ni,T,pk1[jj,])
    K[I] = Ki

    # draw Y2, Y3
    Y1[I]  = A1[l1,Ki] + S1[l1,Ki] * rnorm(ni)
    Y2[I]  = A2[l2,Ki] + S2[l2,Ki] * rnorm(ni)

    i = i + NNm[l1,l2]
  }

  jdatae = data.table(k=K,y1=Y1,y2=Y2,j1=J1,j2=J2)
  return(jdatae)
}

#' Using the model, simulates a dataset of stayers.
#' @export
m2.mixt.simulate.stayers <- function(model,NNs) {

  J1 = array(0,sum(NNs))
  J2 = array(0,sum(NNs))
  Y1 = array(0,sum(NNs))
  Y2 = array(0,sum(NNs))
  K  = array(0,sum(NNs))

  A1  = model$A1
  A2  = model$A2
  S1  = model$S1
  S2  = model$S2
  pk0 = model$pk0
  nk  = model$nk
  nf  = model$nf

  # ------  impute K, Y1, Y4 on jdata ------- #
  i =1
  for (l1 in 1:nf) {
    I = i:(i+NNs[l1]-1)
    ni = length(I)
    J1[I] = l1

    # draw k
    Ki = sample.int(nk,ni,T,pk0[1,l1,])
    K[I] = Ki

    # draw Y2, Y3
    Y1[I]  = A1[l1,Ki] + S1[l1,Ki] * rnorm(ni)
    Y2[I]  = A2[l1,Ki] + S2[l1,Ki] * rnorm(ni)

    i = i + NNs[l1]
  }

  sdatae = data.table(k=K,y1=Y1,y2=Y2,j1=J1,j2=J1,x=1)

  return(sdatae)
}

#' Using the model, simulates a dataset of stayers.
#' @export
m2.mixt.simulate.stayers.withx <- function(model,NNsx) {

  J1 = array(0,sum(NNsx))
  J2 = array(0,sum(NNsx))
  Y1 = array(0,sum(NNsx))
  Y2 = array(0,sum(NNsx))
  K  = array(0,sum(NNsx))
  X  = array(0,sum(NNsx))

  A1  = model$A1
  A2  = model$A2
  S1  = model$S1
  S2  = model$S2
  pk0 = model$pk0
  nk  = model$nk
  nf  = model$nf
  nx  = nrow(NNsx)

  # ------  impute K, Y1, Y4 on jdata ------- #
  i =1
  for (l1 in 1:nf) for (x in 1:nx) {
    I = i:(i+NNsx[x,l1]-1)
    ni = length(I)
    J1[I] = l1

    # draw k
    Ki = sample.int(nk,ni,T,pk0[x,l1,])
    K[I] = Ki
    X[I] = x

    # draw Y2, Y3
    Y1[I]  = A1[l1,Ki] + S1[l1,Ki] * rnorm(ni)
    Y2[I]  = A2[l1,Ki] + S2[l1,Ki] * rnorm(ni)

    i = i + NNsx[x,l1]
  }

  sdatae = data.table(k=K,y1=Y1,y2=Y2,j1=J1,j2=J1,x=X)
  return(sdatae)
}


#' @export
m2.mixt.impute.movers <- function(jdatae,model) {

  A1  = model$A1
  S1  = model$S1
  pk1 = model$pk1
  A2  = model$A2
  S2  = model$S2
  nk  = model$nk
  nf  = model$nf

  # ------  impute K, Y1, Y4 on jdata ------- #
  jdatae.sim = copy(jdatae)
  jdatae.sim[, c('k_imp','y1_imp','y2_imp') := {
    ni = .N
    jj = j1 + nf*(j2-1)
    Ki  = sample.int(nk,.N,prob = pk1[jj,],replace=T)
    # draw Y1, Y4
    Y1 = rnorm(ni)*S1[j1,Ki] + A1[j1,Ki]
    Y2 = rnorm(ni)*S2[j2,Ki] + A2[j2,Ki]
    list(Ki,Y1,Y2)
  },list(j1,j2)]

  return(jdatae.sim)
}
#' @export
m2.mixt.impute.stayers <- function(sdatae,model) {

  A1  = model$A1
  S1  = model$S1
  pk0 = model$pk0
  A2  = model$A2
  S2  = model$S2
  nk  = model$nk
  nf  = model$nf

  # ------  impute K, Y1, Y4 on jdata ------- #
  sdatae.sim = copy(sdatae)
  sdatae.sim[, c('k_imp','y1_imp','y2_imp') := {
    ni = .N
    Ki  = sample.int(nk,.N,prob = pk0[x,j1,],replace=T)
    # draw Y2, Y3
    Y1  = A1[j1,Ki] + S1[j1,Ki] * rnorm(ni)
    Y2  = A2[j1,Ki] + S2[j1,Ki] * rnorm(ni) # false for movers
    list(Ki,Y1,Y2)
  },list(j1,x)]

  return(sdatae.sim)
}

#' Simulates data (movers and stayers) and attached firms ids. Firms have all same expected size.
#' @export
m2.mixt.simulate.sim <- function(model,fsize,smult=1,mmult=1) {
  jdata       = m2.mixt.simulate.movers(model,model$NNm*mmult)
  sdata       = m2.mixt.simulate.stayers(model,model$NNs*smult)

  # create some firm ids
  sdata <- sdata[,f1 := paste("F",j1 + model$nf*(sample.int(.N/fsize,.N,replace=T)-1),sep=""),j1]
  sdata <- sdata[,j1b:=j1]
  sdata <- sdata[,j1true := j1]
  jdata <- jdata[,j1true := j1][,j2true := j2]
  jdata <- jdata[,j1c:=j1]
  jdata <- jdata[,f1:=sample( unique(sdata[j1b %in% j1c,f1]) ,.N,replace=T),j1c]
  jdata <- jdata[,j2c:=j2]
  jdata <- jdata[,f2:=sample( unique(sdata[j1b %in% j2c,f1])  ,.N,replace=T),j2c]
  jdata$j2c=NULL
  jdata$j1c=NULL
  sdata$j1b=NULL
  sdata[,f2:=f1]

  sim = list(sdata=sdata,jdata=jdata)
  return(sim)
}


# -------------------- Estimating functions -----------------------------


#' Estimates the static model parameters for movers
#'
#' @export
m2.mixt.movers <- function(jdatae,model,ctrl) {

  start.time <- Sys.time()
  tic <- tic.new()

  dprior = ctrl$dprior
  model0 = ctrl$model0
  taum   = ctrl$tau

  ### ----- GET MODEL  ---
  nk  = model$nk
  nf  = model$nf
  A1  = model$A1
  S1  = model$S1
  A2  = model$A2
  S2  = model$S2
  pk1 = model$pk1

  # ----- GET DATA
  # movers
  Y1m = jdatae$y1
  Y2m = jdatae$y2
  J1m = jdatae$j1
  J2m = jdatae$j2
  JJm = J1m + nf*(J2m-1)
  Nm = jdatae[,.N]

  # get the constraints
  CS1  = cons.pad(cons.get(ctrl$cstr_type[1],ctrl$cstr_val[1],nk,nf),nk*nf*0, nk*nf*1)
  CS2  = cons.pad(cons.get(ctrl$cstr_type[2],ctrl$cstr_val[2],nk,nf),nk*nf*1,0)
  # combine them
  CS = cons.bind(CS1,CS2)

  # create the stationary contraints
  if (ctrl$fixb==T) {
    CS2 = cons.fixb(nk,nf,2)
    CS  = cons.bind(CS2,CS)
  }

  # create a constraint for the variances
  if (ctrl$model_var==T) {
    CSw = cons.none(nk,nf*2)
  } else{
    CS1  = cons.pad(cons.mono_k(nk,nf),nk*nf*0, nk*nf*3)
    CS2  = cons.pad(cons.mono_k(nk,nf),nk*nf*1, nk*nf*2)
    CSw  = cons.bind(CS1,CS2)
    CSw$meq = length(CSw$H)
  }

  # prepare matrices aggregated at the type level
  Dkj1f     = diag(nf)  %x% rep(1,nf) %x% diag(nk)            # A[k,l] coefficients for j1
  Dkj2f     = rep(1,nf) %x% diag(nf)  %x% diag(nk)            # A[k,l] coefficients for j2

  # regression matrix for the variance
  XX = rBind(
    cBind(    Dkj1f, 0*Dkj2f),
    cBind(  0*Dkj1f,   Dkj2f)
  )

  ## --- prepare regressions covariates --- #
  # create the depend variables

  lik_old  = -Inf
  lik      = -Inf
  lik_best = -Inf
  liks = 0
  likm=0

  lpt1 = array(0,c(Nm,nk))
  lpt2 = array(0,c(Nm,nk))
  lp   = array(0,c(Nm,nk))

  tic("prep")
  stop = F;
  for (step in 1:ctrl$maxiter) {

    model1 = list(nk=nk,nf=nk,A1=A1,A2=A2,S1=S1,S2=S2,
                  pk1=pk1,dprior=dprior)

    ### ---------- E STEP ------------- #
    # compute the tau probabilities and the likelihood
    if (is.na(taum[1]) | (step>1)) {

      # for efficiency we want to group by (l1,l2)
      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        ll = l1 + nf*(l2-1)
        if (length(I)==0) next;

        for (k in 1:nk) {
          lpt1[I,k] = lognormpdf(Y1m[I] , A1[l1,k], S1[l1,k])
          lpt2[I,k] = lognormpdf(Y2m[I] , A2[l2,k], S2[l2,k])

          # sum the log of the periods
          lp[I,k] = log(pk1[ll,k]) + lpt1[I,k] + lpt2[I,k]
        }
      }

      liks     = sum(logRowSumExp(lp))
      taum     = exp(lp - spread(logRowSumExp(lp),2,nk)) # normalize the k probabilities Pr(k|Y1,Y2,Y3,Y4,l)

      # compute prior
      lik_prior = (dprior-1) * sum(log(pk1))
      lik = liks + lik_prior

    } else {
      cat("skiping first max step, using supplied posterior probabilities\n")
    }

    tic("estep")

    if (stop) break;

    # ---------- MAX STEP ------------- #
    # taum = makePosteriorStochastic(tau = taum,m = ctrl$stochastic) # if we want to implement stochastic EM

    # we start by recovering the posterior weight, and the variances for each term
    rwm  = c(t(taum + ctrl$posterior_reg))

    if (ctrl$fixm==F) {
      DYY  = array(0,c(nk,nf,nf,2))
      WWT  = array(1e-7,c(nk,nf,nf,2))

      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        if (length(I)==0) next;
        for (k in 1:nk) {
          # compute the posterior weight, it's not time specific
          ww = sum(taum[I,k] + ctrl$posterior_reg)

          # construct dependent for each time period k,l2,l1,
          DYY[k,l2,l1,1] = sum(  Y1m[I] * (taum[I,k] + ctrl$posterior_reg) )/ww
          DYY[k,l2,l1,2] = sum(  Y2m[I] * (taum[I,k] + ctrl$posterior_reg) )/ww

          # Scaling the weight by the time specific variance
          WWT[k,l2,l1,1] = ww/pmax(ctrl$sd_floor,S1[l1,k]^2)
          WWT[k,l2,l1,2] = ww/pmax(ctrl$sd_floor,S2[l2,k]^2)
        }
      }

      WWT = WWT/sum(WWT)
      fit = slm.wfitc(XX,as.numeric(DYY),as.numeric(WWT),CS)$solution
      is  = 1
      A1  = t(rdim(fit[is:(is + nk*nf-1)],nk,nf)); is = is+nk*nf
      A2  = t(rdim(fit[is:(is + nk*nf-1)],nk,nf)); is = is+nk*nf

      # compute the variances!!!!
      DYY_bar   = array(0,c(nk,nf,nf,2))
      DYY_bar[] = XX%*%fit
      DYYV      = array(0,c(nk,nf,nf,2))

      for (l1 in 1:nf) for (l2 in 1:nf) {
        I = which( (J1m==l1) & (J2m==l2))
        if (length(I)==0) next;
        for (k in 1:nk) {
          # construct dependent for each time period k,l2,l1,
          ww = sum(taum[I,k] + ctrl$posterior_reg)
          DYYV[k,l2,l1,1] = sum(  (Y1m[I]  - DYY_bar[k,l2,l1,1])^2 * (taum[I,k] + ctrl$posterior_reg) )/ww
          DYYV[k,l2,l1,2] = sum(  (Y2m[I]  - DYY_bar[k,l2,l1,2])^2 * (taum[I,k] + ctrl$posterior_reg) )/ww
        }
      }

      fitv  = slm.wfitc(XX,as.numeric(DYYV),as.numeric(WWT),CSw)$solution
      is    = 1
      S1    = sqrt(t(rdim(fitv[is:(is + nk*nf-1)],nk,nf))); is = is+nk*nf
      S2    = sqrt(t(rdim(fitv[is:(is + nk*nf-1)],nk,nf))); is = is+nk*nf
      S1[S1<ctrl$sd_floor]=ctrl$sd_floor # having a variance of exacvtly 0 creates problem in the likelihood
      S2[S2<ctrl$sd_floor]=ctrl$sd_floor
    }
    tic("mstep-ols")

    ## -------- PK probabilities ------------ #
    ## --- movers --- #
    for (l1 in 1:nf) for (l2 in 1:nf) {
      jj = l1 + nf*(l2-1)
      I = which(JJm == jj)
      if (length(I)>1) {
        pk1[jj,] = colSums(taum[I,])
      } else if (length(I)==0) { # this deals with the case where the cell is empty
        pk1[jj,] = 1/nk
      } else {
        pk1[jj,] = taum[I,]
      }
      pk1[jj,] = (pk1[jj,] + dprior-1 )/(sum(pk1[jj,] + dprior -1 ))
    }

    #check_lik = computeLik(Y1m,Y2m,Y3m,Y4m,A12,B12,S12,A43,B43,S43,A2ma,A2mb,S2m,A3ma,A3mb,B32m,S3m)
    #if (check_lik<lik) cat("lik did not go down on pk1 update\n")

    tic("mstep-pks")

    # checking model fit
    if ((!any(is.na(model0))) & ((step %% ctrl$nplot) == (ctrl$nplot-1))) {
      I1 = order(colSums(A1))
      I2 = order(colSums(model0$A1))
      rr = addmom(A2[,I1],model0$A2[,I2],"A2")
      rr = addmom(A1[,I1],model0$A1[,I2],"A1",rr)
      rr = addmom(S2[,I1], model0$S2[,I2], "S2", rr,type="var")
      rr = addmom(S1[,I1], model0$S1[,I2], "S1", rr,type="var")
      rr = addmom(pk1,model0$pk1,"pk1",rr,type="pr")

      print(ggplot(rr,aes(x=val2,y=val1,color=type)) + geom_point() + facet_wrap(~name,scale="free") + theme_bw() + geom_abline(linetype=2))
    } else {
      if ((step %% ctrl$nplot) == (ctrl$nplot-1)) {
        wplot(A1)
      }
    }

    # -------- check convergence ------- #
    dlik = (lik - lik_old)/abs(lik_old)
    lik_old = lik
    lik_best = pmax(lik_best,lik)
    if ( (step %% ctrl$ncat) == 0) flog.info("[%3i][%s] lik=%4.4f dlik=%4.4e liks=%4.4e likm=%4.4e",step,ctrl$textapp,lik,dlik,liks,likm);
    if (step>10) if (abs(dlik)<ctrl$tol) break;

    tic("loop-wrap")
  }
  flog.info("[%3i][%s][final] lik=%4.4f dlik=%4.4e liks=%4.4e likm=%4.4e",step,ctrl$textapp,lik,dlik,liks,likm);

  # Y1 | Y2
  model$A1  = A1
  model$S1  = S1
  model$A2  = A2
  model$S2  = S2
  ## --movers --
  model$pk1  = pk1

  model$NNm  = acast(jdatae[,.N,list(j1,j2)],j1~j2,fill=0,value.var="N")
  model$likm = lik

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  return(list(tic = tic(), model=model,lik=lik,step=step,dlik=dlik,time.taken=time.taken,ctrl=ctrl,liks=liks,likm=likm))
}

m2.mixt.rdim.pk1 <-function(pk1) {





}



#' use the marginal distributions to extract type distributions
#' within each cluster and observable characteristics
#' @export
m2.mixt.stayers <- function(sdata,model,ctrl) {

  # we set a linear programing problem to maximize likelihood subject
  # to non negetivity and summing to one

  # the objective weights are the the density evaluated at each k
  nk  = model$nk
  nf  = model$nf
  Y1  = sdata$y1   # firm id in period 1
  J1  = sdata$j1   #  wage in period 1
  X   = sdata$x    # observable category
  # @todo add code in case X is missing, just set it to one
  nx  = length(unique(X))
  N   = length(Y1)
  Wmu = t(model$A1)
  Wsg = t(model$S1)

  # we create the index for the movement
  # this needs to take into account the observable X
  J1x = X + nx*(J1-1) # joint in index for movement
  J1s <- Matrix(0, nrow = N, ncol = nf * nx, sparse = TRUE)
  II = 1:N + N*( J1x -1 ); J1s[II]=1
  tot_count = t(spread(Matrix::colSums(J1s),2,nk))
  empty_cells = (tot_count[1,]==0)

  #PI = rdirichlet(nf*nx,rep(1,nk))
  PI = rdim(model$pk0,nf*nx,nk)
  PI_old = PI

  lik_old = Inf
  iter_start =1

  for (count in iter_start:ctrl$maxiter) {
    # the coeffs on the pis are the sum of the norm pdf
    norm1 = dnorm(spread(Y1,2,nk),t(Wmu[,J1]),t(Wsg[,J1]))
    tau   = PI[J1x,]*norm1
    tsum  = Matrix::rowSums(tau)
    tau   = tau / spread( tsum ,2,nk  )
    lik   = - sum(log(tsum))

    PI    = t.default( as.array( t(tau) %*% J1s / tot_count  ))
    PI[empty_cells,] = array(1/nk,c(sum(empty_cells),nk))

    dPI = abs(PI - PI_old)
    max_change  = max(dPI)
    mean_change = mean(dPI)
    PI_old = PI

    if (!is.finite(lik)) { status = -5; break; }

    prg = (lik_old - lik)/lik
    lik_old = lik

    if ((count %% ctrl$ncat)==(ctrl$ncat-1)) {
      flog.info("[%3i][%s] lik=%4.4e inc=%4.4e max-pchg=%4.4e mean-pchg=%4.4e",count,ctrl$textapp,lik,prg,max_change,mean_change)
      flush.console()
    }

    if (max_change<ctrl$tol) {
      status = 1;
      msg = "converged";
      break;
    }

  }

  model$pk0  = rdim(PI,nx,nf,nk)
  model$liks = lik
  model$NNs  = sdata[,.N,j1][order(j1)][,N]

  return(model)
}


#' Estimates the static mixture model on 2 periods
#'
#' This estimator uses multiple starting values to try to find the global maxima.
#'
#' @export
m2.mixt.estimate.all <- function(sim,nk=6,ctrl,cl=NA) {
  start.time <- Sys.time()
  sdata = sim$sdata
  jdata = sim$jdata

  mm = mean(sdata$y1)
  ms = 2*sd(sdata$y1)

  # check that sdata has an x column
  if (!("x" %in% names(sdata))) {
    flog.info("creating an x column in sdata and set it to 1")
    sdata$x=1
  } else if (length(unique(sdata$x)) >= 50 ) {
    stop("likely too many values in the x column of sdata")
  }

  nf = max(sdata$j1);
  model_start = m2.mixt.new(nk,nf)

  res_para = m2.mixt.movers(jdata,model_start,ctrl=em.control(ctrl,cstr_type="para",textapp="para0",fixb=F))

  # use cluster if available
  if (!any(is.na(cl))) {
    flog.info("cluster -- exporting objects to nodes")
    # export environment to nodes
    clusterExport(cl,c("res_para","jdata","ctrl"),environment())
    mylapply <- function(...) parLapply(cl,...)
    nnodes=length(cl)
  } else {
    mylapply <- function(...) lapply(...)
    nnodes=1
  }

  flog.info("starting repetitions with %i nodes",nnodes)
  rr  = mylapply(1:ctrl$est_rep, function(i) {
    res_mixt = list()
    tryCatch({
      res_para$model$A1 = spread(sort(rnorm(nk))*ms+mm,1,nf)
      res_para$model$A2 = res_para$model$A1
      res_para_fixm = m2.mixt.movers(jdata,res_para$model,ctrl=em.control(ctrl,cstr_type="para",textapp=sprintf("paraf (%i/%i)",i,ctrl$est_rep),fixm=T,fixb=F))
      res_para_new  = m2.mixt.movers(jdata,res_para_fixm$model,ctrl=em.control(ctrl,textapp=sprintf("para1 (%i/%i)",i,ctrl$est_rep),cstr_type="para",fixm=F,fixb=F))
      res_mixt      = m2.mixt.movers(jdata,res_para_new$model,ctrl=em.control(ctrl,textapp=sprintf("move1 (%i/%i)",i,ctrl$est_rep)))

      # ------ compute connectedness ----- #
      res_mixt$connectedness = model.connectiveness(res_mixt$model)
      res_mixt$rep_id = i
    }, error = function(e) {catf("error in rep %i!\n",i);print(e);})
    flog.info("done with reptitions %i/%i",i,ctrl$est_rep)
    res_mixt
  })

  # backing up to disk
  save(rr,ctrl,file=paste(ctrl$file_backup_prefix,"data",sep="."))

  # extract likelihoods and connectedness
  rrd = ldply(rr,function(r) {
    data.frame(lik_mixt = r$model$likm,connectedness = r$connectedness,i=r$rep_id)
  })

  # selecting best starting value
  rrd     = data.table(rrd)
  rrd[, sel:=-1]
  rrd.sub = rrd[order(-lik_mixt)][1:ctrl$est_nbest]
  rrd[i %in% rrd.sub$i, sel:=0]
  Ibest = rrd.sub[order(-connectedness)][1,i]
  res_mixt = rr[[Ibest]]
  rrd[i==Ibest, sel:=1]

  # sub-sample the stayers for computational reasons (if too large)
  if (ctrl$sdata_subredraw==TRUE) {
    sim$sdata[,sample := rank(runif(.N))/.N<=ctrl$sdata_subsample,j1]
    flog.info("drawing %f from the stayers",ctrl$sdata_subsample)
  }

  res_mixt$model = m2.mixt.stayers(sim$sdata[sample==1],res_mixt$model,ctrl = em.control(ctrl,textapp="stayers"))
  res_mixt$second_stage_reps = rrd
  res_mixt$second_stage_reps_all = rr

  # ------ compute linear decomposition ------- #
  NNm = res_mixt$model$NNm
  NNs = res_mixt$model$NNs/ctrl$sdata_subsample
  NNm[!is.finite(NNm)]=0
  NNs[!is.finite(NNs)]=0
  share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
  share_m  = sum(NNm)/(sum(NNm) + sum(NNs))

  NNs = round(NNs*ctrl$vdec_sim_size*share_s/sum(NNs))
  NNm = round(NNm*ctrl$vdec_sim_size*share_m/sum(NNm))

  # we simulate from the model both movers and stayers
  sdata.sim = m2.mixt.simulate.stayers(res_mixt$model,NNs)
  jdata.sim = m2.mixt.simulate.movers(res_mixt$model,NNm)
  sdata.sim = rbind(sdata.sim[,list(j1,k,y1)],jdata.sim[,list(j1,k,y1)])
  vdec  = lin.proj(sdata.sim,y_col = "y1",k_col="k",j_col = "j1")

  res_mixt$vdec = vdec
  res_mixt$ctrl = ctrl
  end.time <- Sys.time()

  res_mixt$time.taken <- end.time - start.time

  return(res_mixt)
}


#' Computes the variance decomposition by simulation
#' @export
m2.mixt.vdec <- function(model,nsim,stayer_share=1,ydep="y2") {

  if (ydep!="y1") flog.warn("ydep other than y1 is not implemented, using y1")

  # simulate movers/stayers, and combine
  NNm = model$NNm
  NNs = model$NNs
  NNm[!is.finite(NNm)]=0
  NNs[!is.finite(NNs)]=0

  NNs = round(NNs*nsim*stayer_share/sum(NNs))
  NNm = round(NNm*nsim*(1-stayer_share)/sum(NNm))
  flog.info("computing var decomposition with ns=%i nm=%i",sum(NNs),sum(NNm))


  # we simulate from the model both movers and stayers
  sdata.sim = m2.mixt.simulate.stayers(model,NNs)
  jdata.sim = m2.mixt.simulate.movers(model,NNm)
  sdata.sim = rbind(sdata.sim[,list(j1,k,y1)],jdata.sim[,list(j1,k,y1)])
  proj_unc  = lin.proj(sdata.sim,"y1","k","j1");

  return(proj_unc)
}


#' Compute mean effects
#' @export
m2.mixt.meaneffect <- function(model) {
  NNs = model$NNs*10 # used 10% sample
  NNm = model$NNm
  share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
  share_m  = sum(NNm)/(sum(NNm) + sum(NNs))

  NNs = round(NNs*1e6*share_s/sum(NNs))
  NNm = round(NNm*1e6*share_m/sum(NNm))

  # we simulate from the model both movers and stayers
  sdata = m2.mixt.simulate.stayers(model,NNs)
  jdata = m2.mixt.simulate.movers(model,NNm)
  sdata = rbind(sdata[,list(j1,k,y1)],jdata[,list(j1,k,y1)])

  # compute decomposition
  #vdec  = lin.proj(sdata,y_col = "y1",k_col="k",j_col = "j1")
  #res_bs$mixt_all[[nn]]$vdec_1m = vdec
  rt    = sample.stats(sdata,"y1","pk0")

  # then we set the distribution to uniform
  model_nosort = copy(model)
  model_nosort$pk0[1,,] = spread(colSums(model_nosort$pk0[1,,] * spread(NNs/(sum(NNs)),2,model_nosort$nk)),1,model_nosort$nf)

  # movers
  dpk1 = m2.get.pk1(model)
  pk   = dpk1[,pr_k[1],k][,V1]
  model_nosort$pk1      = spread(pk,1,model$nf * model$nf)

  # simulate from uniform
  sdata = m2.mixt.simulate.stayers(model_nosort,NNs)
  jdata = m2.mixt.simulate.movers(model_nosort,NNm)
  sdata = rbind(sdata[,list(j1,k,y1)],jdata[,list(j1,k,y1)])
  rt2   = sample.stats(sdata,"y1","pku")

  return(rbind(rt,rt2))
}


# ------------- Testing functions ---------------------

# for more tests, look at tests/testthat/test_model_mixt2.R

# here we want to check a bunch of properties for the EM steps
# model1 and model2 should be 2 consecutive steps
m2.mixt.check <- function(Y1,Y2,J1,J2,JJ,nk,Nm,model1,...) {

  change = list(...)

  # compute posterior for model1
  res1 = with(model1,{
    taum = array(0,c(Nm,nk))
    lpm  = array(0,c(Nm,nk))
    likm = 0
    for (i in 1:Nm) {
      ltau = log(pk1[JJ[i],])
      lnorm1 = lognormpdf(Y1[i], A1[J1[i],], S1[J1[i],])
      lnorm2 = lognormpdf(Y2[i], A2[J2[i],], S2[J2[i],])
      lall = ltau + lnorm2 + lnorm1
      lpm[i,] = lall
      likm = likm + logsumexp(lall)
      taum[i,] = exp(lall - logsumexp(lall))
    }

    # compute prior
    lik_prior = (dprior-1) * sum(log(pk1)) # dirichlet distribution
    lik = likm + lik_prior
    list(taum = taum, lpm =lpm, lik=likm,lik_prior=lik_prior,post=lik)
  })

  model2 = copy(model1)
  model2[names(change)] = change[names(change)]

  # compute posterior for model2
  res2 = with(model2,{
    taum = array(0,c(Nm,nk))
    lpm  = array(0,c(Nm,nk))
    likm = 0
    for (i in 1:Nm) {
      ltau = log(pk1[JJ[i],])
      lnorm1 = lognormpdf(Y1[i], A1[J1[i],], S1[J1[i],])
      lnorm2 = lognormpdf(Y2[i], A2[J2[i],], S2[J2[i],])
      lall = ltau + lnorm2 + lnorm1
      lpm[i,] = lall
      likm = likm + logsumexp(lall)
      taum[i,] = exp(lall - logsumexp(lall))
    }

    # compute prior
    lik_prior = (dprior-1) * sum(log(pk1)) # dirichlet distribution
    lik = likm + lik_prior
    list(taum = taum, lpm =lpm, lik=likm,lik_prior=lik_prior,post=lik)
  })

  # do the analysis, Evaluate Q(theta | theta^t) , Q(theta^t | theta^t), H(theta | theta^t) and H(theta^t | theta^t)
  Q1 = sum( ( (res1$taum) * res1$lpm ))
  Q2 = sum( ( (res1$taum) * res2$lpm ))
  H1 = - sum( (res1$taum) * log(res1$taum))
  H2 = - sum( (res1$taum) * log(res2$taum))

  warn_str=""
  test = TRUE
  if (( Q2<Q1) | (H2<H1)) {
    warn_str = "!!!!!!!!!";
    test=FALSE
  }
  catf("[emcheck] %s Qd=%4.4e Hd=%4.4e %s\n",paste(names(change),collapse = ","),  Q2-Q1,H2-H1,warn_str)
  return(test)
}


m2.mixt.test <- function() {
  nf = 10
  nk = 6
  model = m2.mixt.new(nk,nf)
  NNm = floor(array(30000/(nf^2),c(nf,nf)))
  jdata = m2.mixt.simulate.movers(model,NNm)

  ctrl   = em.control(nplot=10,check_lik=F,fixb=F,est_rho=F,model0=model,dprior=1.05,maxiter=100)
  ctrl$posterior_reg=0
  ctrl$fixm=FALSE
  ctrl$ncat=5
  ctrl$check_lik=FALSE

  res = m2.mixt(jdata,model,ctrl)

  # trying to do the no from there with 3 components
  ctrl$model0=NA
  model_np = step2.static.em.np.new.from.ns(res$model,nm=3)
  res =   model_np = step2.static.em.np.movers.estimate(jdata,model_np,ctrl)

  # try to plot the outcome.

  res$model$W1


  res = m2.mixt.fixed(jdata,model)
}

em.endo.simulatebest <- function() {

  # load the grid
  load("inst/figures/src/em-endo-full_rhogrid-halton-6x10.dat",verbose=F)

  # find the best
  dd = data.frame()
  for (r in rr) {
    dd = rbind(dd,data.frame(rho=r$model$B32m,lik=r$lik,time=r$time.taken,step=r$step,dlik=r$dlik))
  }
  rbest = rr[[which.max(dd$lik)]]
  cat(sprintf("%i evaluations, best value is %f\n",length(rr),rbest$lik))

  # get number of movers
  load("../figures/src/em-endo-info.dat",verbose=F)

  # reweight the statyers to 30,0000
  tot = NNs[,sum(ni)]
  NNs[,ni := round(ni * 30000 /sum(ni)) ]
  setkey(NNs,j)
  NNs = NNs[,ni]

  # get the movers matrix
  NNm = acast(NNm,j1~j2,value.var="ni")

  # ----- simulate ------ #
  nk = rbest$model$nk;
  nf = rbest$model$nf;
  model = rbest$model
  jdatae = em.endo.full.simulate.movers(model,NNm)
  sdatae = em.endo.full.simulate.stayers(model,NNs)
  jdatae[,m:=1][,w:=1]
  sdatae[,m:=0][,w:=tot/.N]
  sdatae[,j2:=j1]
  adatae = rbind(sdatae,sdatae)
  cat(sprintf("simulated data with %i stayers and %i movers \n",sdatae[,.N],jdatae[,.N]))
  return(adatae)
}

#' @export
m2.get.pk1 <- function(model) {
  pk1 = rdim(model$pk1,model$nf,model$nf,model$nk)
  dd_post = data.table(melt(pk1,c('j1','j2','k')))
  pp = model$NNm/sum(model$NNm)
  dd_post <- dd_post[, pr_j1j2 := pp[j1,j2],list(j1,j2)  ]
  dd_post <- dd_post[, pr_j1j2k := pr_j1j2*value]
  dd_post <- dd_post[, pr_k := sum(pr_j1j2k),k]
  dd_post
}


#' Returns the uconditional type probability in the crossection
#' @export
m2.get.pk_unc <- function(model,supersample=0.1) {
  dpk1     = m2.get.pk1(model)
  pk_m     = acast(dpk1[,sum(pr_j1j2k),list(j1,k)],j1~k,value.var = "V1")
  NNs      = model$NNs*round(1/supersample) # used 10% sample
  NNm      = model$NNm
  share_s  = sum(NNs)/(sum(NNm) + sum(NNs))
  pk_unc   = share_s*rdim(res_main$model$pk0[,,I],res_main$model$nf,res_main$model$nk) +
    (1- share_s) * pk_m
  pk_unc
}

#' check the fit in the movers/stayers using imputed data
#' @export
m2.movers.checkfit <- function(jdata) {
  dd = jdata[, {
    d1=data.frame(src="data",
                  m1=mean(y1),m2=mean(y2),
                  d12=mean(y1-y2),
                  cov12=cov(y1,y2),v1=var(y1),v2=var(y2))
    d2=data.frame(src="imp",
                  m1=mean(y1_imp),m2=mean(y2_imp),
                  d12=mean(y1_imp-y2_imp),
                  cov12=cov(y1_imp,y2_imp),v1=var(y1_imp),v2=var(y2_imp))
    rbind(d1,d2)
  },list(j1,j2)]

  ddm = melt(dd,id.vars = c("j1","j2","src"))
  ddm = cast(ddm,j1+j2+variable~src,value = "value")
  ggplot(ddm,aes(x=data,y=imp)) + geom_point() + facet_wrap(~variable,scales = "free") + theme_bw() +
    geom_abline(linetype=2)
  ddm
}

#' check the fit in the movers/stayers using imputed data
#' @export
m2.stayers.checkfit <- function(sdata,r1,r4) {
  dd = jdata[, {
    d1=data.frame(src="data",
                  m1=mean(y1),m2=mean(y2),
                  d12=mean(y1-y2),
                  cov12=cov(y1,y2),v1=var(y1),v2=var(y2))
    d2=data.frame(src="imp",
                  m1=mean(y1_imp),m2=mean(y2_imp),
                  d12=mean(y1_imp-y2_imp),
                  cov12=cov(y1_imp,y2_imp),v1=var(y1_imp),v2=var(y2_imp))
    rbind(d1,d2)
  },list(j1)]

  ddm = melt(dd,id.vars = c("j1","src"))
  ddm = cast(ddm,j1+variable~src,value = "value")
  ggplot(ddm,aes(x=data,y=imp)) + geom_point() + facet_wrap(~variable,scales = "free") + theme_bw() +
    geom_abline(linetype=2)
  ddm
}
