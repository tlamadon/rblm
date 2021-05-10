test_that("mini model with lots of groups", {
  set.seed(324313)
  model = m2.mini.new(50,serial = F,fixb=T)

  # we set the parameters to something simple
  model$A1 = seq(0,2,l=model$nf) # setting increasing intercepts
  model$B1 = seq(1,2,l=model$nf) # adding complementarity (increasing interactions)
  model$Em = seq(0,1,l=model$nf) # adding sorting (mean type of workers is increasing in k)

  # we make the model stationary (same in both periods)
  model$A2 = model$A1
  model$B2 = model$B1

  # setting the number of movers and stayers
  model$Ns   = array(300000/model$nf,model$nf)
  model$Nm   = 2*toeplitz(ceiling(pmax(seq(100,-50,l=model$nf),0)  ))

  # creating a simulated data set
  ad =  m2.mini.simulate(model,vdec=FALSE)

  # clsutering
  ms    = grouping.getMeasures(ad,"ecdf",Nw=10,y_var = "y1")
  # then we group we choose k=10
  grps  = grouping.classify.once(ms,k = 100,nstart = 500,iter.max = 200,step=250)

  # finally we append the results to adata
  ad   = grouping.append(ad,grps$best_cluster,drop=T)

  acast(ad$jdata[,.N,list(j1,j2)],j1~j2,fill=0)
  # remove all of a given cluster in period in jdata
  ad$jdata[j1==4,j1:=5]

  # try the mini-model
  res = m2.mini.estimate(ad$jdata,ad$sdata,model0 = model,method = "prof",bigk = 1,msub=0.1)

  expect_that(floor_base("year"), is_time("2009-01-01 00:00:00"))
})


