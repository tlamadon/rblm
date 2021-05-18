# test_that("testing trace code", {
#   set.seed(324313)
#   model = m2.mini.new(10,serial = F,fixb=T)

#   # we set the parameters to something simple
#   model$A1 = seq(0,2,l=model$nf) # setting increasing intercepts
#   model$B1 = seq(1,2,l=model$nf) # adding complementarity (increasing interactions)
#   model$Em = seq(0,1,l=model$nf) # adding sorting (mean type of workers is increasing in k)

#   # we make the model stationary (same in both periods)
#   model$A2 = model$A1
#   model$B2 = model$B1

#   # setting the number of movers and stayers
#   model$Ns   = array(300000/model$nf,model$nf)
#   model$Nm   = 10*toeplitz(ceiling(pmax(seq(100,-50,l=model$nf),0)  ))

#   # creating a simulated data set
#   ad =  m2.mini.simulate(model)

#   # clsutering
#   ms    = grouping.getMeasures(ad,"ecdf",Nw=10,y_var = "y1")
#   grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=250)
#   ad   = grouping.append(ad,grps$best_cluster,drop=T)

#   # try the mini-model
#   res = m2.mini.estimate(ad$jdata,ad$sdata,model0 = model,method = "linear.ss")

#   # testing the penalized AKM
#   for (lambda in 10^(-5:4)) {
#     m2.firmfe.pen(ad,model,lambda=0.00001) # FIXME changed from m2.firmfe.pen(sim,mode,lambda=0.00001)
#   }
#   m2.trace.estimate(ad) # FIXME changed from m2.trace.estimate(sim)

#   expect_equal(0, 0)
#   # expect_that(floor_base("year"), is_time("2009-01-01 00:00:00"))
# })


