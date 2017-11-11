require(rblm)

# ====== simulate some data ==========
model = model.mini2.new(10,serial = F)
NNs   = array(300000/model$nf,model$nf)
NNm   = array(30000/model$nf^2,c(model$nf,model$nf))
sdata = model.mini2.simulate.stayers(model,NNs)
jdata = model.mini2.simulate.movers(model,NNm)

# randomly assign firm IDs
sdata[,f1:=paste("F",j1 + model$nf*(sample.int(.N/50,.N,replace=T)-1),sep=""),j1]
sdata[,j1b:=j1]
jdata[,j1c:=j1]
jdata[,f1:=sample( unique(sdata[j1b==j1c,f1])    ,.N,replace=T),j1c]
jdata[,f2:=sample( unique(sdata[j1b==j2,f1])    ,.N,replace=T),j2]

# combine the movers and stayers, ad stands for all data:
ad = list(sdata=sdata,jdata=jdata)

# plot firm size distribution
ad$sdata[,.N,f1][,hist(N)]

# ========== clustering firms ================

# we start by extracting the measures that will be used to cluster
ms    = grouping.getMeasures(ad,"ecdf",Nw=20,y_var = "y1")
# in the previous command we tell rblm that we want to use the firm
# specific empirical measure "ecdf" with 20 points of supports and that
# the dependent variable is "y1". The firm identifier should be "f1"

# then we group we choose k=10
grps  = grouping.classify.once(ms,k = 10,nstart = 1000,iter.max = 200,step=100)

# finally we append the results to adata
ad   = grouping.append(ad,grps$best_cluster)

# ==== ESTIMATE THE MODEL (OLD CODE) ========
res = model.mini2.estimate(ad$jdata,ad$sdata,model0 = model)

# ==== Extract all necesserary moments for off-line estimation ===
mstats = ad$jdata[,list(m1=mean(y1),sd1=sd(y1),
                         m2=mean(y2),sd2=sd(y2),
                         v12 = cov(y1,y2),.N),list(j1,j2)]
cstats = ad$sdata[,list(m1=mean(y1),sd1=sd(y1),
                         m2=mean(y2),sd2=sd(y2),
                         v12 = cov(y1,y2),.N),list(j1)]

# save these for later estimation


