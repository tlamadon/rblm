require(rblm)

# ====== simulate some data ==========
model = model.mini2.new(10,serial = F)
NNs   = array(300000/model$nf,model$nf)
sdata = model.mini2.simulate.stayers(model,NNs)
jdata = model.mini2.simulate.movers(model,NNm)

# randomly assign firm IDs
sdata[,f1:=j1 + model$nf*(sample.int(.N/50,.N,replace=T)-1),j1]
sdata[,j1b:=j1]
jdata[,j1c:=j1]
jdata[,f1:=sample( unique(sdata[j1b==j1c,f1])    ,.N,replace=T),j1c]
jdata[,f2:=sample( unique(sdata[j1b==j2,f1])    ,.N,replace=T),j2]

# ========== clustering firms ================

#plot firm size distribution
sdata[,.N,f1][,hist(N)]

# rename to match cluster function names: fid, y
# we cluster using the wage in period 1, we need to call firm id "fid" and the wage "lw" hence we use sdata[,list(fid=f1,lw=y1)]
clus = cluster.firms.data(sdata[,list(fid=f1,lw=y1)],ncluster=10,nw=40,nstart=300,step=20,merge = F)

# reattach the clus to the orginal data in period
sdata[,j1t:=j1] # save true type
setkey(sdata,f1)
setkey(clus,fid)
sdata[,j1:=clus[sdata,clus]]

# order by mean wage
JJ = sdata[,mean(y1),j1][order(j1)][,rank(V1)]
sdata[,j1 := JJ[j1]]

# check that the firms are well classified, you should see a strong diagonal
rr = sdata[,list(j1t=j1t[1],j1=j1),f1][,.N,list(j1,j1t)]
ggplot(rr,aes(x=j1,y=j1t,size=N)) + geom_point()

# attach clsuter to movers
jdata[,j1t:=j1] # save true type
setkey(jdata,f1)
jdata[,j1:=clus[jdata,clus]]
jdata[,j2t:=j2] # save true type
setkey(jdata,f2)
jdata[,j2:=clus[jdata,clus]]



# ==== ESTIMATE THE MODEL ========

res = model.mini2.estimate(jdata,sdata,model0 = model)


# ==== simulating and estimating the 4 period model ========

# not working just yet!

model  = model.mini4.new(10,r1=0.4,r4=0.8)
NNm    = array(50000/(model$nf^2),c(model$nf,model$nf))
NNs    = array(100000/model$nf,model$nf)
jdata  = model.mini4.simulate.movers(model,NNm)
sdata  = model.mini4.simulate.stayers(model,NNs)
model1 = model.mini4.estimate(jdata,sdata,model$r1,model$r4,model0 = model);
