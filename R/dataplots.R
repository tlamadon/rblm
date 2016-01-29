# pre-processing of the data
# 1) keep only stayers and movers
# 2) regress wages on agedum + agedum * period2 , remove interaction only


#' plot the mean wage at origin conditional on where it is coming from
plot.wage <- function(jdata) {
  
  rr = jdata[, {
    fit = lm(y2 ~ bs(y1,df=10), .SD)
    xnew = quantile(y1,seq(0.05,0.95,l=100))
    list(y1=y1,y2=predict(fit))
  },list(j1,j2)]
  
  rr = jdata[,mean(y2),list(j1,j2)]
  ggplot(rr,aes(x=j1,y=V1,color=j2 , group=j2)) + geom_line() + theme_bw()
  
  rr = jdata[,mean(y2-y1),list(j1,j2)]
  ggplot(rr,aes(x=j1,y=V1,color=j2 , group=j2)) + geom_line() + theme_bw()
  
  # check education / gender / age per cluster
  
  rr = jdata[,.N,]
  
  
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}