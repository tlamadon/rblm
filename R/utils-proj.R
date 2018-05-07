#' Generate a linear projection decomposition for the model
#' with continuous worker hetergoneity
#'
#' @export
lin.proja <- function(sdata,y_col="y",k_col="k",j_col="j") {
  rr = list()

  sdata2 = copy(data.table(sdata))
  sdata2[,y_imp := get(y_col)]
  sdata2[,k_imp := get(k_col)]
  sdata2[,j     := get(j_col)]

  fit = lm(y_imp ~ k_imp + factor(j),sdata2)
  sdata2$res = residuals(fit)
  pred = predict(fit,type = "terms")
  sdata2$k_hat = pred[,1]
  sdata2$l_hat = pred[,2]

  rr$cc = sdata2[,cov.wt( data.frame(y_imp,k_hat,l_hat,res))$cov]
  rr$rsq1 = summary(fit)$r.squared

  fit2 = lm(y_imp ~ 0+  k_imp:factor(j) + factor(j),sdata2)
  rr$rsq2 = 1-mean(resid(fit2)^2)/var(sdata2$y_imp)

  get.stats <- function(cc) {
    r=list()
    den = cc[2,2] + cc[3,3] + 2 * cc[2,3]
    r$cor_kl = round(cc[2,3]/sqrt(cc[2,2]*cc[3,3]),4)
    r$cov_kl = 2*round(cc[2,3]/den,4)
    r$var_k  = round(cc[2,2]/den,4)
    r$var_l  = round(cc[3,3]/den,4)
    r$rsq    = round((cc[1,1] - cc[4,4])/cc[1,1],4)
    return(r)
  }

  rr$stats = get.stats(rr$cc)
  print(data.frame(rr$stats))
  rr$NNs = sdata[,.N,j1][order(j1)][,N]

  return(rr)
}



#' @export
lin.proj <- function(sdata,y_col="y",k_col="k",j_col="j",usex=FALSE,do.unc=TRUE) {
  rr = list()

  sdata2 = copy(data.table(sdata))
  sdata2[,y_imp := get(y_col)]
  sdata2[,k_imp := get(k_col)]
  sdata2[,j     := get(j_col)]

  fit = lm(y_imp ~ factor(k_imp) + factor(j),sdata2)
  sdata2$res = residuals(fit)
  pred = predict(fit,type = "terms")
  sdata2$k_hat = pred[,1]
  sdata2$l_hat = pred[,2]

  rr$cc = sdata2[,cov.wt( data.frame(y_imp,k_hat,l_hat,res))$cov]
  rr$rsq1 = summary(fit)$r.squared

  if (do.unc) {
    fit2 = lm(y_imp ~factor(k_imp) * factor(j),sdata2)
    rr$rsq2 = summary(fit2)$r.squared
  }

  get.stats <- function(cc) {
    r=list()
    den = cc[2,2] + cc[3,3] + 2 * cc[2,3]
    r$cor_kl = round(cc[2,3]/sqrt(cc[2,2]*cc[3,3]),4)
    r$cov_kl = 2*round(cc[2,3]/den,4)
    r$var_k  = round(cc[2,2]/den,4)
    r$var_l  = round(cc[3,3]/den,4)
    r$rsq    = round((cc[1,1] - cc[4,4])/cc[1,1],4)
    return(r)
  }

  rr$stats = get.stats(rr$cc)
  rr$NNs = sdata2[,.N,j][order(j)][,N]
  print(data.frame(rr$stats))

  return(rr)
}

#' Computes the linear projection using X
#' @export
lin.projx <- function(sdata,y_col="y",k_col="k",j_col="j") {
  rr = list()

  sdata2 = copy(data.table(sdata))
  sdata2[,y_imp := get(y_col)]
  sdata2[,k_imp := get(k_col)]
  sdata2[,j     := get(j_col)]

  fit = lm(y_imp ~ factor(k_imp) + factor(j) +factor(x),sdata2)
  sdata2$res = residuals(fit)
  pred = predict(fit,type = "terms")
  sdata2$k_hat = pred[,1]
  sdata2$l_hat = pred[,2]
  sdata2$x_hat = pred[,3]

  rr$cc = sdata2[,cov.wt( data.frame(y_imp,k_hat,l_hat,res,x_hat))$cov]

  fit2 = lm(y_imp ~factor(k_imp) * factor(j),sdata2)

  rr$rsq1 = summary(fit)$r.squared
  rr$rsq2 = summary(fit2)$r.squared

  get.stats <- function(cc) {
    r=list()
    den = cc[2,2] + cc[3,3] + 2 * cc[2,3]
    r$cor_kl = round(cc[2,3]/sqrt(cc[2,2]*cc[3,3]),4)
    r$cov_kl = 2*round(cc[2,3]/den,4)
    r$var_k  = round(cc[2,2]/den,4)
    r$var_l  = round(cc[3,3]/den,4)
    r$rsq    = round((cc[1,1] - cc[4,4])/cc[1,1],4)
    return(r)
  }

  rr$stats = get.stats(rr$cc)
  print(data.frame(rr$stats))

  return(rr)
}

#' Computes the linear projection using X
#' @export
lin.projax <- function(sdata,y_col="y",k_col="k",j_col="j") {
  rr = list()

  sdata2 = copy(data.table(sdata))
  sdata2[,y_imp := get(y_col)]
  sdata2[,k_imp := get(k_col)]
  sdata2[,j     := get(j_col)]

  fit = lm(y_imp ~ k_imp + factor(j) +factor(x),sdata2)
  sdata2$res = residuals(fit)
  pred = predict(fit,type = "terms")
  sdata2$k_hat = pred[,1]
  sdata2$l_hat = pred[,2]
  sdata2$x_hat = pred[,3]

  rr$cc = sdata2[,cov.wt( data.frame(y_imp,k_hat,l_hat,res,x_hat))$cov]
  rr$rsq1 = summary(fit)$r.squared

  get.stats <- function(cc) {
    r=list()
    den = cc[2,2] + cc[3,3] + 2 * cc[2,3]
    r$cor_kl = round(cc[2,3]/sqrt(cc[2,2]*cc[3,3]),4)
    r$cov_kl = 2*round(cc[2,3]/den,4)
    r$var_k  = round(cc[2,2]/den,4)
    r$var_l  = round(cc[3,3]/den,4)
    r$rsq    = round((cc[1,1] - cc[4,4])/cc[1,1],4)
    return(r)
  }

  rr$stats = get.stats(rr$cc)
  print(data.frame(rr$stats))

  return(rr)
}
