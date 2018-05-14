
#' Create a control structure for running EM algorithms
#'
#' @param nplot how often to plot wages of comparaison to model0
#' @param ncat how often to log information
#' @param maxiter maximum number of iterations
#' @param model_var whether to allow for flexible variance (default is TRUE)
#' @param est_Amb TBD
#' @param cstr_type defines the type of constraints on the means for estimation
#' @param cstr_val TBD
#' @param tol tolerance for stopping the EM
#' @param tau posterior likelihood to use instead of computing them using initial parameters
#' @param model0 model to compare estimation to when plotting
#' @param presp=1e-9 TBD
#' @param rel_weight setting a particular weight for stayers versus movers (default is 1)
#' @param est_rho  vector of TRUE/FALSE that states which rho should be estimated
#' @param rho_in_diff whether to estimate rho in level of in differences
#' @param dprior    Dirichlet prior for proportions (default 1.01)
#' @param nfirms_to_update=10, number firms to update in probabilistic approach
#' @param proba_include what terms to include in the liklihood for probabilistic approach (default to = c(1,1,1,1))
#' @param check_lik whether to check the likelihood at each updating parameter
#' @param stochastic whether to use stochastic EM instead of straight EM
#' @param fixb   when TRUE, imposes fixed interactions in different time periods
#' @param fixm when TRUE, levels are not updated, only variances and proportions
#' @param deps=1e-50 TBD
#' @param file_backup_prefix TBD
#' @param sd_floor= floor imposed on all standard deviations (default is 1e-10)
#' @param posterior_reg term added to posterior probablities (this is to deal with numerical issues, default is 1e-9)
#' @param textapp  text to show in logging
#' @param sdata_subsample share of the stayers to use in estimation
#' @param sdata_subredraw whether to redraw the subsample of stayers
#' @param vdec_sim_size  size to use in simulation
#' @param stayer_weight weight attributed to stayers in joint estimation
#' @param est_rep= number of starting values for EM
#' @param est_nbest number of best starting values to select from using connectedness
#'
#' @export
em.control <- function(ctrl=NULL,...) {
  args = list(...)
  if (is.null(ctrl)) {
    ctrl = list(
      nplot=5,
      ncat=25,
      maxiter=2000,
      model_var=TRUE,
      est_Amb=TRUE,
      cstr_type = "none",
      cstr_val = 0,
      tol=1e-9,
      tau=NA,
      model0=NA,
      presp=1e-9,
      rel_weight = 1,
      est_rho=FALSE,       # vector of TRUE/FALSE that states which rho should be estimated
      rho_in_diff=FALSE,   # whether to estimate rho in level of in differences
      dprior = 1.01,       # Dirichlet prior for proportions
      nfirms_to_update=10, # number firms to update in probabilistic approach
      proba_include = c(1,1,1,1), # what terms to include in the liklihood for probabilistic approach
      check_lik=FALSE,
      stochastic=0,
      fixb=FALSE,          # when TRUE, imposes fixed interactions in different time periods
      fixm=FALSE,          # when TRUE, levels are not updated, only variances and proportions
      deps=1e-50,
      file_backup_prefix = "estimation-bu-tmp",
      sd_floor=1e-10,      # floor imposed on all standard deviations
      posterior_reg=1e-9,  # term added to posterior probablities (this is to deal with numerical issues)
      textapp="",          # text to show in logs
      sdata_subsample=1.0, # share of the stayers to use in estimation
      sdata_subredraw=TRUE,# whether to redraw the subsample of stayers
      vdec_sim_size=1e6,        # size to use in simulation
      stayer_weight=1,     # weight attributed to stayers in joint estimation
      est_rep=10,          # number of starting values for EM
      est_nbest=5)         # number of best starting values to select from using connectedness
  }
  ctrl[names(args)]  = args[names(args)]

  # for the following argument, if only one, we repeat them
  for (nn in c('cstr_type','cstr_val','est_rho','est_Amb')) {
    if (length(ctrl[[nn]])==1) ctrl[[nn]]=rep(ctrl[[nn]],6);
  }

  return(ctrl)
}


#' functions for em
#' @export
lognormpdf <- function(Y,mu=0,sigma=1)   -0.5 * (  (Y-mu) / sigma )^2   - 0.5 * log(2.0*pi) - log(sigma)

#' logsumexp function
#' @export
logsumexp <- function(v) {
  vm = max(v)
  log(sum(exp(v-vm))) + vm
}

#' logsumexp function by Row
#' @export
logRowSumExp <- function(M) {
  if (is.null(dim(M))) {return(M)}
  vms = apply(M,1,max)
  log(rowSums(exp(M-spread(vms,2,dim(M)[2])))) + vms
}

