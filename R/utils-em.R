
#' Create a control structure for running EM algorithms
#'
#' fixb will make the interactions stationary
#' fixm will keep the means fixed to the starting value
#' stochastic=n will draw latent variables instead of computing all probabilities
#' check_lik will compute H and Q function at every step, for each parameter update
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

