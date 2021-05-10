test_that("testing EM algorithm Event Study", {
  # Set parameters
  set.seed(1234)
  nn = 50000
  nl = 3
  nk = 4

  tol = 1e-4
  
  # Simulate the model  
  model = m2.mixt.new(nk=nl,nf=nk,fixb=F,stationary=F)

  # Reduce variances
  model$S1 = model$S1 / 10
  model$S2 = model$S2 / 10

  # Simulate movers
  movers = m2.mixt.simulate.movers(model)

  # Estimate the model using the simulated mover data, and the truth as starting values
  model_copy = duplicate(model, shallow=FALSE)

  res = m2.mixt.movers(jdata, model_start, ctrl=em.control(ctrl=NULL, maxiter=1, cstr_type="para", textapp="para0", fixb=F))
  res_model = res$model

  expect_gt(mean((model$A1 - res_model$A1) ^ 2), 1e-20)
  expect_lt(mean((model$A1 - res_model$A1) ^ 2), tol)

  expect_gt(mean((model$S1 - res_model$S1) ^ 2), 1e-20)
  expect_lt(mean((model$S1 - res_model$S1) ^ 2), tol)

  expect_gt(mean((model$A2 - res_model$A2) ^ 2), 1e-20)
  expect_lt(mean((model$A2 - res_model$A2) ^ 2), tol)

  expect_gt(mean((model$S2 - res_model$S2) ^ 2), 1e-20)
  expect_lt(mean((model$S2 - res_model$S2) ^ 2), tol)

  expect_gt(mean((model$pk1 - res_model$pk1) ^ 2), 1e-20)
  expect_lt(mean((model$pk1 - res_model$pk1) ^ 2), tol)
})
