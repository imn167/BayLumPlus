#Mock data creations helpers


test_that("Full function test", {
dt <- list(D = OSLJingbian$D, sD = OSLJingbian$sD, ddot = OSLJingbian$ddot)
  expect_s3_class(suppressWarnings(Compute_AgeS_D(DATA = dt, Nb_sample = OSLJingbian$Nb_Sample,
                                  SampleNames = OSLJingbian$SampleNames, ThetaMatrix = OSLJingbian$ThetaMatrix,PriorAge = rep(c(1, 1400), OSLJingbian$Nb_Sample), jags_method = "rjags",
                                 prior = "unconstrained_Jeffrey", Iter = 100, adapt = 100, burnin = 50, n.chains = 2, quiet = TRUE)),
                  class = "BayLum.list")

  expect_s3_class(suppressWarnings(Compute_AgeS_D(DATA = dt, Nb_sample = OSLJingbian$Nb_Sample,
                                                  SampleNames = OSLJingbian$SampleNames, ThetaMatrix = OSLJingbian$ThetaMatrix,PriorAge = rep(c(1, 1400), OSLJingbian$Nb_Sample), jags_method = "rjags",
                                                  prior = "constrained_Jeffrey", Iter = 100, adapt = 100, burnin = 50, n.chains = 2, quiet = TRUE)),
                  class = "BayLum.list")

  })



