test_that("hybrid function OSL - C14 - UTh", {
  D = c(17.72383, 104.95460 )
  sD = c(0.6, 4)
  M = c(D[1], 25, DATA_C14$C14[1,1], D[2]) # measures
  DATA3 <- combine_DataFiles(L1 = DATA2, L2 = DATA1)
  Nb_sample = 4
  SampleNames = c("GDB5", "UTh", DATA_C14$Names[1], "GDB3")
  Theta = diag(c(DATA3$ddot_env[1,1]**2 * .1, DATA3$ddot_env[1,2]**2 * .1))
  encoding = c(1, 2,0,1)

  ## DATA input
  dt = list(M = M, sD = sD, ddot = DATA3$ddot_env[1,], sigmaC14Cal =  DATA_C14$C14[1,2], sigma_UTh = 3)

  expect_s3_class(suppressWarnings(Compute_AgeS_DC14(dt, Nb_sample, SampleNames, encoding, Theta,
                                                     prior = "Unconstrained", Iter = 200, n.chains = 2,
                                                     jags_method = "rjags")), class = "BayLum.list")
})


test_that("hybrid function C14 - UTh", {
  M = c( 25, DATA_C14$C14[1,1]) # measures
  Nb_sample = 2
  SampleNames = c("UTh", DATA_C14$Names[1])
  Theta = NULL
  encoding = c(2,0)

  ## DATA input
  dt = list(M = M, sigmaC14Cal =  DATA_C14$C14[1,2], sigma_UTh = 3)

  expect_s3_class(suppressWarnings(Compute_AgeS_DC14(dt, Nb_sample, SampleNames, encoding, Theta,
                                                     prior = "Unconstrained", Iter = 200, n.chains = 2,
                                                     jags_method = "rjags")), class = "BayLum.list")
})

test_that("hybrid function OSL - UTh", {
  D = c(17.72383, 104.95460 )
  sD = c(0.6, 4)
  M = c(D[1], 25, D[2]) # measures
  DATA3 <- combine_DataFiles(L1 = DATA2, L2 = DATA1)
  Nb_sample = 3
  SampleNames = c("GDB5", "UTh",  "GDB3")
  Theta = diag(c(DATA3$ddot_env[1,1]**2 * .1, DATA3$ddot_env[1,2]**2 * .1))
  encoding = c(1, 2,1)

  ## DATA input
  dt = list(M = M, sD = sD, ddot = DATA3$ddot_env[1,], sigma_UTh = 3)

  expect_s3_class(suppressWarnings(Compute_AgeS_DC14(dt, Nb_sample, SampleNames, encoding, Theta,
                                                     prior = "Unconstrained", Iter = 200, n.chains = 2,
                                                     jags_method = "rjags")), class = "BayLum.list")
})


test_that("hybrid function OSL - C14", {
  D = c(17.72383, 104.95460 )
  sD = c(0.6, 4)
  M = c(D[1], DATA_C14$C14[1,1], D[2]) # measures
  DATA3 <- combine_DataFiles(L1 = DATA2, L2 = DATA1)
  Nb_sample = 3
  SampleNames = c("GDB5", DATA_C14$Names[1], "GDB3")
  Theta = diag(c(DATA3$ddot_env[1,1]**2 * .1, DATA3$ddot_env[1,2]**2 * .1))
  encoding = c(1,0,1)

  ## DATA input
  dt = list(M = M, sD = sD, ddot = DATA3$ddot_env[1,], sigmaC14Cal =  DATA_C14$C14[1,2])

  expect_s3_class(suppressWarnings(Compute_AgeS_DC14(dt, Nb_sample, SampleNames, encoding, Theta,
                                                     prior = "Unconstrained", Iter = 200, n.chains = 2,
                                                     jags_method = "rjags")), class = "BayLum.list")
})








