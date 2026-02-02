test_that("multiplication works", {
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
                                                     prior = "NichollsJones", Iter = 200, n.chains = 2,
                                                     jags_method = "rjags")), class = "BayLum.list")
})
