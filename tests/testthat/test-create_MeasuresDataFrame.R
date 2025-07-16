

test_that("Covariance Matrix for D", {
  DATA4 <-  combine_DataFiles(DATA1, DATA2)
  PalaeoObject <- Palaeodose_Computation(DATA4, c("samp1","samp2"), 2, Iter = 10000, n.chains = 2)
  # expect_equal(length(create_MeasuresDataFrame(PalaeoObject, DATA4, .2, c(.001, .058))), 2)
})
