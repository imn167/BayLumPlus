test_that("output_of_the_measuresData ", {
  DATA4 <- combine_DataFiles(L1 = DATA2, L2 = DATA1)
  expect_s3_class(suppressWarnings(Palaeodose_Computation(DATA4,
                                 Nb_sample = 2, Iter = 1000, SampleNames = DATA4$SampleNames)),
                  class = "BayLum.list")
  create_MeasuresDataFrame(P, DATA4, .01, c(.2, .3))
})
