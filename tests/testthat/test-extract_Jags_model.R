
skip_if(!interactive(), "skipping interactive test")

test_that("No running if not interactive", {
 expect_equal(extract_Jags_model(), NULL)
})


