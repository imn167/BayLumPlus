

if (capabilities("cairo")) {
  options(bitmapType = "cairo")
}

# For headless testing
if (!interactive()) {
  options(device = "pdf")
}

#Mock data creations helpers
create_mock_object <-  function(n_samples=5) {
  mock_samples <- matrix(abs(rnorm(n_samples*100)), ncol = n_samples)
  sample_names <- paste0("Sample", 1:n_samples)
  colnames(mock_samples) <- sample_names

  #returned
  .list_BayLum(
    "Sampling" = coda::as.mcmc.list(coda::as.mcmc(mock_samples)),
    "Ages" = data.frame(SAMPLE = sample_names, stringsAsFactors = F, HPD95.MIN = rep(10, n_samples), HPD95.MAX = rep(20, n_samples),
                        HPD68.MIN = rep(12, n_samples), HPD68.MAX = rep(18, n_samples),
                        AGE = rep(15, n_samples), Unit = 1:n_samples) %>% dplyr::mutate(`Bayes sd` =rep(1, n_samples) ),
    "Summary" = as.data.frame(replicate(2, rep(1, n_samples))) %>% dplyr::mutate(`Bayes sd` =rep(1, n_samples) ) ,
    prior = "unconstrained_Jeffrey"
  )
}

create_mock_constraints <- function(n = 5) {
  # Create simple stratigraphic constraint matrix
  matrix(c(rep(1, n), upper.tri(matrix(1, ncol = n, nrow = n)) * 1),
         nrow = 6, byrow = TRUE)
}

create_mock_isoObject <- function(n_samples = 5) {
  List <- create_mock_object(n_samples)
  List[[3]] <-  buildNetwork(create_mock_constraints(n_samples))
  return(List)
}


test_that("IsotonicCurve runs with valid inputs", {
  mock_obj = create_mock_object()
  mock_constraints = create_mock_constraints()

  expect_no_error({result <- IsotonicCurve(mock_constraints, mock_obj, levels = c(.95), interactive = F)})

  ## Null object

  expect_error(IsotonicCurve(mock_constraints, NULL, interactive = F), class = "error")
})


test_that("function handles different StratiConstraints formats", {

  mock_obj <- create_mock_object()
  matrix_constraint <- create_mock_constraints()
  expect_no_error(IsotonicCurve(c(), mock_obj, levels = c(.95), interactive =F))


  temp_csv <- tempfile(fileext = ".csv")
  write.csv(matrix_constraint, temp_csv, row.names = FALSE)

  expect_no_error(
    IsotonicCurve(temp_csv, mock_obj, levels = 0.95, interactive = FALSE)
  )
})

test_that("function handles different levels correctly", {
  skip_if_not_installed("runjags")

  mock_obj <- create_mock_object()
  mock_constraints <- create_mock_constraints()

  # Single level
  expect_no_error(
    IsotonicCurve(mock_constraints, mock_obj, levels = 0.95, interactive = FALSE)
  )

  # Multiple levels
  expect_no_error(
    IsotonicCurve(mock_constraints, mock_obj, levels = c(0.65, 0.95), interactive = FALSE)
  )

  # Invalid levels (outside 0-1 range)
  expect_warning(expect_error(
    IsotonicCurve(mock_constraints, mock_obj, levels = 1.5, interactive = FALSE)
  ))


})

### Using other prior than unconstrained_jeffrey
test_that("restriction to unconstrained priors", {
  skip_if_not_installed("runjags")

  mock_obj <- create_mock_object()
  mock_constraints <- create_mock_constraints()

  mock_obj$prior = "other_prior"
  expect_error(
    IsotonicCurve(mock_constraints, mock_obj, interactive = FALSE)
  )
})


### PlotIsotonicCurve tests

test_that("Restriction to one HPD region / Credible Interval level only", {
  mock_object <- create_mock_isoObject()

  # Test with a calculated level
    expect_no_error(PlotIsotonicCurve(mock_object, .95))

  # Test with an uncalculated level

    expect_error(PlotIsotonicCurve(mock_object, .62))

    # Test with several levels

    expect_error(PlotIsotonicCurve(mock_object, c(.95, .68)))

})























