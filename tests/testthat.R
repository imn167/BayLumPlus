library(testthat)
library(BayLumPlus)

# Set number of parallel workers to available cores minus 1
library(parallel)
num_cores <- max(1, parallel::detectCores() - 1)
Sys.setenv("TESTTHAT_PARALLEL" = "true")
Sys.setenv("TESTTHAT_CPUS" = num_cores)

test_check("BayLumPlus")
