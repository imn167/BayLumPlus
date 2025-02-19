############################################################@
#       Script for rading Jags Models
# Author : Bouafia Imene
# Date : february 18 2025
####################################################@#

##cleaning
rm(list = ls())

library(here)
allModels <-  list.files(here("data"))[grep("Model", list.files(here("data"))) ]

for (f in allModels) {
  path <- paste0("data/", f)
load(here(path))
}

## Age Estimation
cat(Model_AgeS$AgesMultiCS2_EXPZO$gaussian)

cat(Model_Palaeodose$PalaeodosesMultiBF_EXPLIN$cauchy)

## BayLum Prior with Jeffreys priors ?
cat(ModelPrior$OSL)



