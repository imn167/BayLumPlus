############################################################@
#       Script for rading Jags Models
# Author : Bouafia Imene
# Date : february 18 2025
####################################################@#

##cleaning
rm(list = ls())

library(here)
allModels <- list.files(here("data"))
allModels <- allModels[grep(".rda", allModels)]

for (f in allModels) {
  path <- paste0("data/", f)
load(here(path))
}

cat(Model_AgeS$AgesMultiCS2_EXPZO$cauchy)

cat(Model_Palaeodose$PalaeodosesMultiBF_EXPLIN$cauchy)
