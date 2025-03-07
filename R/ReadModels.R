############################################################@
#       Script for rading Jags Models
# Author : Bouafia Imene
# Date : february 18 2025
####################################################@#
## global variables
Item <- NULL
#' @title Extract and View Bayesian Models in BayLum
#'
#' @description
#' This function interacts with the user to display a list of all available Bayesian models in the **BayLum** package.
#' The user can select a model to view its details, including its type (if multiple types exist)
#' and the distribution/configuration used.
#'
#' @details
#' This function is designed to facilitate the exploration of Bayesian models within **BayLum**.
#' A selection menu is displayed in the console, allowing the user to choose a model for inspection.
#' Once selected, the model's script is shown in a text window.
#' The user also has the option to save the model script for further analysis.
#'
#' For more details on model usage, see the documentation for datasets.
#'@md


#' @export
extract_Jags_model <- function() {
    #getting names and purpose of the Models
  AllData <- as.data.frame(utils::data(package = "BayLum")$results[, c("Item", "Title")]) %>%
    dplyr::filter(substr(Item, 1, 5) == "Model")

  ## interactive consol

  cli::cat_rule(" Available Models ")
  cli::cat_print(dplyr::tibble(AllData))
  cli::cat_rule()
  cli::cat_line()
  con <- readline(prompt = "Enter the model name : ")

  #If the model does not exist
  if (! (con %in% AllData$Item)) {
    stop(paste(con, "not found in package BayLum"))
  }

  temp_env <- new.env()
    ## Global variable carefull
    data(list = con,  package = "BayLum", envir = temp_env)

    model <- get(con, envir = temp_env)
    cli::cat_rule(paste("Available type in the ", con))
    for (name in names(model)) {
    cli::cat_bullet(name)
    }
  cli::cat_rule()

    type <- readline(prompt = "Enter the model type : ")

    if (class(model[[type]]) == "list") {

      cli::cat_rule(paste(" Available distributions for the", con, type ) )
      for (name in names(model[[type]])) {
        cli::cat_bullet(name)
      }
      cli::cat_rule()

      distribution <- readline(prompt = "Enter the model distribution : ")
      final_model <- model[[type]][[distribution]]

    }
    else {
          final_model <- model[[type]]
    }

    ## printing the file
    temp_file <- tempfile(fileext = ".txt")
    writeLines(final_model, con = temp_file)
    file.show(temp_file)
}



