 
#' Load Real Lung Cancer Data
#'
#' This function loads the real lung cancer data from the package's 'lung_cancer_data' folder.
#' It returns a list containing the four datasets:
#'   - M: Mediator data (from 'M_0.05_theshold_1974.rds')
#'   - X: Predictor data (from 'X_thirdorder.rds')
#'   - A: Treatment (Smoking) indicator (from 'A_Smoking.rds')
#'   - Y: Outcome (Survival Time) data (from 'Y_SurvivalTime.rds')
#'
#' @return A list with elements M, X, A, Y.
#'
#' @examples
#' \dontrun{
#'   data_list <- real_data()
#'   attach(data_list)  # Now you can directly use M, X, A, Y
#' }
#'
#' @export
real_data <- function() {
  data_path <- system.file('lung_cancer_data', package = 'UHDmedi')
  if(data_path == '') {
    stop('Lung cancer data folder not found in package. Please ensure the folder and data files exist.')
  }
  M <- readRDS(file.path(data_path, 'M_0.05_theshold_1974.rds'))
  X <- readRDS(file.path(data_path, 'X_thirdorder.rds'))
  A <- readRDS(file.path(data_path, 'A_Smoking.rds'))
  Y <- readRDS(file.path(data_path, 'Y_SurvivalTime.rds'))
  
  data_list <- list(M = M, X = X, A = A, Y = Y)
  return(data_list)
}

