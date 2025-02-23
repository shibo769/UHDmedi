 
#' Exploratory Data Analysis Report
#'
#' This function performs exploratory data analysis (EDA) for the input data.
#' It reports basic information such as dimensions, missing values, and summary statistics for predictors (X) and mediators (M).
#' If outcome (Y) is provided, density and box plots are generated, along with outlier information.
#' If treatment assignment (A) is provided, its distribution is also reported.
#'
#' @param X A data.frame or matrix containing predictors.
#' @param M A data.frame or matrix containing mediators.
#' @param A (optional) A vector of treatment assignments.
#' @param Y (optional) A vector of outcomes.
#' @param plot_corr Logical. If TRUE, a correlation matrix of X is plotted. Default is TRUE.
#' @param plot_Y Logical. If TRUE, density and box plots for Y are generated (if Y is provided). Default is TRUE.
#' @param verbose Logical. If TRUE, summary information is printed to the console. Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{summary}{A list with dimensions, missing values, and summary statistics for X, M, (and Y, A if provided).}
#'   \item{plots}{A list of ggplot objects, including the correlation matrix plot for X and density and box plots for Y (if provided).}
#' }
#'
#' @examples
#' \dontrun{
#' eda_results <- eda_report(X = sim_data$X,
#'                           M = sim_data$M,
#'                           A = sim_data$A,
#'                           Y = sim_data$Y,
#'                           plot_corr = TRUE,
#'                           plot_Y = TRUE,
#'                           verbose = TRUE)
#' 
#' # View EDA output
#' print(eda_results$summary)
#' # To view plots in an interactive environment, run:
#' print(eda_results$plots$scatter_plot_minvar)
#' print(eda_results$plots$scatter_plot_maxvar)
#' print(eda_results$plots$corr_plot)
#' print(eda_results$plots$Y_density)
#' print(eda_results$plots$Y_boxplot)
#' }
#'
#' @export
eda_report <- function(X, M, A = NULL, Y = NULL, plot_corr = TRUE, plot_Y = TRUE, verbose = TRUE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" is required for eda_report(). Please install it.")
  }
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Package \"reshape2\" is required for eda_report(). Please install it.")
  }
  library(ggplot2)
  library(reshape2)
  
  summary_list <- list()
  
  if (!is.null(Y)) {
    summary_list$Y_summary <- summary(Y)
  }
  
  if (!is.null(A)) {
    summary_list$A_summary <- summary(A)
  }
  
  # Convert X and M to data.frame if they are matrices
  if (is.matrix(X)) {
    X <- as.data.frame(X)
  }
  if (is.matrix(M)) {
    M <- as.data.frame(M)
  }
  
  # Dimensions
  summary_list$X_dim <- dim(X)
  summary_list$M_dim <- dim(M)
  if (!is.null(Y)) {
    summary_list$Y_length <- length(Y)
  }
  if (!is.null(A)) {
    summary_list$A_length <- length(A)
  }
  
  # Missing values for X
  missing_X <- colSums(is.na(X))
  has_na_X <- which(missing_X > 0)  # search for columns with missing values
  if (length(has_na_X) > 0) {
    # keep the columns with missing values and store the number of columns
    missing_X_df <- data.frame(Column = names(missing_X)[has_na_X],
                               MissingCount = missing_X[has_na_X])
  } else {
    missing_X_df <- data.frame(Column = character(0), MissingCount = numeric(0))
  }
  
  # Missing values for M
  missing_M <- colSums(is.na(M))
  has_na_M <- which(missing_M > 0)
  if (length(has_na_M) > 0) {
    missing_M_df <- data.frame(Column = names(missing_M)[has_na_M],
                               MissingCount = missing_M[has_na_M])
  } else {
    missing_M_df <- data.frame(Column = character(0), MissingCount = numeric(0))
  }
  
  summary_list$missing_X <- missing_X_df
  summary_list$missing_M <- missing_M_df
  
  # For Y and A, present the number of missing values
  if (!is.null(Y)) {
    y_na_count <- sum(is.na(Y))
    summary_list$missing_Y <- y_na_count
  }
  if (!is.null(A)) {
    a_na_count <- sum(is.na(A))
    summary_list$missing_A <- a_na_count
  }
  
  # categorical or numerical?
  is_categorical_column <- function(x, threshold = 5) {
    # if factor or character，then categorical
    if (is.factor(x) || is.character(x)) {
      return(TRUE)
    }
    # numerical
    if (is.numeric(x) || is.integer(x)) {
      uniq_vals <- unique(x[!is.na(x)])
      if (length(uniq_vals) <= threshold) {
        return(TRUE)
      }
    }
    return(FALSE)
  }
  
  # Convert X to data.frame if it's a matrix
  if (is.matrix(X)) {
    X <- as.data.frame(X)
  }
  
  # summarize categorical and continuous
  cat_cols_indices <- sapply(X, is_categorical_column, threshold = 5)
  cat_cols <- sum(cat_cols_indices)
  num_cols <- ncol(X) - cat_cols
  
  summary_list$X_cat_columns <- cat_cols
  summary_list$X_num_columns <- num_cols
  
  # Range for X and M
  summary_list$X_range <- list(min = min(X, na.rm = TRUE), max = max(X, na.rm = TRUE))
  summary_list$M_range <- list(min = min(M, na.rm = TRUE), max = max(M, na.rm = TRUE))
  
  # Correlation for first 10 columns of M and X with Y
  if(!is.null(Y)) {
    # For M: calculate marginal correlation between M and Y and pick the first ten
    corr_M <- sapply(1:ncol(M), function(i) cor(M[, i], Y, use = "complete.obs"))
    indices_M <- order(abs(corr_M), decreasing = TRUE)[1:min(10, length(corr_M))]
    top10_corr_M <- corr_M[indices_M]
    summary_list$cor_M_Y <- top10_corr_M
    
    # For X: calculate marginal correlation between X and Y and pick the first ten
    corr_X <- sapply(1:ncol(X), function(i) cor(X[, i], Y, use = "complete.obs"))
    indices_X <- order(abs(corr_X), decreasing = TRUE)[1:min(10, length(corr_X))]
    top10_corr_X <- corr_X[indices_X]
    summary_list$cor_X_Y <- top10_corr_X
  }
  
  # Minimum and maximum variance in M
  var_M <- apply(M, 2, var, na.rm = TRUE)
  min_var <- min(var_M, na.rm = TRUE)
  min_var_col <- which(var_M == min_var)[1]
  summary_list$min_var_M <- list(column = min_var_col, variance = min_var)
  
  max_var <- max(var_M, na.rm = TRUE)
  max_var_col <- which(var_M == max_var)[1]
  summary_list$max_var_M <- list(column = max_var_col, variance = max_var)
  
  plots <- list()
  
  # Scatterplot for the M column with minimum variance versus Y
  if(!is.null(Y)) {
    scatter_plot_minvar <- ggplot2::ggplot(data = data.frame(M_col = M[, min_var_col], Y = Y), 
                                    ggplot2::aes(x = M_col, y = Y)) +
      ggplot2::geom_point(color = "darkblue", alpha = 0.7) +
      ggplot2::geom_smooth(method = "lm", col = "red") +
      ggplot2::labs(title = sprintf("Scatterplot of M[, %d] (min variance) vs Y", min_var_col),
                    x = sprintf("M[, %d]", min_var_col),
                    y = "Y") +
      ggplot2::theme_minimal()
    plots$scatter_plot_minvar <- scatter_plot_minvar
  }
  
  if(!is.null(Y)) {
    scatter_plot_maxvar <- ggplot2::ggplot(data = data.frame(M_col = M[, max_var_col], Y = Y), 
                                    ggplot2::aes(x = M_col, y = Y)) +
      ggplot2::geom_point(color = "darkblue", alpha = 0.7) +
      ggplot2::geom_smooth(method = "lm", col = "red") +
      ggplot2::labs(title = sprintf("Scatterplot of M[, %d] (max variance) vs Y", max_var_col),
                    x = sprintf("M[, %d]", min_var_col),
                    y = "Y") +
      ggplot2::theme_minimal()
    plots$scatter_plot_maxvar <- scatter_plot_maxvar
  }
  
  
  # Correlation matrix for X (only numeric columns)
  if (plot_corr) {
    if (is.data.frame(X)) {
      X_num <- X[, sapply(X, is.numeric), drop = FALSE]
    } else {
      X_num <- X
    }
    if (ncol(X_num) == 0) {
      stop("No numeric columns found in X for correlation plot.")
    }
    corr_mat <- cor(X_num, use = "pairwise.complete.obs")
    corr_long <- melt(corr_mat)
    plots$corr_plot <- ggplot(corr_long, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Correlation Matrix of X", x = "", y = "", fill = "Correlation")
  }
  
  # Density and boxplot for Y
  if (plot_Y && !is.null(Y)) {
    Y_df <- data.frame(Y = Y)
    plots$Y_density <- ggplot(Y_df, aes(x = Y)) +
      geom_density(fill = "skyblue", alpha = 0.5) +
      labs(title = "Density Plot of Y", x = "Y", y = "Density") +
      theme_minimal()
    
    # plot provided by A
    if (!is.null(A)) {
      Y_df$A <- as.factor(A)
      plots$Y_boxplot <- ggplot(Y_df, aes(x = A, y = Y, fill = A)) +
        geom_boxplot() +
        labs(title = "Boxplot of Y by Treatment", x = "Treatment (A)", y = "Y") +
        coord_flip() +
        theme_minimal() +
        scale_fill_manual(values = c("0" = "lightgreen", "1" = "skyblue"))
    } else {
      plots$Y_boxplot <- ggplot(Y_df, aes(y = Y)) +
        geom_boxplot(fill = "lightgreen") +
        labs(title = "Boxplot of Y", y = "Y") +
        theme_minimal()
    }
  }
  
  
  if (verbose) {
    cat("\n==========================================================\n")
    cat("                  EDA REPORT SUMMARY\n")
    cat("==========================================================\n\n")
    
    # ----------------- Dimensions ------------------
    cat("----------------------------------------------------------\n")
    cat("                     DIMENSIONS\n")
    cat("----------------------------------------------------------\n")
    cat("X dimension: ")
    print(summary_list$X_dim)
    cat("M dimension: ")
    print(summary_list$M_dim)
    if (!is.null(Y)) {
      cat("Length of Y:", summary_list$Y_length, "\n")
    }
    if (!is.null(A)) {
      cat("Length of A:", summary_list$A_length, "\n")
    }
    
    # ----------------- Missing Values ------------------
    cat("\n----------------------------------------------------------\n")
    cat("                    MISSING VALUES\n")
    cat("----------------------------------------------------------\n")
    
    cat("Missing values in X:\n")
    if (nrow(summary_list$missing_X) == 0) {
      cat("  No missing values in X.\n")
    } else {
      print(summary_list$missing_X)
    }
    
    cat("\nMissing values in M:\n")
    if (nrow(summary_list$missing_M) == 0) {
      cat("  No missing values in M.\n")
    } else {
      print(summary_list$missing_M)
    }
    
    if (!is.null(Y)) {
      cat("\nMissing values in Y:", summary_list$missing_Y, "\n")
    }
    if (!is.null(A)) {
      cat("Missing values in A:", summary_list$missing_A, "\n")
    }
    
    # ----------------- Y Summary (if provided) ------------------
    if (!is.null(Y)) {
      cat("\n----------------------------------------------------------\n")
      cat("                      SUMMARY OF Y\n")
      cat("----------------------------------------------------------\n")
      print(summary_list$Y_summary)
    }
    
    # ----------------- A Summary (if provided) ------------------
    if (!is.null(A)) {
      cat("\n----------------------------------------------------------\n")
      cat("                    DISTRIBUTION OF A\n")
      cat("----------------------------------------------------------\n")
      print(summary_list$A_summary)
    }
    
    # ----------------- X Characteristics ------------------
    cat("\n----------------------------------------------------------\n")
    cat("                 CHARACTERISTICS FOR X\n")
    cat("----------------------------------------------------------\n")
    cat("Number of numeric columns in X:", summary_list$X_num_columns, "\n")
    cat("Number of categorical columns in X:", summary_list$X_cat_columns, "\n")
    cat("Range of X: [", summary_list$X_range$min, ", ", summary_list$X_range$max, "]\n")
    
    # ----------------- M Characteristics ------------------
    cat("\n----------------------------------------------------------\n")
    cat("                 CHARACTERISTICS FOR M\n")
    cat("----------------------------------------------------------\n")
    cat("Range of M: [", summary_list$M_range$min, ", ", summary_list$M_range$max, "]\n")
    
    # ----------------- Correlations with Y ------------------
    if (!is.null(Y)) {
      cat("\n----------------------------------------------------------\n")
      cat("     TOP 10 CORRELATIONS (BY ABS. VALUE) WITH Y\n")
      cat("----------------------------------------------------------\n")
      cat("M with Y:\n")
      print(summary_list$cor_M_Y)
      
      cat("\nX with Y:\n")
      print(summary_list$cor_X_Y)
    }
    
    # ----------------- Variance in M ------------------
    cat("\n----------------------------------------------------------\n")
    cat("               VARIANCE INFORMATION FOR M\n")
    cat("----------------------------------------------------------\n")
    cat("Minimum variance in M:\n")
    cat("  Column:", summary_list$min_var_M$column, 
        "with variance:", summary_list$min_var_M$variance, "\n")
    cat("Maximum variance in M:\n")
    cat("  Column:", summary_list$max_var_M$column, 
        "with variance:", summary_list$max_var_M$variance, "\n")
    
    cat("\nScatterplot for the M column with minimum variance vs Y is stored in plots$scatter_plot.\n")
    
    cat("\n==========================================================\n")
    cat("                  END OF EDA REPORT\n")
    cat("==========================================================\n\n")
  }
  
  
  
  return(list(summary = summary_list, plots = plots))
}
