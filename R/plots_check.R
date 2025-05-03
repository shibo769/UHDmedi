 
#' Plot Bootstrap Samples: Histograms and QQ Plots
#'
#' This function generates diagnostic plots for bootstrap estimates.
#' It creates histograms (with vertical red dashed lines marking the 2.5\% and 97.5\% quantiles)
#' and QQ plots (comparing the sample quantiles with a normal distribution) for each effect.
#'
#' @param boot_samples A data.frame containing bootstrap estimates. It should have columns named
#'   TE, NIE, NDE, Y11, Y10, and Y00.
#'
#' @return A list with two ggplot objects:
#'   \describe{
#'     \item{histogram}{Histogram plots for each effect with 95\% CI indicated.}
#'     \item{qqplot}{QQ plots for each effect.}
#'   }
#'
#' @examples
#' \dontrun{
#'   # Suppose boot_results$boot_samples is the bootstrap samples data.frame from bootstrap_est_effects()
#'   plots <- plots_check(boot_results$boot_samples)
#'   print(plots$histogram)
#'   print(plots$qqplot)
#' }
#'
#' @export
plots_check <- function(boot_samples) {
  if (!requireNamespace('ggplot2', quietly = TRUE)) {
    stop('Package "ggplot2" is required for plots_check(). Please install it.')
  }
  if (!requireNamespace('reshape2', quietly = TRUE)) {
    stop('Package "reshape2" is required for plots_check(). Please install it.')
  }
  library(ggplot2)
  library(reshape2)
  
  # Convert bootstrap samples to long format
  boot_long <- reshape2::melt(boot_samples, variable.name = 'Effect', value.name = 'Estimate')
  
  # Compute CI for each effect
  ci_df <- aggregate(Estimate ~ Effect, data = boot_long, FUN = function(x) {
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  })
  ci_df <- do.call(data.frame, ci_df)
  names(ci_df)[2:3] <- c('Lower', 'Upper')
  
  # Histogram with vertical lines for CI
  hist_plot <- ggplot2::ggplot(boot_long, ggplot2::aes(x = Estimate)) +
    ggplot2::geom_histogram(bins = 30, fill = 'skyblue', color = 'black') +
    ggplot2::facet_wrap(~ Effect, scales = 'free') +
    ggplot2::geom_vline(data = ci_df, ggplot2::aes(xintercept = Lower), color = 'red', linetype = 'dashed', size = 1) +
    ggplot2::geom_vline(data = ci_df, ggplot2::aes(xintercept = Upper), color = 'red', linetype = 'dashed', size = 1) +
    ggplot2::labs(title = 'Bootstrap Estimates Histogram with 95% CI', x = 'Estimate', y = 'Frequency') +
    ggplot2::theme_minimal()
  
  # QQ plot for each effect
  qq_plot <- ggplot2::ggplot(boot_long, ggplot2::aes(sample = Estimate)) +
    ggplot2::stat_qq(color = 'blue') +
    ggplot2::stat_qq_line(color = 'red') +
    ggplot2::facet_wrap(~ Effect, scales = 'free') +
    ggplot2::labs(title = 'QQ Plots of Bootstrap Estimates', x = 'Theoretical Quantiles', y = 'Sample Quantiles') +
    ggplot2::theme_minimal()
  
  return(list(histogram = hist_plot, qqplot = qq_plot))
}

