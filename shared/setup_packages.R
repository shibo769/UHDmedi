################################################################################
### setup_packages.R
###
### One-time installation of every R package the simulation project needs.
### Run this ONCE (not inside the simulation scripts):
###
###     Rscript setup_packages.R
###
### The simulation scripts themselves only CHECK for packages and stop with a
### clear message if something is missing; they never auto-install.
################################################################################

## Pick a CRAN mirror if none configured.
repos <- getOption("repos")
if (is.null(repos) || is.na(repos["CRAN"]) || repos["CRAN"] == "@CRAN@") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

## Use all but a few cores to compile faster.
ncpu <- max(1L, parallel::detectCores(logical = TRUE) - 2L)
options(Ncpus = ncpu)

################################################################################
### CRAN packages
################################################################################

cran_pkgs <- c(
  "foreach",
  "doParallel",
  "doRNG",
  "pracma",
  "MASS",
  "glmnet",
  "quadprog",
  "ggplot2",
  "gridExtra",
  "dplyr",
  "readr",
  "stringr",
  "tidyr",
  "purrr",
  "ncvreg",
  "remotes"    # needed to install pogs from GitHub
)

installed <- rownames(installed.packages())
to_install <- setdiff(cran_pkgs, installed)

if (length(to_install) > 0) {
  cat("Installing CRAN packages:\n  ", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install)
} else {
  cat("All CRAN packages already installed.\n")
}

################################################################################
### pogs (non-CRAN)
###
### The residual-balancing debiased estimator uses optimizer = "pogs.new",
### which calls the pogs ADMM solver (foges/pogs). Without pogs, approx.balance
### silently falls back to a DIFFERENT optimizer (quadprog), which changes the
### estimator. We therefore install pogs from GitHub here.
################################################################################

if (!requireNamespace("pogs", quietly = TRUE)) {

  ## Preferred: install from the vendored source tarball checked into this repo
  ## (vendor/pogs_1.0.tar.gz). This is the foges/pogs R interface (commit
  ## 6e172d7, the last self-consistent CPU R binding) with an idempotent
  ## configure, so it builds offline against a modern toolchain.
  script_dir <- tryCatch({
    args <- commandArgs(trailingOnly = FALSE)
    fa <- grep("^--file=", args, value = TRUE)
    if (length(fa) > 0) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
  }, error = function(e) getwd())

  vendored <- file.path(script_dir, "vendor", "pogs_1.0.tar.gz")

  ok <- FALSE
  if (file.exists(vendored)) {
    cat("Installing pogs from vendored tarball:", vendored, "\n")
    ok <- tryCatch({
      install.packages(vendored, repos = NULL, type = "source")
      requireNamespace("pogs", quietly = TRUE)
    }, error = function(e) {
      message("pogs vendored install failed: ", conditionMessage(e)); FALSE
    })
  }

  ## Fallback: try GitHub (only works with network access to api.github.com).
  if (!ok) {
    cat("Vendored install unavailable; trying GitHub (foges/pogs)...\n")
    ok <- tryCatch({
      remotes::install_github("foges/pogs/src/interface_r", upgrade = "never")
      requireNamespace("pogs", quietly = TRUE)
    }, error = function(e) {
      message("pogs GitHub install failed: ", conditionMessage(e)); FALSE
    })
  }

  if (!ok || !requireNamespace("pogs", quietly = TRUE)) {
    cat("\n*** WARNING: pogs is NOT installed. ***\n",
        "The debiased estimator (optimizer = 'pogs.new') requires it.\n",
        "Without pogs the code falls back to quadprog, a DIFFERENT estimator.\n",
        "Set MC_ALLOW_QUADPROG_FALLBACK=TRUE only for quick smoke tests.\n",
        sep = "")
  }
} else {
  cat("pogs already installed.\n")
}

################################################################################
### Final report
################################################################################

all_pkgs <- c(cran_pkgs, "pogs")
final_installed <- rownames(installed.packages())
still_missing <- setdiff(all_pkgs, final_installed)

if (length(still_missing) == 0) {
  cat("\nAll packages installed successfully.\n")
} else {
  cat("\nStill missing after setup:\n  ", paste(still_missing, collapse = ", "), "\n")
}
