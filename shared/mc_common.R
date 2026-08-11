################################################################################
### mc_common.R
###
### Shared plumbing for all Monte Carlo scripts:
###   - portable script-location paths
###   - environment-variable parameter parsing
###   - Linux-optimised parallel cluster setup (FORK, single BLAS thread)
###   - safe cv.glmnet wrapper controlled by MC_NFOLDS
###
### This file intentionally does NOT touch any statistical estimator logic.
### It is meant to be source()'d near the top of every simulation script,
### AFTER the BLAS/OpenMP thread environment variables have been set.
################################################################################

################################################################################
### 1. Thread / BLAS oversubscription guard
################################################################################

### Set single-thread BLAS/OpenMP so that the outer foreach parallelism does
### not fight with multi-threaded linear algebra. Scripts should ALSO set these
### at their very top (before loading glmnet), this is a belt-and-suspenders
### copy for anything sourced later.
set_single_thread_blas <- function() {
  Sys.setenv(
    OMP_NUM_THREADS        = "1",
    OPENBLAS_NUM_THREADS   = "1",
    MKL_NUM_THREADS        = "1",
    BLIS_NUM_THREADS       = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )
  ## Env vars are read by OpenBLAS only at library-load time, which for the
  ## conda-forge R (multi-threaded libopenblasp) has usually already happened.
  ## RhpcBLASctl actually re-sets the thread count at runtime -- this is what
  ## prevents each foreach worker from spawning its own BLAS thread pool and
  ## oversubscribing the machine. Call it in the master AND inside each worker.
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
  }
  invisible(TRUE)
}
set_single_thread_blas()

## Lightweight per-worker BLAS cap: call at the top of each foreach iteration.
## FORK workers inherit the master, but calling here is cheap insurance in case
## a worker re-initialised its BLAS thread pool.
worker_single_thread_blas <- function() {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
  }
  invisible(TRUE)
}

################################################################################
### 2. Portable paths
################################################################################

### Directory containing the running Rscript file (falls back to getwd()).
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  normalizePath(getwd())
}

### Candidate directories where shared/util files may live, given that the
### repository is split into Main_simu/ and Comparison/ sub-folders.
util_search_dirs <- function(base_dir) {
  parent <- dirname(base_dir)
  cand <- c(
    base_dir,
    parent,
    file.path(parent, "Main_simu"),
    file.path(parent, "Comparison"),
    getwd(),
    dirname(getwd())
  )
  unique(cand[dir.exists(cand)])
}

### Locate a file by name across the candidate directories. Returns the first
### existing path or stops with a clear message.
locate_file <- function(name, base_dir, must_work = TRUE) {
  dirs <- util_search_dirs(base_dir)
  cand <- file.path(dirs, name)
  hit <- cand[file.exists(cand)]
  if (length(hit) == 0) {
    if (must_work) {
      stop(
        "Could not locate '", name, "'. Searched:\n  ",
        paste(cand, collapse = "\n  ")
      )
    }
    return(NA_character_)
  }
  normalizePath(hit[1])
}

### Source a utility file located anywhere in the search dirs.
source_util <- function(name, base_dir) {
  path <- locate_file(name, base_dir, must_work = TRUE)
  source(path)
  invisible(path)
}

################################################################################
### 3. Environment-variable parsing helpers
################################################################################

parse_int_vec  <- function(x) as.integer(strsplit(x, ",")[[1]])
parse_char_vec <- function(x) trimws(strsplit(x, ",")[[1]])

env_chr <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") default else v
}

env_int <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") as.integer(default) else as.integer(v)
}

env_num <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") as.numeric(default) else as.numeric(v)
}

env_int_vec <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") parse_int_vec(default) else parse_int_vec(v)
}

env_lgl <- function(name, default) {
  v <- Sys.getenv(name, unset = NA_character_)
  if (is.na(v) || v == "") return(isTRUE(default))
  tolower(v) %in% c("1", "true", "t", "yes", "y")
}

################################################################################
### 4. Package handling (no auto-install)
################################################################################

### Verify that every required package is installed. Does NOT install anything.
### Stops with an actionable message pointing at setup_packages.R.
require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(
    pkgs,
    function(p) requireNamespace(p, quietly = TRUE),
    logical(1)
  )]
  if (length(missing) > 0) {
    stop(
      "Missing required R packages: ", paste(missing, collapse = ", "), "\n",
      "Install them once with:\n",
      "    Rscript setup_packages.R\n",
      "or in R:\n",
      "    install.packages(c(", paste0("\"", missing, "\"", collapse = ", "), "))",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

################################################################################
### 5. Parallel cluster setup (Linux optimised)
################################################################################

### Build and register a parallel cluster. FORK on unix (cheap, shares memory),
### PSOCK elsewhere. Returns the cluster object; caller is responsible for
### stopCluster (typically via on.exit).
make_mc_cluster <- function(iter, seed) {
  n.cores <- as.integer(Sys.getenv(
    "MC_CORES",
    unset = max(1L, parallel::detectCores(logical = TRUE) - 4L)
  ))
  n.cores <- max(1L, min(n.cores, iter))

  cluster_type <- if (.Platform$OS.type == "unix") "FORK" else "PSOCK"

  cat("Detected cores:", parallel::detectCores(logical = TRUE),
      "| Using cores:", n.cores,
      "| Cluster type:", cluster_type, "\n")

  cl <- parallel::makeCluster(n.cores, type = cluster_type)
  doParallel::registerDoParallel(cl = cl)
  doRNG::registerDoRNG(seed)
  cl
}

### foreach scheduling options: disable prescheduling so that replications with
### uneven runtimes get load-balanced. multicore option covers FORK, snow covers
### PSOCK; foreach ignores whichever does not apply.
mc_foreach_options_multicore <- list(preschedule = FALSE)
mc_foreach_options_snow      <- list(preschedule = FALSE)

################################################################################
### 6. Safe cv.glmnet
################################################################################

### Choose a valid number of CV folds: at least 3 (glmnet minimum), at most n,
### capped by MC_NFOLDS.
cv_nfolds <- function(n, nfolds = getOption("MC_NFOLDS", 3L)) {
  n <- as.integer(n)
  nf <- as.integer(nfolds)
  max(3L, min(nf, n))
}

### Estimate a B0 coefficient matrix (regressing each column of M on X) using a
### single lasso path per column with the lambda chosen by BIC -- NO K-fold
### cross-validation. This is the "C3" estimator: ~5-6x faster than per-column
### cv.glmnet at high dimension because it fits one glmnet path (not nfolds+1)
### per response.
###
### Returns a matrix with one row per column of M (ncol(M) x ncol(X)), matching
### the shape produced by the original per-column cv.glmnet loop
### (B0_hat[j, ] = coefficients for M[, j]).
###
### On any per-column failure that column's row is left at 0 and a flag is set
### via attr(out, "n_failed_cols"); callers may inspect it but need not.
b0_bic_matrix <- function(Xsub, Msub, nlambda = 50L) {
  Xsub <- as.matrix(Xsub)
  Msub <- as.matrix(Msub)
  n <- nrow(Xsub)
  pX <- ncol(Xsub)
  qM <- ncol(Msub)

  B <- matrix(0, nrow = qM, ncol = pX)
  n_failed <- 0L

  for (j in seq_len(qM)) {
    y <- Msub[, j]
    ok <- tryCatch({
      g <- glmnet::glmnet(Xsub, y, alpha = 1, intercept = FALSE, nlambda = nlambda)
      pred <- predict(g, Xsub)
      rss  <- colSums((y - pred)^2)
      ## BIC = n * log(RSS/n) + df * log(n); guard RSS = 0.
      bic  <- n * log(pmax(rss, .Machine$double.eps) / n) + g$df * log(n)
      k    <- which.min(bic)
      B[j, ] <- as.vector(g$beta[, k])
      TRUE
    }, error = function(e) FALSE)
    if (!ok) n_failed <- n_failed + 1L
  }

  attr(B, "n_failed_cols") <- n_failed
  B
}

### Safe cv.glmnet: clamps nfolds, guards against degenerate inputs, and returns
### NULL on failure instead of throwing. Callers must handle NULL.
safe_cv_glmnet <- function(x, y, nfolds = getOption("MC_NFOLDS", 3L), ...) {
  x <- as.matrix(x)
  n <- nrow(x)

  if (n < 3L) return(NULL)
  if (all(!is.finite(y))) return(NULL)
  if (stats::sd(as.numeric(y), na.rm = TRUE) == 0) {
    ## glmnet cannot cross-validate a constant response
    return(NULL)
  }

  nf <- cv_nfolds(n, nfolds)

  tryCatch(
    glmnet::cv.glmnet(x, y, nfolds = nf, ...),
    error = function(e) NULL
  )
}

################################################################################
### 7. Tags and output helpers
################################################################################

### e.g. sigma = 1 -> "sigma_1", sigma = 0.5 -> "sigma_0p5"
make_sigma_tag <- function(sigma) {
  paste0("sigma_", gsub("\\.", "p", format(sigma, trim = TRUE, scientific = FALSE)))
}

### Build an all-NA named row matching a template of column names, marking it as
### a failed replication. Keeps rbind shapes consistent across replications.
na_result_row <- function(cols) {
  v <- rep(NA_real_, length(cols))
  names(v) <- cols
  if ("failed" %in% cols) v["failed"] <- 1
  v
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}
