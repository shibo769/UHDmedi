################################################################################
### Verify the original (pre-cleaning) mediator matrix against the cleaned one.
###
### The cleaned matrix keeps its ORIGINAL column indices in its colnames
### ("V1", "V2", "V5", ...), with 85,222 indices missing out of 350,019. If the
### cleaning rule was "drop any column containing an NA", then the set of
### columns of M_raw that contain an NA must coincide exactly with those missing
### indices. This script checks that, and reports the alternative rules if not.
################################################################################

data_dir <- Sys.getenv("DATA_DIR", unset = "real_data")
raw_path <- file.path(data_dir, "M_raw.rds")
if (!file.exists(raw_path)) stop("M_raw.rds not found in ", data_dir)

cat("Loading M_raw.rds ...\n")
t0 <- Sys.time()
M_raw <- readRDS(raw_path)
cat(sprintf("loaded in %.0f s | class = %s | dim = %s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs")),
            paste(class(M_raw), collapse = "/"),
            paste(dim(M_raw), collapse = " x ")))

M_cl <- readRDS(file.path(data_dir, "M_raw_cleaned.rds"))
kept_idx <- as.integer(sub("^V", "", colnames(M_cl)))   # original positions retained
cat(sprintf("cleaned: dim = %s | original indices span %d..%d | %d missing\n",
            paste(dim(M_cl), collapse = " x "), min(kept_idx), max(kept_idx),
            max(kept_idx) - length(kept_idx)))

## ---- dimensions -------------------------------------------------------------
cat("\n== dimensions ==\n")
cat(sprintf("rows: raw %d vs cleaned %d %s\n", nrow(M_raw), nrow(M_cl),
            if (nrow(M_raw) == nrow(M_cl)) "(match: no samples dropped)" else "(DIFFER: rows were filtered too)"))
cat(sprintf("cols: raw %d vs max kept index %d %s\n", ncol(M_raw), max(kept_idx),
            if (ncol(M_raw) == max(kept_idx)) "(match)" else "(MISMATCH)"))

## ---- NA structure -----------------------------------------------------------
cat("\n== NA structure of the raw matrix ==\n")
na_per_col <- colSums(is.na(M_raw))
has_na  <- which(na_per_col > 0)
cat(sprintf("total NAs: %s | columns with >=1 NA: %d | columns fully NA: %d\n",
            format(sum(na_per_col), big.mark = ","), length(has_na),
            sum(na_per_col == nrow(M_raw))))
if (length(has_na)) {
  cat("NA count per affected column -- quantiles:\n")
  print(quantile(na_per_col[has_na], c(0, .25, .5, .75, .9, 1)))
}

## ---- does "drop any column with an NA" reproduce the cleaned matrix? --------
cat("\n== testing the cleaning rule ==\n")
dropped_idx <- setdiff(seq_len(ncol(M_raw)), kept_idx)
cat(sprintf("columns dropped by the cleaning: %d\n", length(dropped_idx)))
cat(sprintf("columns of raw containing an NA: %d\n", length(has_na)))
if (identical(sort(dropped_idx), sort(as.integer(has_na)))) {
  cat(">>> EXACT MATCH: the rule is complete-case-by-column (drop a CpG if it has any NA).\n")
} else {
  cat(">>> NOT an exact match. Breaking it down:\n")
  cat(sprintf("    dropped but NA-free      : %d\n", length(setdiff(dropped_idx, has_na))))
  cat(sprintf("    has NA but was kept      : %d\n", length(setdiff(has_na, dropped_idx))))
  ## a plausible alternative: an NA-fraction threshold
  fr <- na_per_col / nrow(M_raw)
  for (thr in c(0, 0.01, 0.05, 0.1, 0.2, 0.5)) {
    pred <- which(fr > thr)
    cat(sprintf("    threshold NA-frac > %-4g -> %d dropped (%s)\n", thr, length(pred),
                if (identical(sort(as.integer(pred)), sort(dropped_idx))) "EXACT MATCH" else "no"))
  }
}

## ---- do the retained values actually agree? --------------------------------
cat("\n== value agreement on retained columns (random sample) ==\n")
set.seed(1)
chk <- sort(sample(seq_along(kept_idx), min(200, length(kept_idx))))
ok <- vapply(chk, function(j)
  isTRUE(all.equal(as.numeric(M_raw[, kept_idx[j]]), as.numeric(M_cl[, j]),
                   tolerance = 1e-12)), logical(1))
cat(sprintf("%d / %d sampled columns identical\n", sum(ok), length(ok)))
if (any(!ok)) cat("first mismatching cleaned-column indices: ",
                  paste(head(chk[!ok], 10), collapse = ", "), "\n")

## ---- annotation alignment ---------------------------------------------------
g_path <- file.path(data_dir, "gene_sites_for_M_raw_cleaned.rds")
if (file.exists(g_path)) {
  g <- readRDS(g_path)
  cat(sprintf("\n== annotation ==\nlength(gene_sites) = %d | ncol(cleaned) = %d %s\n",
              length(g), ncol(M_cl), if (length(g) == ncol(M_cl)) "(aligned)" else "(MISALIGNED)"))
  if (!is.null(colnames(M_raw)) && all(grepl("^cg|^ch\\.|^rs", head(colnames(M_raw), 20)))) {
    cat("raw matrix carries CpG ids in its colnames -- checking against the annotation:\n")
    cat(sprintf("  identical on retained positions: %s\n",
                identical(as.character(colnames(M_raw)[kept_idx]), as.character(g))))
  } else {
    cat("raw colnames (first 5): ", paste(head(colnames(M_raw), 5), collapse = ", "), "\n")
  }
}
cat("\ndone.\n")
