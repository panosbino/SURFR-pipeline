# =============================================================================
# SURFR Pipeline — install_r_packages.R
#
# Installs the R packages used by FindCancerSpecificRNAs.r from a date-pinned
# CRAN snapshot (env var CRAN, set in the Dockerfile), so every rebuild of the
# image gets the same package versions. Fails the build on any problem and
# records the installed versions in /opt/surfr/R_package_versions.tsv.
# =============================================================================

expected_r <- "4.4.1"
pkgs <- c("tidyverse", "paletteer", "arrow", "ggvenn", "MASS")

# --- Guard: correct R version -------------------------------------------------
if (getRversion() != expected_r) {
  stop("Expected R ", expected_r, ", found ", getRversion())
}

# --- Guard: dated snapshot, never a moving repository -------------------------
repo <- Sys.getenv("CRAN")
if (!grepl("/\\d{4}-\\d{2}-\\d{2}$", repo)) {
  stop("CRAN must point to a dated snapshot (…/YYYY-MM-DD), got: '", repo, "'")
}
options(repos = c(CRAN = repo), Ncpus = parallel::detectCores())

# --- Install ------------------------------------------------------------------
# MASS is a recommended package shipped with R 4.4.1; setdiff() keeps that
# bundled version instead of replacing it with the snapshot's.
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install)

# install.packages() only warns on failure, so verify explicitly.
ok <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
if (!all(ok)) {
  stop("Failed to install: ", paste(pkgs[!ok], collapse = ", "))
}

# --- Record exact versions ----------------------------------------------------
ip <- installed.packages()[, c("Package", "Version", "LibPath"), drop = FALSE]
write.table(ip, "/opt/surfr/R_package_versions.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("R", as.character(getRversion()), "| snapshot:", repo, "\n")
for (p in pkgs) cat(sprintf("  %-10s %s\n", p, as.character(packageVersion(p))))
