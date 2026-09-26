#!/usr/bin/env Rscript

# Step 0: fetch the raw dataset.
# This script downloads and unpacks the UCI dataset into data/
# so the rest of the pipeline can run from a fresh clone.

library(utils)

ROOT <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."))
DATA <- file.path(ROOT, "data")
URL <- paste0(
  "https://archive.ics.uci.edu/static/public/235/",
  "individual+household+electric+power+consumption.zip"
)
ZIP <- file.path(DATA, "household_power_consumption.zip")
TXT <- file.path(DATA, "household_power_consumption.txt")

EXPECTED_ROWS <- 2075260   # includes header
EXPECTED_BYTES <- 132960755

# Helper: print progress
download_with_progress <- function(url, dest) {
  cat("downloading", url, "\n")
  tryCatch({
    utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE)
  }, error = function(e) {
    if (file.exists(dest)) file.remove(dest)
    stop(
      "download failed: ", e$message,
      "\nFetch it manually from https://archive.ics.uci.edu/dataset/235/ into data/"
    )
  })
}

main <- function() {
  dir.create(DATA, showWarnings = FALSE, recursive = TRUE)
  
  # If TXT already exists and matches expected size, skip
  if (file.exists(TXT) && file.info(TXT)$size == EXPECTED_BYTES) {
    cat("already present and correct size:", TXT, "\n")
    return(invisible(NULL))
  }
  
  # Download ZIP if missing
  if (!file.exists(ZIP)) {
    download_with_progress(URL, ZIP)
  }
  
  # Unpack ZIP
  cat("unpacking", basename(ZIP), "\n")
  unzip(ZIP, exdir = DATA)
  
  if (!file.exists(TXT)) {
    stop("expected ", TXT, " inside the archive but it is not there")
  }
  
  # Validate size and row count
  size <- file.info(TXT)$size
  n_rows <- R.utils::countLines(TXT)
  
  cat("  size", format(size, big.mark = ","), "bytes (expected",
      format(EXPECTED_BYTES, big.mark = ","), ")\n")
  cat("  rows", format(n_rows, big.mark = ","), "(expected",
      format(EXPECTED_ROWS, big.mark = ","), ")\n")
  
  if (size != EXPECTED_BYTES || n_rows != EXPECTED_ROWS) {
    stop("raw file does not match the documented dataset; stopping so that ",
         "downstream results are not silently wrong")
  }
  
  cat("\nok. Next: Rscript src/01_aggregate.R\n")
}

main()

