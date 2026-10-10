#!/usr/bin/env Rscript
# Step 0: fetch the raw dataset.
#
# The raw file is not tracked in git: it is 127 MB, above GitHub's 100 MB
# per-file limit. This script downloads and unpacks it into data/ so the rest of
# the pipeline can run from a fresh clone.
#
#   Rscript src/00_download.R

script_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 0) return(normalizePath("."))
  script <- sub("^--file=", "", file_arg[[1]])
  normalizePath(file.path(dirname(normalizePath(script)), ".."))
}

count_lines <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  n <- 0
  last_nl <- TRUE
  repeat {
    buf <- readBin(con, "raw", 1000000L)
    if (length(buf) == 0) break
    n <- n + sum(buf == as.raw(10L))
    last_nl <- buf[length(buf)] == as.raw(10L)
  }
  if (!last_nl && n > 0) n <- n + 1
  n
}

main <- function() {
  ROOT <- script_root()
  DATA <- file.path(ROOT, "data")
  URL <- paste0(
    "https://archive.ics.uci.edu/static/public/235/",
    "individual+household+electric+power+consumption.zip"
  )
  ZIP <- file.path(DATA, "household_power_consumption.zip")
  TXT <- file.path(DATA, "household_power_consumption.txt")
  EXPECTED_ROWS <- 2075260L
  EXPECTED_BYTES <- 132960755

  dir.create(DATA, showWarnings = FALSE, recursive = TRUE)
  if (file.exists(TXT) && file.info(TXT)$size == EXPECTED_BYTES) {
    cat(sprintf("already present and correct size: %s\n", TXT))
    return(invisible(NULL))
  }
  
  options(timeout = max(600, getOption("timeout")))

  if (!file.exists(ZIP)) {
    cat(sprintf("downloading %s\n", URL))
    ok <- tryCatch({
      download.file(URL, ZIP, mode = "wb", quiet = FALSE)
      TRUE
    }, error = function(e) {
      if (file.exists(ZIP)) unlink(ZIP)
      cat(sprintf(
        "download failed: %s\nFetch it manually from https://archive.ics.uci.edu/dataset/235/ into data/\n",
        conditionMessage(e)
      ))
      FALSE
    })
    if (!ok) quit(status = 1)
  }

  cat(sprintf("unpacking %s\n", basename(ZIP)))
  unzip(ZIP, exdir = DATA)
  if (!file.exists(TXT)) {
    cat(sprintf("expected %s inside the archive but it is not there\n", TXT))
    quit(status = 1)
  }

  size <- file.info(TXT)$size
  n_rows <- count_lines(TXT)
  cat(sprintf("  size %s bytes (expected %s)\n", format(size, big.mark = ",", scientific = FALSE),
              format(EXPECTED_BYTES, big.mark = ",", scientific = FALSE)))
  cat(sprintf("  rows %s (expected %s)\n", format(n_rows, big.mark = ",", scientific = FALSE),
              format(EXPECTED_ROWS, big.mark = ",", scientific = FALSE)))
  if (size != EXPECTED_BYTES || n_rows != EXPECTED_ROWS) {
    cat("raw file does not match the documented dataset; stopping so that downstream results are not silently wrong\n")
    quit(status = 1)
  }
  cat("\nok. Next: Rscript src/01_aggregate.R\n")
}

# if (sys.nframe() == 0L) main()
main()
