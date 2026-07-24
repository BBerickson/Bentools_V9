# One-time dependency installer for BenTools v9.
# Run once in a fresh R session:  source("setup.R")
#
# Prefer renv for reproducible versions:  renv::restore()  (uses renv.lock).
# This script is the plain-install fallback and mirrors kRequiredPackages in
# Ben_Tools.v9.R.

cran_packages <- c(
  "tidyverse", "shiny", "shinydashboard", "shinydashboardPlus",
  "shinycssloaders", "shinyWidgets", "shinyjs", "RColorBrewer",
  "colourpicker", "colorspace", "DT", "patchwork", "zip", "ggpubr",
  "ggtext", "fastcluster", "dendextend", "data.table", "matrixStats"
)

installed <- rownames(installed.packages())

to_install <- setdiff(cran_packages, installed)
if (length(to_install)) {
  message("Installing from CRAN: ", paste(to_install, collapse = ", "))
  install.packages(to_install, dependencies = TRUE)
}

# valr depends on Bioconductor's GenomicRanges
if (!requireNamespace("valr", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("GenomicRanges", update = FALSE, ask = FALSE)
  install.packages("valr", dependencies = TRUE)
}

message("Setup complete. Launch the app with:  shiny::runApp()")
