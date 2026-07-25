# One-time dependency installer for BenTools v9.
# Run once in a fresh R session:  source("setup.R")
#
# Prefer renv for reproducible versions:  renv::restore()  (uses renv.lock).
# This script is the plain-install fallback and mirrors kAttachedPackages +
# kNamespacedPackages in app.R.
#
# NOTE: valr is a plain CRAN package (as of 0.8.x/0.10.0 its only hard deps are
# CRAN: cpp11bigwig, dplyr, ggplot2, ...). Earlier versions were sometimes paired
# with Bioconductor's GenomicRanges; that is no longer required, so we just
# install valr from CRAN with everything else.

cran_packages <- c(
  "tidyverse", "shiny", "shinydashboard", "shinydashboardPlus",
  "shinycssloaders", "shinyWidgets", "shinyjs", "RColorBrewer",
  "colourpicker", "colorspace", "DT", "patchwork", "zip", "ggpubr",
  "ggtext", "fastcluster", "dendextend", "valr", "data.table", "matrixStats"
)

installed <- rownames(installed.packages())

to_install <- setdiff(cran_packages, installed)
if (length(to_install)) {
  message("Installing from CRAN: ", paste(to_install, collapse = ", "))
  install.packages(to_install, dependencies = TRUE)
} else {
  message("All required packages already installed.")
}

message("Setup complete. Launch the app with:  shiny::runApp()")
