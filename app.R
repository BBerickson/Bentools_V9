# Created by Benjamin Erickson BBErickson@gmail.com
#
# BenTools v9 entry point. Launch with shiny::runApp() or the RStudio "Run App"
# button. This file only loads packages, sources the app's pieces, and starts the
# app; the real code lives in R/ (globals, functions, ui, server).

# Packages required by the app. See setup.R for a one-time install helper, or use
# renv::restore() if you have the committed renv.lock. Packages are NOT installed
# automatically at launch (that silently mutated the user's library); if any are
# missing the app stops with a clear message telling you how to install them.
kRequiredPackages <- c(
  "tidyverse",
  "shiny",
  "shinydashboard",
  "shinydashboardPlus",
  "shinycssloaders",
  "shinyWidgets",
  "shinyjs",
  "RColorBrewer",
  "colourpicker",
  "colorspace",
  "DT",
  "patchwork",
  "zip",
  "ggpubr",
  "ggtext",
  "fastcluster",
  "dendextend",
  "valr",
  "data.table",
  "matrixStats"
)

# load packages, or stop with actionable guidance if any are missing ----
local({
  missing <- kRequiredPackages[
    !vapply(kRequiredPackages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "Missing required package(s): ", paste(missing, collapse = ", "), ".\n",
      "Install them once with:  source(\"setup.R\")\n",
      "or, if using renv:       renv::restore()",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(
    invisible(lapply(kRequiredPackages, library, character.only = TRUE))
  )
})

# By default, the file size limit is 5MB. Raise it to 500MB for large matrices.
options(shiny.maxRequestSize = 500 * 1024 ^ 2)

# Source the app in dependency order: constants -> helpers -> ui -> server.
# NOTE: helpers live in R_scripts/ (not R/) on purpose — Shiny auto-sources any
# top-level R/ directory at startup, before this file's library() calls run, which
# would build the UI before its packages are loaded. Explicit sourcing here keeps
# load order deterministic.
source("R_scripts/globals.R", local = TRUE)
source("R_scripts/functions.R", local = TRUE)
source("R_scripts/ui.R", local = TRUE)
source("R_scripts/server.R", local = TRUE)

shinyApp(ui = ui, server = server)
