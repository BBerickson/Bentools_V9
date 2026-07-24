# One-time renv adoption for BenTools v9.
#
# Run this ONCE, interactively, in RStudio from the project root. It is kept as a
# separate opt-in script (not run at launch) because it reorganizes the project's
# package library and adds an renv autoloader to .Rprofile — a change to how every
# future R session in this folder starts. After running it, commit renv.lock so
# collaborators can reproduce your exact package versions with renv::restore().
#
#   1. install.packages("renv")
#   2. source("renv_setup.R")
#   3. commit the generated renv.lock
#
# Thereafter, on any machine:  renv::restore()  then  shiny::runApp()

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize renv for this project. Uses your currently-installed packages as the
# basis for the project library (no fresh downloads needed for what you already have).
renv::init(bare = FALSE, restart = FALSE)

# Record the exact versions of every package the app uses into renv.lock.
renv::snapshot(prompt = FALSE)

message("renv.lock written. Commit it to version control.")
