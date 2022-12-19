# program for loading packages ----
my_packages <- function(x) {
  for (i in x) {
    #  require returns TRUE invisibly if it was able to load package
    if (!require(i , character.only = TRUE)) {
      #  If package was not able to be loaded then re-install
      if(i == "valr"){
        if (!require("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        
        BiocManager::install("rtracklayer")
      }
      install.packages(i , dependencies = TRUE,)
      print(paste("installing ", i, " : please wait"))
    }
    #  Load package after installing
    require(i , character.only = TRUE)
  }
}

# Brewer color sets to be available ----
kBrewerList <-
  c("Set1","Paired","Dark2","Spectral")

# lines and labels types ----
kLinesandlabels <- c(
  "543",
  "5",
  "3",
  "5L",
  "PI"
)

# math options available ----
kMathOptions <- c("mean", "sum", "median")

# hotkeys
hotkeys <- c("enter")
  