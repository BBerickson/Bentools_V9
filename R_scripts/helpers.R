# fix factor lableing up and down includeing what is saved in gene list
# TSS i showing up red for 5' plots, ratio picker fix label
#  filter % on % crashes?
# cdf box color fix
#loading gene list has all lists in data table
# lines and lables, figuring out body/ put body in 
# dendogram sub sampling/plotting
# color change crashes on other lists? if name is*-1
# change color does not update drop downs(it seems to do it for the first change but no others), add check box on common to update all lists
# test making a version w/o global LIST_DATA and put things in more reactive values (pre delare and fix undeclared) with more observations? ... faster?, more control?, less memory?
# program for loading packages ----
my_packages <- function(x) {
  for (i in x) {
    #  require returns TRUE invisibly if it was able to load package
    if (!require(i , character.only = TRUE)) {
      #  If package was not able to be loaded then re-install
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

