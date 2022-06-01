# fix common LinesLabelsPreSetGuess, and test,
# observe print("observe line and labels") ... put in logical limits of what numbers can be, add steps to num boxes? 
# for loops ... preallocate the output container ... out <- vector("list", length(x))  ... seq_along(x)
# why does filter plot not always update? same number of genes?
# add info to top line of table
# dendogram sub sampling/plotting
# test making a version w/o global LIST_DATA and put things in more reactive values with more observations? ... faster?, more control?, less memory?
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

