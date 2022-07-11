#  
# 5' only plots have red line for TSS ? why is it using Pa settings?
# CDF plot gene list n/ brakes
# why does filter plot not always update? same number of genes?
# dendogram sub sampling/plotting
# do I need  theme(axis.title.x = element_text(size =  line_list$mysize[3], vjust = .5)) + 
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

