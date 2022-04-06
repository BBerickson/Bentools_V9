# #### 
##### test for crashing of what I have set up ### tool usage and tab switching picker testing/working as intended?
# 
# # QC tab
# #     stats on files, quantile plots low end and broad  ranges (1, 2.5, 5, 7.5, 10)(10, 25, 50, 75, 90) 
    # with note on signal being above , % 0,s per bin, 
#   number of genes that are all 0s, with option to remove these
# in QC add gene size and separation, deep tools inputs/bin ranges

# load old style table files (mod matrix?)

###### long term plans ####
# brake up functions into smaller functions so I can more easly add setprogress recycle common fun for more consistency 
#     Ratio :smarter labels for single file verses 2 files, CDF smarter PI/EI in info, and better info other tools
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
  "543 bins 20,20,40",
  "5' 1k 1k 80bins",
  "3' 1k 9k 100bins",
  "5' .25k 10k 205bins",
  "5'",
  "4",
  "3",
  "generic 543",
  "PI"
)

# math options available ----
kMathOptions <- c("mean", "sum", "median")

