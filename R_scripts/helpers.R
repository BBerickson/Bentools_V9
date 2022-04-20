# brake up functions into smaller functions so I can more easy add set progress recycle common fun for more consistency
# LoadTableFile - test - read in file/URL - load gene list before/after - speed - withspinner on DTs
#  add withProgress updates, reset loading file between each file as well, add examples of files that can be loaded,
# LoadGeneFile - count_fields and file type test with compatibility - read in file - intersect gene test -
#     add to LIST_DATA
# remove setProgress from functions - FindClusters, ClusterNumList, CumulativeDistribution, 
##### test for crashing of what I have set up ### tool usage and tab switching picker testing/working as intended?
# don't have the plot button show on first plot
# GGplotLineDot tidy up code if posible 
# clean up CompareRatios and consolidate ... function for gene_file$info?
# move QC to functions
# load old style table files (mod matrix?)
#     Ratio :smarter labels for single file verses 2 files, CDF smarter PI/EI in info, and better info other tools
# valuebox for peak filter (mabey the others) has total gene count not distinct, have all filters update the gene list picker, don't remove list being filtered on
# for preview plots have show with spinner before calculations
# CDF with PI proper num denom selection range, add controls for the number of observations on the dot plot
# fix save file help text, CDF save gather so sets are column names?, fix cluster select sample label  
# add AUC tool
# norm -1 add making negitive for antisense plots
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
  "5",
  "4",
  "3",
  "generic 543",
  "PI"
)

# math options available ----
kMathOptions <- c("mean", "sum", "median")

