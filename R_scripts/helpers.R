# 
#  
# save file needs work on output format, 
    # make file name auto fill better, less details compair, filter sum, per, peak, ratio, cluster, cdf, 
    # gene_list$info up date description on all created gene lists
    # add a line to discribe columns in file
    # remove any unneeded columns, add any needed ones
    # make sure saved file will load in correctly
# loading files with complex names does not show in DT summery?
# cdf legend 20 character new line sub ... normaly happens in Active_list_data()
#    
#  filter per reactive on numeric text boxes fix crash (single/multiple) ... change behavior of all user dependent sliders and boxes? denies running + message box?
#  try and have filter plot show in better postion
# if loading a url file have progress bar count number of files completed
# # DT tab 
#    just for showing/sorting/filtering/clustering/grouping gene lists, plots:CDF,Violin,Bar,scatter,cdf
#  get big spinner running while thinkig- loading files, clustinging,   
# filter per have min and max starting numbers adjusted based on total number of genes, 
      # add spinner for waiting, include step for looking for first bit of signal (test with big data set and adjust amout going up)
#    auto looking down conditions?
# 
# #### 
##### test for crashing of what I have set up ### tool usage and tab switching picker testing/working as intended?
# size of memory used limit/warning
# have ability to change the order of legend in plots (DT sort?)
# add benefits for loading matrix files, separation and gene size?
# lines and labels TSS = bin 1 fix
# checks on pickers left side?, have then stay within view (word wrap)?
# a plot to show delta PI, MA/scatter/dot?
# # 5'  and 3' files  focused sliders optimize
# # data options and QC tab
# #     stats on files, distribution of mean signals (5,4,3), number of 0,s na's, peaks, ???
# #     ability to find and remove outliers and filter (add to server R RNAseq filtering outliers )
# in QC add gene size and seperation, deep tools inputs/bin ranges
# "factoextra" for dendogram plot fviz_dend()
# load old style table files (mod matrix?)

###### long term plans ####
# brake up functions into smaller functions so I can more easly add setprogress recycle common fun for more consistency 

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

