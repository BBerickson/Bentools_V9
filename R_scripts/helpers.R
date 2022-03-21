# speed up plot ... why is it running so slow?
# test groups with ttest... if groups only restrict it to ttest on groups not on all active
  #  
  # 
#  making norm file add count to gene_info
# filter peak ... clear up which is peak and which is filter area, keep/loose peak correctely labels, overlapping test fix
# 
# 
 # 
# 
# if filter has same number of genes will it update, what happens if all genes are filtered out?
# intersecting more than 2 items is not compleatly inclusive (check on other checkbox filters), update lists on run
# work on TSS bin 1 lines and labels, RF norm y axis settings fix
# save file needs work on output format (including a just gene list output option), file name, and loading files with complex names does not show in DT summury?
# cdf legend 20 character new line sub ... normaly happens in Active_list_data()
# filter tools: have search for last include type name? is that the replot problem?    
#  filter per reactive on numeric text boxes fix crash (single/multiple, and clarify min/max profile kept?) ... 
#color set update pickers  
# if loading a url file have progress bar count number of files completed
# # DT tab 
#    just for showing/sorting/filtering/clustering/grouping gene lists, plots:CDF,Violin,Bar,scatter,cdf
#  get big spinner running while thinkig- loading files, clustinging,   
# filter per have min and max numbers adjusted based on total number of genes, and set to relative 1%  filter to start, on 1 file show up and down a few steps making picked a wider line on plot
#    add note to if all bins == 0 or add auto go to next per up or down?
# smart filter tool ... set high med and low areas have filter % move # gene at a time until pattern found 
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

