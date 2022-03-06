#  activate gene list if common have $use remade with distinct
# loading a file after a gene list?, are my gene n numbers correct(inclusive/exclusive, add option)?, filter have if 0 check        
#  
# 
# # DT tab 
#    just for showing/sorting/filtering/clustering/grouping gene lists, plots:CDF,Violen,Bar,scatter,cdf
#    
# filter per have min and max numbers ajusted based on total number of genes, and set to relitive 1%  filter to start, on 1 file show up and down a few steps making picked a wider line on plot
# make filter tool to remove spike data (set peak area and check if any bin has a higher value)
# smart filter tool ... set high med and low areas have filter % move # gene at a time until pattern found 
# #### test RF and other norms to not change the sign of the data
##### test for crashing of what I have set up ### 

# test loading in grep gene lists ... updtate to like V8a
# add a test to detect _ to add regex for better pattern matching
# map(paste0("\\|" , mRNAs_up$name,"$"), str_subset, string = hg19_gene$name) 
  # MatchGenes <- function(common_list, gene_list){
    # for(g in seq_along(gene_list$gene)){
    #   gene_list$gene[g] <- str_subset(
    #     common_list$gene, gene_list$gene[g]
    #   )[1]
    # }
    # gene_list <- filter(gene_list, !is.na(gene))
    # return(gene_list)
  # }

# have ability to change the order of lists (DT function?)
# if loading a url file have progress bar count number of files completed

# a plot to show delta PI, MA/scatter/dot?
# #
# # data options and QC tab
# #     stats on files, distribution of mean signals (5,4,3), number of 0,s na's, peaks, ???
# #     ability to find and remove outliers and filter (add to server R RNAseq filtering outliers )
# 
# 
# ### change observereactive for ggplot to eventReactive??? (direct output)
# 
# 
# in QC add gene size and seperation, deep tools inputs/bin ranges
# "factoextra" for dendogram plot fviz_dend()

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

