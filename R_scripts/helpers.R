# 
# # filter/sort tool tab
# # 
# 
# hidesorttable - switch tab hide, sort/filter show, (hide switch tab goes to start)
# 
# clustering/grouping, add abilty to use #s from gene_list$full (fix lists so selection of proper data can be done)
# plot the dendeogram? and allow for more than 4 clusters.
# test using the mean of the bin range for faster clustering
# hclust.vector(df %>% filter(bin %in% c(15:16)) %>% group_by(gene) %>% summarise(value=mean(score),.groups = "drop") %>% select(value), method = "ward")
# #
# # data options and QC tab
# #     stats on files, distribution of mean signals (5,4,3), number of 0,s na's, peaks, ???
# #     ability to find and remove outliers and filter (add to server R RNAseq filtering outliers )
# 
# 
# ### change observereactive for ggplot to eventReactive??? (direct output)
# 
# do I use LIST_DATA$gene_info$count? if so fix filter
# in QC add gene size and seperation, deep tools inputs/bin ranges

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
  "generic 543"
)

# math options available ----
kMathOptions <- c("mean", "sum", "median")

