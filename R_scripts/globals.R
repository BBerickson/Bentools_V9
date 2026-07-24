# Global constants and the (legacy) shared data store.
# LIST_DATA is migrated to Shiny reactivity in a later refactor phase.

# Brewer color sets to be available ----
kBrewerList <- c("Set1","Paired","Dark2","Spectral")

LIST_DATA <<- list(
  table_file = NULL,
  # gene bin score set
  gene_file = NULL,
  # holds $Complete genes from files and $gene file(s)
  meta_data = NULL,
  # for holding meta data gene file(s) [c("gene_list", "count", "set", "color", plot?, "legend", "plot_legend")]
  ttest = NULL,
  # t.test results $full is for numbers $meta_data for holding plotting options
  meta_data_plot = list(
    binning = c(543,100,1500,3500,2000,500,500,500), # type, bp/bin, before, after, body, un5, un3, spacing
    binning2 = c(543,100,1500,3500,2000,500,500,500), # save for reset
    rnaseq = FALSE, # T/F rnaseq data type?
    landmarks = c(15, 45, 20, 40,  5), # tssbin, tesbin, body1bin, body2bin, bin spacing
    tss_tes = c("TSS", "pA"), # tss and tes labels
    x_plot_range = c(0, 0) # number of bins
  ),
  # info of matrix and lines and lables settings
  STATE = c(0, 0) # flow control
  # [1] 1 = at least one file has been loaded and lets reactive fill in info
  #
  # [2] 0 = first time switching tab auto plotting
  #     1 = hidden plot button, reactive for plot enabled
  #     2 = on/off reactive picker changed, shows plot button, reactive for plot disabled
)
