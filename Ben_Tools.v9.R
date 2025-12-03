# Created by Benjamin Erickson BBErickson@gmail.com

# program for loading packages ----
my_packages <- function(x) {
  for (i in x) {
    #  require returns TRUE invisibly if it was able to load package
    if (!require(i , character.only = TRUE)) {
      #  If package was not able to be loaded then re-install
      # if(i == "valr"){
      #   if (!require("BiocManager", quietly = TRUE))
      #     install.packages("BiocManager")
      #   BiocManager::install("SparseArray")
      #   BiocManager::install("rtracklayer")
      # }
      install.packages(i , dependencies = TRUE,)
      print(paste("installing ", i, " : please wait"))
    }
    #  Load package after installing
    require(i , character.only = TRUE)
  }
}

# run load needed packages using my_packages(x) ----
suppressPackageStartupMessages(my_packages(
  c(
    "tidyverse",
    "shiny",
    "shinydashboard",
    "shinydashboardPlus",
    "shinycssloaders",
    "shinyWidgets",
    "shinyjs",
    "RColorBrewer",
    "colourpicker",
    "colorspace",
    "DT",
    "patchwork",
    "zip",
    "ggpubr",
    "ggtext",
    "fastcluster",
    "dendextend",
    "valr"
  )
))

source("R_scripts/functions.R", local = TRUE)

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 500MB. ----
options(shiny.maxRequestSize = 500 * 1024 ^ 2)

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

# server ----
server <- function(input, output, session) {
  # remove on non-local deployment
  session$onSessionEnded(stopApp)
  tt <- tibble(gene="chr1:10-100-;NM_Name|YFG",bin=1:3,score=c(.1,2,2.2))
  dt2 <- datatable(
    tt,
    caption = "example .table file",
    colnames = c("","",""),
    rownames = FALSE,
    filter = "none",
    class = 'cell-border stripe compact',
    options = list(
      scrollX = TRUE,
      scrollY = TRUE,
      autoWidth = T,
      width = 30,
      sDom  = '<"top">lrt<"bottom">ip',
      info = FALSE,
      paging = FALSE,
      lengthChange = FALSE,
      ordering = FALSE,
      columnDefs = list(
        list(className = 'dt-center ', targets = "_all"),
        list(
          targets = 0,
          render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 25 ?",
            "'<span title=\"' + data + '\">' + data.substr(0, 27) + '...</span>' : data;",
            "}"
          )
        )
      )
    )
  )
  output$loadedfilestotaltable <- DT::renderDataTable(dt2)
  # reactive values ----
  reactive_values <- reactiveValues(
    Apply_Math = NULL,
    Plot_Options = NULL,
    Y_Axis_numbers = c(0,0),
    Lines_Labels_List = list(mybrakes="",mylabels="",
                             myline = tibble(use_virtical_line_color=c("green","red","black","black"),
                                             use_virtical_line_type=c("dotted","dotted","solid","solid")),
                             mysize = c(2.0, 2.5, 13.0, 13.0, 10.0, 0.8, 20, 5),
                             myset = c(20, 40, 15, 45, 100, 5)),
    Picker_controler = NULL,
    mymath = c("mean", "none", "0", "FALSE", "FALSE", "1", "80", "none"),
    ttest = NULL,
    ttest_values = c("none", "wilcox.test", "two.sided", "FALSE", "FALSE", "-log10", "fdr"),
    ttest_options = c(0, 0, 1, "select sample", 0.05),
    clustergroups = NULL,
    groupiesgroups = NULL,
    cluster_control = NULL,
    groupies_control = NULL,
    setsliders = NULL
  )
  
  output$user <- renderUser({
    dashboardUser(
      name = "BenTools V9.m",
      image = "ben head.jpg",
      title = "Benjamin Erickson",
      subtitle = "BBErickson@gmail.com"
    )
  })
  
  # disables tabs on start ----
  addCssClass(selector = "a[data-value='mainplot']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='qcOptions']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='filenorm']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='grouptab']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='genelists']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='sorttool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='ratiotool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='clustertool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='groupiestool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='cdftool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='DataTableTool']", class = "inactiveLink")
  
  # loads data file(s) ----
  observeEvent(input$filetable, {
    # print("load file")
    shinyjs::disable("startoff")
    shinyjs::disable("startoff2")
    shinyjs::show("hidespiners")
    output$loadedfilestotaltable <- DT::renderDataTable(datatable())
    withProgress(message = 'loading file(s) in progress',
                 detail = 'This may take a sec...',
                 value = 0,
                 {
                   # preps/makes meta data: group, nickname,  rgb/hex plot color
                   meta_data <- PrepMetaFile(input$filetable$datapath,
                                             input$filetable$name)
                   if (!is_empty(meta_data)) {
                     for (i in seq_along(meta_data$filepath)) {
                       setProgress(i/length(meta_data$filepath), 
                                   detail = paste("Gathering info on",meta_data$nick[i]))
                       # check unique nicknames and fix
                       if (!is_empty(LIST_DATA$table_file$set) && meta_data$nick[i] %in% distinct(LIST_DATA$table_file,set)$set) {
                         meta_data$nick[i] <- paste0(meta_data$nick[i],"_dup")
                       }
                       # check file can be accessed and has is a matrix file
                       bin_colname <- tableTestbin(meta_data[i, ])
                       if(is_empty(bin_colname)){
                         showModal(
                           modalDialog(
                             title = "Information message",
                             paste("File", meta_data$nick[i], "not there, or wrong file type"),
                             size = "s",
                             easyClose = TRUE
                           )
                         )
                         if(i == length(meta_data$filepath)){
                           return()
                         } else {
                           next()
                         }
                       }
                       
                       setProgress(i/length(meta_data$filepath), 
                                   detail = paste("Testing compatiblity",meta_data$nick[i]))
                       # if first file loaded, save info else test file compatibility 
                       if (LIST_DATA$meta_data_plot$x_plot_range[2] == 0) {
                         LIST_DATA$meta_data_plot$x_plot_range <<- c(1,bin_colname$num_bins-6)
                         LIST_DATA$meta_data_plot$landmarks <<- LinesLableLandmarks(bin_colname$binning)
                         # for testing if files are loaded later and resetting lines and labels
                         LIST_DATA$meta_data_plot$binning <<- bin_colname$binning
                         LIST_DATA$meta_data_plot$binning2 <<- bin_colname$binning
                         LIST_DATA$meta_data_plot$rnaseq <<- bin_colname$rnaseq
                       } else if (max(bin_colname$num_bins-6) != LIST_DATA$meta_data_plot$x_plot_range[2] | 
                                  all(c(bin_colname$binning) != LIST_DATA$meta_data_plot$binning2)) {
                         showModal(
                           modalDialog(
                             title = "Information message",
                             "Can't load file, different number of bins or bin sizes",
                             size = "s",
                             easyClose = TRUE
                           )
                         )
                         if(i == length(meta_data$filepath)){
                           return()
                         } else {
                           next()
                         }
                       }
                       shinyjs::reset("filetable")
                       setProgress(i/length(meta_data$filepath), 
                                   detail = paste("Loading file",meta_data$nick[i]))
                       LD <- LoadTableFile(meta_data[i, ], bin_colname)
                       if(is.na(distinct(LD,set))){
                         if(i == length(meta_data$filepath)){
                           return()
                         } else {
                           next()
                         }
                       }
                       if (!is_empty(LIST_DATA$gene_file)){
                         gene_names <-
                           semi_join(LD, LIST_DATA$gene_file$Complete$full, by = "gene") %>% distinct(gene)
                         if (n_distinct(gene_names$gene) == 0) {
                           showModal(
                             modalDialog(
                               title = "Information message",
                               " No genes in common ",
                               size = "s",
                               easyClose = TRUE
                             )
                           )
                         }
                         setProgress(i/length(meta_data$filepath), 
                                     detail = paste("Finishing up",meta_data$nick[i]))
                         # make complete gene list
                         LIST_DATA$gene_file$Complete$full <<-
                           full_join(LD, LIST_DATA$gene_file$Complete$full, by = c("gene","chrom", "start", "end","strand")) %>%
                           distinct(gene, chrom, start, end, strand)
                         LIST_DATA$meta_data <<- LIST_DATA$meta_data %>%
                           dplyr::mutate(count = if_else(gene_list == "Complete",
                                                         paste("n =", n_distinct(LIST_DATA$gene_file$Complete$full$gene, na.rm = T)), count))
                         
                       } else {
                         LIST_DATA$gene_file$Complete$full <<- distinct(LD, gene, chrom, start, end, strand)
                       }
                       if (LIST_DATA$STATE[2] == 0 &
                           n_distinct(LIST_DATA$table_file$set) < 2) {
                         oo <- meta_data$nick[i]
                       } else {
                         oo <- "0"
                       }
                       LIST_DATA$table_file <<- distinct(bind_rows(LIST_DATA$table_file, LD)) %>% filter(!is.na(set))
                       LIST_DATA$gene_file$Complete$info <<- tibble(loaded_info = paste("all loaded genes",
                                                                                        Sys.Date()))
                       LIST_DATA$meta_data <<- distinct(bind_rows(LIST_DATA$meta_data,tibble(
                         gene_list = "Complete",
                         count = paste("n =", n_distinct(LIST_DATA$gene_file$Complete$full$gene, na.rm = T)),
                         set = meta_data$nick[i],
                         mycol = meta_data$color[i],
                         onoff = oo,
                         sub = " ",
                         plot_legend = " ",
                         group = meta_data$group[i]
                       )))
                       ### check if grep of gene has occurred
                       for(gg in seq_along(LIST_DATA$gene_file)[-1]){
                         if("org_gene" %in% names(LIST_DATA$gene_file[[gg]]$full)){
                           setProgress(i/length(meta_data$filepath), 
                                       detail = paste("re-matching gene list",meta_data$nick[i]))
                           if(str_detect(LIST_DATA$meta_data$gene_list[gg], "_intersected")){
                             new_gene_match <- MatchGenes(LIST_DATA$gene_file[[1]]$full,
                                                          LIST_DATA$gene_file[[gg]]$full %>% select(-gene) %>% 
                                                            dplyr::rename(gene = org_gene), bedfile = T)
                           } else {
                             new_gene_match <- MatchGenes(LIST_DATA$gene_file[[1]]$full,
                                                          LIST_DATA$gene_file[[gg]]$full %>% select(org_gene) %>%
                                                            dplyr::rename(gene = org_gene))
                           }
                           
                           if (n_distinct(new_gene_match$gene, na.rm = T) != 0) {
                             # fix name, fix info
                             LIST_DATA$gene_file[[gg]]$full <<- new_gene_match
                             listname <- names(LIST_DATA$gene_file)[gg]
                             LIST_DATA$meta_data <<- LIST_DATA$meta_data %>%
                               dplyr::mutate(count=if_else(gene_list == listname,paste("n =",n_distinct(new_gene_match$gene, na.rm = T)),
                                                           count))
                           }
                         }
                         setProgress(i/length(meta_data$filepath), 
                                     detail = paste("Updating gene lists",meta_data$nick[i]))
                         # add missing data
                         LIST_DATA$meta_data <<- LIST_DATA$meta_data %>%
                           dplyr::filter(gene_list == names(LIST_DATA$gene_file)[gg]) %>%
                           dplyr::mutate(set = meta_data$nick[i],
                                         count = count[i],
                                         mycol = meta_data$color[i],
                                         onoff = "0",
                                         sub = " ",
                                         plot_legend = " ") %>%
                           bind_rows(LIST_DATA$meta_data, .)
                       }
                     }
                   } else {
                     return()
                   }
                 })
    # filter overlapping regions, binning = (type, bp/bin, before, after, body, un5, un3, spacing)
    sep_bins <- 1
    gene_length <- 10
    if(LIST_DATA$meta_data_plot$binning2[1] == "543"){
      sep_bins <- max(LIST_DATA$meta_data_plot$binning2[c(3,4)],1)
      gene_length <- max(sum(LIST_DATA$meta_data_plot$binning2[c(6,7)]),gene_length)
    } else if(LIST_DATA$meta_data_plot$binning2[1] == "5"){
      sep_bins <- max(LIST_DATA$meta_data_plot$binning2[3],1)
      gene_length <- max(sum(LIST_DATA$meta_data_plot$binning2[4]),gene_length)
    } else if(LIST_DATA$meta_data_plot$binning2[1] == "3"){
      sep_bins <- max(LIST_DATA$meta_data_plot$binning2[4],1)
      gene_length <- max(sum(LIST_DATA$meta_data_plot$binning2[3]),gene_length)
    } 
    updateNumericInput(session, "geneSizeMin", value = gene_length) 
    updateNumericInput(session, "geneSeparation", value = sep_bins)
    LD <- FilterSepSize(distinct(LIST_DATA$table_file,gene,chrom,start,end,strand),
                        sep_bins,
                        gene_length,
                        0,
                        LIST_DATA$meta_data_plot$rnaseq)
    
    if(n_distinct(LD$gene) > 0) {
      # adds full n count to nickname
      
      # preps meta data
      LIST_DATA$meta_data <<- LIST_DATA$meta_data %>% 
        dplyr::mutate(count = if_else(gene_list == "Complete",
                                      paste("n =", n_distinct(LD$gene, na.rm = T)),
                                      count),
                      sub = if_else(gene_list == "Complete",
                                    paste0("Complete: sep > ",sep_bins,", length > ", gene_length),
                                    sub)
                                    )
      # saves data in list of lists
      LIST_DATA$gene_file[["Complete"]]$full <<- distinct(LD)
    }
    # first time starting
    if (LIST_DATA$STATE[1] == 0) {
      shinyjs::show("startoff")
      shinyjs::show("startoff2")
      # sets lines and labels
      # tries to guess lines and labels type
      reactive_values$slider_breaks <- LinesLabelsSet(LIST_DATA$meta_data_plot$binning,
                                                      LIST_DATA$meta_data_plot$landmarks,
                                                      LIST_DATA$meta_data_plot$x_plot_range[2],
                                                      slider = T)
      reactive_values$setsliders <- SlidersSetsInfo(reactive_values$slider_breaks, LIST_DATA$meta_data_plot$binning[1])
      if(LIST_DATA$meta_data_plot$rnaseq){
        updateAwesomeCheckbox(session,inputId = "checkboxfull",value = TRUE)
        updateAwesomeRadio(session,inputId = "checkboxbin",selected = "subtract")
      }
    } 
    # enables tabs after loading file
    shinyjs::enable("startoff")
    if(length(names(LIST_DATA$gene_file))>1){
      shinyjs::enable("startoff2")
    }
    shinyjs::reset("filetable")
    removeCssClass(selector = "a[data-value='qcOptions']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='mainplot']", class = "inactiveLink")
    
    gts <- LIST_DATA$table_file %>% group_by(set) %>%
      summarise(number_of_genes = n_distinct(gene, na.rm = T)) %>%
      filter(!is.na(set)) %>% 
      rename(File = set)
    dt <- datatable(
      gts,
      colnames = names(gts),
      rownames = FALSE,
      filter = "none",
      class = 'cell-border stripe compact',
      options = list(
        scrollX = TRUE,
        scrollY = TRUE,
        autoWidth = TRUE,
        width = 20,
        sDom  = '<"top">lrt<"bottom">ip',
        info = FALSE,
        paging = FALSE,
        lengthChange = FALSE,
        columnDefs = list(
          list(className = 'dt-center ', targets = "_all"),
          list(
            targets = 0,
            render = JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 25 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 27) + '...</span>' : data;",
              "}"
            )
          )
        )
      )
    )
    output$loadedfilestable <- DT::renderDataTable(dt)
    gts2 <-
      LIST_DATA$gene_file$Complete$full %>% summarise(total_number_distinct_genes = n_distinct(gene, na.rm = T))
    dt2 <- datatable(
      gts2,
      colnames = names(gts2),
      rownames = FALSE,
      filter = "none",
      class = 'cell-border stripe compact',
      options = list(
        scrollX = TRUE,
        scrollY = TRUE,
        autoWidth = T,
        width = 30,
        sDom  = '<"top">lrt<"bottom">ip',
        info = FALSE,
        paging = FALSE,
        lengthChange = FALSE,
        ordering = FALSE,
        columnDefs = list(list(
          className = 'dt-center ', targets = "_all"
        ))
      )
    )
    output$loadedfilestotaltable <- DT::renderDataTable(dt2)
    if (length(LIST_DATA$gene_file) > 1) {
      gg <- LIST_DATA$meta_data %>% dplyr::filter(!str_detect(gene_list,"^Complete")) %>%
        select(., gene_list, count) %>% dplyr::rename(Usable = count) %>% 
        distinct()
      ggg <- NULL
      for (i in gg$gene_list) {
        ggg <-
          c(
            ggg,
            sapply(LIST_DATA$gene_file[i], function(x) x$full$gene) %>% bind_cols(.) %>% suppressMessages() %>% n_distinct(1,na.rm = T)
          )
      }
      ggg <- dplyr::mutate(gg, "total_in_file" = ggg)
      dtg <- datatable(
        ggg,
        colnames = names(ggg),
        rownames = FALSE,
        filter = "none",
        class = 'cell-border stripe compact',
        options = list(
          scrollX = TRUE,
          scrollY = TRUE,
          autoWidth = TRUE,
          width = 20,
          sDom  = '<"top">lrt<"bottom">ip',
          info = FALSE,
          paging = FALSE,
          lengthChange = FALSE,
          columnDefs = list(
            list(className = 'dt-center ', targets = "_all"),
            list(targets = 0, width = 20)
          )
        )
      )
      output$loadedgenetable <- DT::renderDataTable(dtg)
    }
  })
  
  
  # loads gene list file ----
  observeEvent(input$filegene1, {
    # print("load gene file")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   # load info, update select boxes, switching works and changing info and plotting
                   LD <- LoadGeneFile(input$filegene1$datapath,
                                      input$filegene1$name,
                                      LIST_DATA,
                                      input$checkboxgenematch)
                 })
    if (!is.null(LD)) {
      LIST_DATA <<- LD
    }
    shinyjs::reset("filegene1")
    
    gg <-
      LIST_DATA$meta_data %>% dplyr::filter(!str_detect(gene_list,"^Complete|^Filter|^Gene_List_|^Ratio_|^Cluster_|^Groups_|^CDF")) %>%
      select(., gene_list, count) %>% dplyr::rename(Usable = count) %>%
      distinct()
    ggg <- NULL
    for (i in gg$gene_list) {
      ggg <-
        c(
          ggg,
          sapply(LIST_DATA$gene_file[i], function(x) x$full$gene) %>% bind_cols(.) %>% suppressMessages() %>% n_distinct(1,na.rm = T)
        )
    }
    ggg <- dplyr::mutate(gg, "total_in_file" = ggg)
    dtg <- datatable(
      ggg,
      colnames = names(ggg),
      rownames = FALSE,
      filter = "none",
      class = 'cell-border stripe compact',
      options = list(
        scrollX = TRUE,
        scrollY = TRUE,
        autoWidth = TRUE,
        width = 20,
        sDom  = '<"top">lrt<"bottom">ip',
        info = FALSE,
        paging = FALSE,
        lengthChange = FALSE,
        columnDefs = list(
          list(className = 'dt-center ', targets = "_all"),
          list(targets = 0, width = 20)
        )
      )
    )
    output$loadedgenetable <- DT::renderDataTable(dtg)
    shinyjs::enable("startoff2")
    updateSelectInput(session, "selectsave",
                      choices = names(LIST_DATA$gene_file)[-1],
                      selected = names(LIST_DATA$gene_file)[2])
    if( LIST_DATA$STATE[1] != 0){
      LIST_DATA$STATE[1] <<- .5
    }
  })
  
  # save functions, gene list ----
  output$downloadGeneList <- downloadHandler(
    filename = function() {
      paste0(LIST_DATA$gene_file[[input$selectsave]]$info$save_name,
             ".tsv")
    },
    content = function(file) {
      new_comments <-
        new_comments <- paste("#", Sys.Date(), "\n")
      new_comments <-
        c(new_comments,  paste("\n#", gsub("\nn = ", " n = ",  input$selectsave)))
      new_comments <-
        c(new_comments, paste("#", gsub(
          "\nn = ", " n = ",
          paste(LIST_DATA$gene_file[[input$selectsave]]$info$loaded_info)
        )))
      new_comments <-
        c(new_comments, paste("#", 
                              paste(LIST_DATA$gene_file[[input$selectsave]]$info$col_info)
        ))
      new_comments2 <-
        inner_join(LIST_DATA$gene_file$Complete$full,LIST_DATA$gene_file[[input$selectsave]]$full)
        LIST_DATA$gene_file[[input$selectsave]]$full
      if(input$selectsave == "CDF Log2 PI Cumulative plot"){
        new_comments2 <- new_comments2 %>% 
          select(-plot_legend, -bin) %>% 
          spread(set,value)
      }
      write_lines(new_comments, file)
      write_tsv(new_comments2,
                file,
                col_names = TRUE,
                append = T)
      
    }
  )
  
  # observe actionmyplot, Lines_Labels_List, myplot update apply_Math ----
  observeEvent(c(input$actionmyplot, reactive_values$Lines_Labels_List, reactive_values$myplot), ignoreInit = TRUE, {
    # print("plot button")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   list_data_frame <- Active_list_data(LIST_DATA,input$mygroup, input$checkboxfull, 
                                                       input$selectlegendnewline, input$selectlegendnewlinespace)
                   if (!is_empty(list_data_frame) & nrow(list_data_frame) !=0) {
                     LIST_DATA$meta_data <<- rows_update(LIST_DATA$meta_data,
                                                         list_data_frame %>% 
                                                           distinct(set,gene_list,plot_legend), 
                                                         by=c("set","gene_list"))
                     reactive_values$Apply_Math <- ApplyMath(
                       list_data_frame,
                       input$myMath,
                       input$selectplotnrom,
                       as.numeric(floor(reactive_values$slider_breaks$mybrakes[
                         reactive_values$slider_breaks$mylabels %in% input$selectplotBinNorm])),
                       input$selectplotBinNorm,
                       input$mygroup,
                       input$checkboxabs,
                       input$checkboxbin
                     )
                     sliderplotBinRange <- floor(reactive_values$slider_breaks$mybrakes[
                       reactive_values$slider_breaks$mylabels %in% input$sliderplotBinRange])
                     if(length(sliderplotBinRange) < 2){
                       sliderplotBinRange <- floor(reactive_values$slider_breaks$mybrakes[
                         reactive_values$slider_breaks$mylabels %in% reactive_values$slider_breaks$myselect])
                     } else {
                       sliderplotBinRange <- range(sliderplotBinRange,na.rm = T)
                     }
                     reactive_values$Y_Axis_numbers <-
                       YAxisValues(
                         reactive_values$Apply_Math,
                         sliderplotBinRange,
                         input$checkboxlog2
                       )
                     reactive_values$Y_Axis_numbers_set <- reactive_values$Y_Axis_numbers
                     if(!is.null(reactive_values$ttest)){
                       if(input$switchttest == "by lists" & 
                          n_distinct(list_data_frame$gene_list, na.rm = T) == 1 | 
                          input$switchttest == "by files" & n_distinct(list_data_frame$set, na.rm = T) == 1) {
                         updatePickerInput(session, "switchttest", selected = "none")
                       }
                       if (input$switchttest != "none"){ 
                         LIST_DATA$ttest <<- ApplyTtest(list_data_frame,
                                                        input$switchttest,
                                                        input$selectttestlog,
                                                        input$switchttesttype,
                                                        input$padjust,
                                                        input$selectttestalt,
                                                        input$selectttestexact,
                                                        input$selectttestpaired,
                                                        input$mygroup,
                                                        input$selectlegendnewline, 
                                                        input$selectlegendnewlinespace)
                       } else {
                         LIST_DATA$ttest <<- NULL
                       }
                     }
                     reactive_values$Plot_Options <- NULL
                     reactive_values$Plot_Options <-
                       MakePlotOptionFrame(LIST_DATA$meta_data)
                     
                   } else {
                     LIST_DATA$STATE[2] <<- 2
                     text = paste("Nothing selected to plot.\n")
                     reactive_values$Plot_controler <- ggplot() +
                       annotate(
                         "text",
                         x = 4,
                         y = 25,
                         size = 8,
                         label = text
                       ) +
                       theme_void()
                     return()
                   }
                 })
    shinyjs::hide("actionmyplotshow")
  })
  
  # updates mymath plot ----
  observeEvent(input$actionMathUpDatePlot, ignoreInit = T, {
    if(between(input$numericsmooth,0.0001,1)){
      numericsmooth <- input$numericsmooth
    } else {
      updateNumericInput(session,"numericsmooth",value = 0.1)
      numericsmooth <- 0.1
    }
    # print("actionMathUpDatePlot")
    mymath <- c(input$myMath,
                input$selectplotnrom,
                input$selectplotBinNorm,
                input$checkboxsmooth,
                numericsmooth,
                input$checkboxlog2,
                input$sliderplotBinRange,
                input$mygroup,
                input$checkboxauc,
                input$checkboxabs
    )
    Y_Axis_numbers <-
      c(input$numericYRangeLow,input$numericYRangeHigh)
    if(!identical(reactive_values$mymath, mymath)){
      reactive_values$mymath <- mymath
      reactive_values$myplot <- mymath
    } else if(!identical(reactive_values$Y_Axis_numbers, Y_Axis_numbers)){
      reactive_values$Y_Axis_numbers <- Y_Axis_numbers
      if (!is_empty(reactive_values$Apply_Math)) {
        reactive_values$Plot_Options <- NULL
        reactive_values$Plot_Options <-
          MakePlotOptionFrame(LIST_DATA$meta_data)
      }
    }
    updateBoxSidebar(id = "sidebarmath")
  })
  
  # reactive Apply_Math, sets Y axis min max ----
  observeEvent(reactive_values$Y_Axis_numbers_set, ignoreInit = T, ignoreNULL = T, {
    # print("updates reactive_values$Y_Axis_numbers")
    my_step <-
      (max(reactive_values$Y_Axis_numbers) - min(reactive_values$Y_Axis_numbers)) /
      20
    updateNumericInput(session,
                       "numericYRangeHigh",
                       value = round(max(reactive_values$Y_Axis_numbers), 4),
                       step = my_step)
    updateNumericInput(session,
                       "numericYRangeLow",
                       value = round(min(reactive_values$Y_Axis_numbers), 4),
                       step = my_step)
  })
  
  # reactive Violin, sets Y axis min max ----
  observeEvent(reactive_values$Y_Axis_numbersViolin,
               ignoreInit = T, ignoreNULL = T, {
    # print("violin min max set")
    my_step <-
      (max(reactive_values$Y_Axis_numbersViolin) - min(reactive_values$Y_Axis_numbersViolin)) /
      20
    updateNumericInput(session,
                       "numericYRangeHighViolin",
                       value = round(max(reactive_values$Y_Axis_numbersViolin), 4),
                       step = my_step)
    updateNumericInput(session,
                       "numericYRangeLowViolin",
                       value = round(min(reactive_values$Y_Axis_numbersViolin), 4),
                       step = my_step)
  })
  
  # keep violin merge bin in check
  observeEvent(input$VmergeBin, ignoreInit = T, ignoreNULL = T, {
    VmergeBin <- input$VmergeBin
    
    # Check if numeric and greater than or equal to 1
    if(!is.numeric(VmergeBin) || VmergeBin < 1) {
      VmergeBin <- 1
      updateNumericInput(session, "VmergeBin", value = 1)
    } 
    # Check if it's a whole number (not a decimal)
    else if(VmergeBin != floor(VmergeBin)) {
      VmergeBin <- floor(VmergeBin)
      updateNumericInput(session, "VmergeBin", value = floor(VmergeBin))
    }
  })
  
  observeEvent(c(reactive_values$Y_Axis_numbers_set_Violin, input$VplotType),
               ignoreInit = T, ignoreNULL = T, {
    # print("updates plot")
    
    Y_Axis_Label <- YAxisLabel(input$myMath,
                               input$selectplotnrom,
                               input$selectplotBinNorm,
                               FALSE,
                               input$checkboxlog2Violin)
    
    sliderplotBinRange <- floor(reactive_values$slider_breaks$mybrakes[
      reactive_values$slider_breaks$mylabels %in% input$sliderplotBinRange])
    if(length(sliderplotBinRange) < 2){
      sliderplotBinRange <- floor(reactive_values$slider_breaks$mybrakes[
        reactive_values$slider_breaks$mylabels %in% reactive_values$slider_breaks$myselect])
    } else {
      sliderplotBinRange <- range(sliderplotBinRange,na.rm = T)
    }
    
    list_data_frame <- Active_list_data(LIST_DATA,input$mygroup, input$checkboxfull, 
                                        input$selectlegendnewline, input$selectlegendnewlinespace) %>% 
      dplyr::mutate(set = plot_legend)
    
    Y_Axis_numbersViolin <-
      c(input$numericYRangeLowViolin,input$numericYRangeHighViolin)
    if(sum(Y_Axis_numbersViolin,na.rm = T)==0){
      Y_Axis_numbersViolin <- reactive_values$Y_Axis_numbersViolin
    }
    LLset <- LinesLabelsSet(LIST_DATA$meta_data_plot$binning,
                             LIST_DATA$meta_data_plot$landmarks,
                             LIST_DATA$meta_data_plot$x_plot_range[2],
                             slider = T)
    
    reactive_values$Plot_controler_Violin <-
      GGplotBoxViolin(list_data_frame,
                      sliderplotBinRange,
                      reactive_values$Plot_Options,
                      Y_Axis_numbersViolin,
                      reactive_values$Lines_Labels_List,
                      LLset,
                      input$checkboxlog2Violin,
                      Y_Axis_Label,
                      bin_step = input$VmergeBin,  # Aggregate bins to reduce clutter
                      plot_type = input$VplotType) 
    
  })
  
  # reactive Plot_Options, triggers plot ----
  observeEvent(reactive_values$Plot_Options, ignoreInit = T, ignoreNULL = T, {
    # print("plot")
    if(!is.null(LIST_DATA$ttest)){
      mm <- reactive_values$ttest_options[1:2]
      Plot_Options_ttest <- MakePlotOptionttest(LIST_DATA$ttest, as.numeric(mm),
                                                input$selectttestlog,input$hlinettest,input$padjust,input$switchttesttype)
    }
    Y_Axis_Label <- YAxisLabel(input$myMath,
                               input$selectplotnrom,
                               input$selectplotBinNorm,
                               input$checkboxsmooth,
                               input$checkboxlog2)
    Lines_Labels_List <- reactive_values$Lines_Labels_List
    
    sliderplotBinRange <- floor(reactive_values$slider_breaks$mybrakes[
      reactive_values$slider_breaks$mylabels %in% input$sliderplotBinRange])
    if(length(sliderplotBinRange) < 2){
      sliderplotBinRange <- floor(reactive_values$slider_breaks$mybrakes[
        reactive_values$slider_breaks$mylabels %in% reactive_values$slider_breaks$myselect])
    } else {
      sliderplotBinRange <- range(sliderplotBinRange,na.rm = T)
    }
    
    reactive_values$Plot_controler <-
      GGplotLineDot(
        reactive_values$Apply_Math,
        sliderplotBinRange,
        reactive_values$Plot_Options,
        reactive_values$Y_Axis_numbers,
        Lines_Labels_List,
        input$checkboxsmooth,
        input$numericsmooth,
        Plot_Options_ttest,
        input$checkboxlog2,
        Y_Axis_Label,
        input$sliderplotOccupancy,
        input$checkboxauc
      )
    
    LIST_DATA$STATE[2] <<- 1
  })
  
  # reactive Violin Plot, triggers plot ----
  observeEvent(input$actionViolinPlot, ignoreInit = T, ignoreNULL = T, {
    # print("Violin action button")
    
    sliderplotBinRange <- floor(reactive_values$slider_breaks$mybrakes[
      reactive_values$slider_breaks$mylabels %in% input$sliderplotBinRange])
    if(length(sliderplotBinRange) < 2){
      sliderplotBinRange <- floor(reactive_values$slider_breaks$mybrakes[
        reactive_values$slider_breaks$mylabels %in% reactive_values$slider_breaks$myselect])
    } else {
      sliderplotBinRange <- range(sliderplotBinRange,na.rm = T)
    }
    
    list_data_frame <- Active_list_data(LIST_DATA,input$mygroup, input$checkboxfull, 
                                        input$selectlegendnewline, input$selectlegendnewlinespace) %>% 
      dplyr::mutate(set = plot_legend)
    
    reactive_values$Y_Axis_numbersViolin <-
      YAxisValues(
        list_data_frame %>% mutate(min=min(score,na.rm=T),max=max(score,na.rm=T)),
        sliderplotBinRange,
        input$checkboxlog2Violin
      )
    
    LIST_DATA$STATE[2] <<- 2
    
    reactive_values$Y_Axis_numbers_set_Violin <- c(reactive_values$Y_Axis_numbersViolin,input$VmergeBin,
                                                   input$numericYRangeLowViolin,input$numericYRangeHighViolin,
                                                   LIST_DATA$STATE[2])
    
  })
  
  # checks that number of names == position ----
  observeEvent(c(input$landlnames, input$landlposition), ignoreInit = TRUE, {
    if(LIST_DATA$STATE[2] > 0){
      # print("names == position")
      my_pos <-
        suppressWarnings(as.numeric(unlist(
          strsplit(input$landlposition, split = "\\s+")
        )))
      my_label <- unlist(strsplit(input$landlnames, split = "\\s+"))
      
      if (length(my_pos) == length(my_label)) {
        shinyjs::enable("actionlineslabels")
        updateActionButton(session, "actionlineslabels", label = "SET and Plot")
      } else {
        updateActionButton(session, "actionlineslabels", label = "Labels must equel # of positions")
        shinyjs::disable("actionlineslabels")
      }
    }
  }) 
  
  # dropcolor opens color select dialog box ----
  observeEvent(c(input$dropcolor), ignoreInit = T, {
    # print("dropcolor")
    showModal(modalDialog(
      title = "Information message",
      " Update Nickname and color of samples",
      size = "l",
      easyClose = F,
      footer = tagList(
        box(
          title = "File Options",
          solidHeader = T,
          width = 12,
          box(
            title =  "Set Plot Color Options",
            width = 4,
            status = "navy",
            solidHeader = T,
            fluidRow(
              box(
                width = 12,
                solidHeader = T,
                status = "info",
                background = "light-blue",
                colourInput("colourhex", "Select color HEX",value = distinct(LIST_DATA$meta_data,mycol)$mycol[1]),
                tags$hr(),
                textInput("textrgbtohex", "RGB", value = RgbToHex(x = distinct(LIST_DATA$meta_data,mycol)$mycol[1], convert = "rgb")),
                tags$hr(),
                actionButton("actionmyrgb", "Update color",width = 100)
              )
            )
          ),
          box(
            title =  "Set Plot Options",
            width = 8,
            status = "navy",
            solidHeader = T,
            pickerInput("selectgenelistoptions", "", 
                        width = 300, 
                        choices = distinct(LIST_DATA$meta_data,gene_list)$gene_list,
                        selected = distinct(LIST_DATA$meta_data,gene_list)$gene_list[1],
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$meta_data,gene_list)$gene_list)
                        )),
            pickerInput("selectdataoption", "", choices = distinct(LIST_DATA$meta_data,set)$set, 
                        selected = distinct(LIST_DATA$meta_data,set)$set[1],
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$meta_data,set)$set)
                        )),
            tags$hr(style = "color: #2e6da4; background-color: #2e6da4; border-color: #2e6da4;"),
            textInput("textnickname", "Update Nickname",value = distinct(LIST_DATA$meta_data,set)$set[1]),
            actionButton("actionoptions", "Set Nickname"),
            helpText("Need to press to update")
          )
        ),
        modalButton("Close")
      )
    ))
  })
  
  # CDF color select dialog box ----
  observeEvent(c(input$actioncdfcolor), ignoreInit = T, {
    gene_info <- LIST_DATA$meta_data %>% 
      dplyr::filter(str_detect(gene_list,"^CDF"))
    # print("actioncdfcolor")
    showModal(modalDialog(
      title = "Information message",
      " Update color of samples",
      size = "l",
      easyClose = F,
      footer = tagList(
        box(
          title = "File Options",
          solidHeader = T,
          width = 12,
          box(
            title =  "Set Plot Color Options",
            width = 4,
            status = "navy",
            solidHeader = T,
            fluidRow(
              box(
                width = 12,
                solidHeader = T,
                status = "info",
                background = "light-blue",
                colourInput("colourhex", "Select color HEX",value = distinct(gene_info,mycol)$mycol[1]),
                tags$hr(),
                textInput("textrgbtohex", "RGB", value = RgbToHex(x = distinct(gene_info,mycol)$mycol[1],
                                                                  convert = "rgb")),
                tags$hr(),
                actionButton("actionmyrgb", "Update color",width = 100)
              )
            )
          ),
          box(
            title =  "Set Plot Options",
            width = 8,
            status = "navy",
            solidHeader = T,
            pickerInput("selectgenelistoptions", "", width = 300, choices = distinct(gene_info,gene_list)$gene_list,
                        selected = distinct(gene_info,gene_list)$gene_list),
            pickerInput("selectdataoption", "", choices = distinct(gene_info,plot_legend)$plot_legend,
                        selected = distinct(gene_info,plot_legend)$plot_legend[1],
                        choicesOpt = list(
                          content = gsub("\\n", "\\1<br>", distinct(gene_info,plot_legend)$plot_legend)
                        ))
          ),
          modalButton("Close")
        )
      )))
  })
  
  # update display selected item info ----
  observeEvent(c(input$selectdataoption, input$selectgenelistoptions),
               ignoreInit = TRUE, ignoreNULL = TRUE,
               {
                 # print("options update")
                 my_sel <- LIST_DATA$meta_data %>% 
                   dplyr::filter(gene_list == input$selectgenelistoptions & 
                                   set == input$selectdataoption | plot_legend == input$selectdataoption)
                 
                 updateColourInput(session, "colourhex", value = paste(my_sel$mycol))
                 
                 updateTextInput(session,
                                 "textnickname",
                                 value = paste(my_sel$set))
               })
  
  # update color based on rgb text input ----
  observeEvent(input$actionmyrgb, {
    # print("color rgb")
    updateColourInput(session, "colourhex", value = RgbToHex(input$textrgbtohex, convert = "hex"))
  })
  
  # save color selected and update plot ----
  observeEvent(input$colourhex, ignoreInit = TRUE, {
    # print("update text color")
    updateTextInput(session,
                    "textrgbtohex",
                    value = RgbToHex(x = input$colourhex, convert = "rgb"))
    if (!is.null(names(LIST_DATA$gene_file))) {
      
      my_sel <- LIST_DATA$meta_data %>% 
        dplyr::filter(gene_list == input$selectgenelistoptions & 
                        set == input$selectdataoption | plot_legend == input$selectdataoption)
      ploton <- if_else(input$selectgenelistoptions == "Complete", T,
                        if_else(any(LIST_DATA$meta_data$gene_list == input$selectgenelistoptions & 
                                      LIST_DATA$meta_data$onoff !=0), T, F))
      
      if (input$colourhex != my_sel$mycol) {
        # print("color new")
        LIST_DATA$meta_data <<- LIST_DATA$meta_data %>% 
          dplyr::mutate(mycol=if_else((gene_list == input$selectgenelistoptions | 
                                         input$selectgenelistoptions == "Complete" |
                                         input$selectgenelistoptions == "Complete_filtered") & 
                                        (set == input$selectdataoption | plot_legend == input$selectdataoption),
                                      input$colourhex, mycol))
        reactive_values$Picker_controler <- 
          c(names(LIST_DATA$gene_file), distinct(LIST_DATA$meta_data, mycol)$mycol)
        if(ploton){
          if (!is.null(reactive_values$Apply_Math) & input$leftSideTabs == "mainplot") {
            reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA$meta_data)
          } 
        }
        if(input$leftSideTabs == "cdftool"){
          reactive_values$df_options <- my_sel
        }
      }
    }
  })
  
  # record new nickname  ----
  observeEvent(input$actionoptions, ignoreInit = TRUE, {
    # sets/resets nickname
    if (nchar(input$textnickname) == 0) {
      updateTextInput(session,
                      "textnickname",
                      value = paste(input$selectdataoption))
    } else if (input$textnickname != input$selectdataoption) {
      # print("new nickname")
      if (any(input$textnickname == distinct(LIST_DATA$meta_data, set)$set)) {
        updateTextInput(session,
                        "textnickname",
                        value = paste0(input$selectdataoption,"-",input$textnickname,
                                       "-dup"))
      }
      LIST_DATA$meta_data <<- LIST_DATA$meta_data %>% 
        dplyr::mutate(set = if_else(set == input$selectdataoption,
                                    input$textnickname, set)) %>% 
        dplyr::mutate(onoff = if_else(onoff == input$selectdataoption,
                                      input$textnickname, onoff)) %>% 
        dplyr::mutate(plot_legend = str_replace(LIST_DATA$meta_data$plot_legend,
                                                paste0("^",input$selectdataoption),input$textnickname)) %>% 
        dplyr::mutate(group = str_replace(LIST_DATA$meta_data$group,
                                          input$selectdataoption,input$textnickname))
      
      LIST_DATA$table_file <<- LIST_DATA$table_file %>%
        dplyr::mutate(set = if_else(set == input$selectdataoption,
                                    input$textnickname, set))
      reactive_values$Picker_controler <- 
        c(names(LIST_DATA$gene_file), distinct(LIST_DATA$meta_data, set)$set)
      ff <- distinct(LIST_DATA$table_file, set)$set
      updatePickerInput(session,
                        "selectdataoption",
                        choices = ff,
                        selected = input$textnickname,
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", ff)
                        ))
      
    }
  })
  
  # droplinesandlabels ----
  observeEvent(c(input$droplinesandlabels, reactive_values$droplinesandlabels), ignoreInit = T, ignoreNULL = T, {
    # print("droplinesandlabels")
    mynames <- LIST_DATA$meta_data_plot$tss_tes
    showModal(modalDialog(
      title = "Information message",
      " Set Lines and Labels for plot ",
      size = "xl",
      easyClose = F,
      footer = tagList(
        box(
          width = 12,
          solidHeader = F,
          collapsible = FALSE,
          collapsed = FALSE,
          column(12,
                 div(
                   helpText(paste("Upstream bp:",LIST_DATA$meta_data_plot$binning[3]))
                 ),
                 div(
                   helpText(paste("unscaled 5prime:",LIST_DATA$meta_data_plot$binning[6]))
                 ),
                 div(
                   helpText(paste("unscaled 3prime:",LIST_DATA$meta_data_plot$binning[7]))
                 ),
                 div(
                   helpText(paste("Downstream bp:",LIST_DATA$meta_data_plot$binning[4]))
                 ),
                 div(
                   helpText(paste("bin size:",LIST_DATA$meta_data_plot$binning[2]))
                 ),
          ),
          column(12,
                 div(
                   style = "padding:2px; display:inline-block; text-align:center;",
                   textInput("numerictssname", value = mynames[1], label = "TSS label",width = "100px"),
                 ),
                 div(
                   style = "padding:2px; display:inline-block; text-align:center;",
                   textInput("numerictesname", value = mynames[2], label = " TES label",width = "100px"),
                 ),
                 div(
                   style = "padding:2px 8px 2px 2px; display:inline-block; text-align:center;",
                   numericInput(
                     "numericlabelspaceing",
                     "label spacing",
                     value = LIST_DATA$meta_data_plot$binning[8],
                     min = 0,
                     max = 1000
                   )
                 )
          ),
          div(
            textInput("landlnames", reactive_values$Lines_Labels_List$mylabels, label = "Yaxis labels"),
            textInput("landlposition", reactive_values$Lines_Labels_List$mybrakes, label = "Yaxis label position (numbers only)")
          ),
          helpText("select buttons for more options"),
          column(
            12,
            div(
              style = "padding-left: 25px; display:inline-block;",
              dropdownButton(
                tags$h3("Set font Options"),
                numericInput(
                  inputId = 'selectlegendnewline',
                  "Set legend character brake",
                  value = reactive_values$Lines_Labels_List$mysize[7],
                  min = 0,
                  max = 50,
                  step = 1
                ),
                numericInput(
                  inputId = 'selectlegendnewlinespace',
                  "Set character brake window",
                  value = reactive_values$Lines_Labels_List$mysize[8],
                  min = 0,
                  max = 50,
                  step = 1
                ),
                numericInput(
                  inputId = 'selectlegendsize',
                  "Set legend size",
                  value = reactive_values$Lines_Labels_List$mysize[5],
                  min = 1,
                  max = 20,
                  step = 1
                ),
                numericInput(
                  inputId = 'selectfontsizex',
                  "Set X axis font size",
                  value = reactive_values$Lines_Labels_List$mysize[3],
                  min = 1,
                  max = 30,
                  step = 1
                ),
                numericInput(
                  inputId = 'selectfontsizey',
                  "Set Y axis font size",
                  value = reactive_values$Lines_Labels_List$mysize[4],
                  min = 1,
                  max = 30,
                  step = 1
                ),
                icon = icon("sliders"),
                status = "warning",
                tooltip = tooltipOptions(title = "Font Options")
              )
            ),
            div(
              style = "padding-left: 25px; display:inline-block;",
              dropdownButton(
                tags$h3("Set line Options"),
                
                numericInput(
                  inputId = 'selectlinesize',
                  "Set plot line size",
                  value = reactive_values$Lines_Labels_List$mysize[2],
                  min = .5,
                  max = 10,
                  step = .5
                ),
                numericInput(
                  inputId = 'selectvlinesize',
                  "Set vertcal line size",
                  value = reactive_values$Lines_Labels_List$mysize[1],
                  min = .5,
                  max = 10,
                  step = .5
                ),
                numericInput(
                  inputId = 'selectalpha',
                  "Set line alpha",
                  value = reactive_values$Lines_Labels_List$mysize[6],
                  min = .01,
                  max = 1,
                  step = .1
                ),
                icon = icon("sliders"),
                status = "warning",
                tooltip = tooltipOptions(title = "Line Options")
              )
            ),
            div(
              style = "padding-left: -5px; display:inline-block;",
              dropdownButton(
                tags$h3("Set TSS Options"),
                
                selectInput(
                  inputId = 'selecttsscolor',
                  label = 'TSS line and label color',
                  choices = c("green", "red", "blue", "brown", "black", "white"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_color[1]
                ),
                selectInput(
                  inputId = 'selecttssline',
                  label = 'TSS line type',
                  choices = c("dotted", "dashed", "solid"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_type[1]
                ),
                icon = icon("sliders"),
                status = "success",
                tooltip = tooltipOptions(title = "TSS Options")
              )
            ),
            div(
              style = "padding-left: 20px; display:inline-block;",
              dropdownButton(
                tags$h3("Set 5|4 Options"),
                
                selectInput(
                  inputId = 'selectbody1color',
                  label = '5|4 line and label color',
                  choices = c("black", "red", "green", "blue", "brown", "white"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_color[3]
                ),
                selectInput(
                  inputId = 'selectbody1line',
                  label = '5|4 line type',
                  choices = c("solid", "dotted"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_type[3]
                ),
                icon = icon("sliders"),
                tooltip = tooltipOptions(title = "5|4 Options")
              )
            ),
            div(
              style = "padding-left: 20px; display:inline-block;",
              dropdownButton(
                tags$h3("Set 4|3 Options"),
                
                selectInput(
                  inputId = 'selectbody2color',
                  label = '4|3 line and label color',
                  choices = c("black", "red", "green", "blue", "brown", "white"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_color[4]
                ),
                selectInput(
                  inputId = 'selectbody2line',
                  label = '4|3 line type',
                  choices = c("solid", "dotted"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_type[4]
                ),
                icon = icon("sliders"),
                tooltip = tooltipOptions(title = "4|3 Options")
              )
            ),
            div(
              style = "padding-left: 25px; display:inline-block;",
              dropdownButton(
                tags$h3("Set TES Options"),
                
                selectInput(
                  inputId = 'selecttescolor',
                  label = 'TES line and label color',
                  choices = c("red", "green", "blue", "brown", "black", "white"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_color[2]
                ),
                selectInput(
                  inputId = 'selecttesline',
                  label = 'TES line type',
                  choices = c("dotted", "dashed", "solid"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_type[2]
                ),
                icon = icon("sliders"),
                status = "danger",
                tooltip = tooltipOptions(title = "TES Options")
              )
            )
          )
        ),
        actionBttn("resetlineslabels", "reset"),
        actionBttn("actionlineslabels", "SET and Plot"),
      )
    ))
  })
  
  # observe switching tabs ----
  observeEvent(input$leftSideTabs, ignoreInit = TRUE, {
    # print("switch tab")
    # load files tab ----
    if (input$leftSideTabs == "loaddata"){
      if(length(names(LIST_DATA$gene_file))>1){
        shinyjs::enable("startoff2")
        updateSelectInput(session, "selectsave",
                          choices = names(LIST_DATA$gene_file)[-1],
                          selected = names(LIST_DATA$gene_file)[2])
      }
      
    }
    # main plot tab ----
    if (input$leftSideTabs == "mainplot") {
      reactive_values$Picker_controler <-
        c(names(LIST_DATA$gene_file), distinct(LIST_DATA$meta_data, set)$set)
      if(LIST_DATA$STATE[1] == 0){
        LIST_DATA$STATE[1] <<- 1
        LIST_DATA$STATE[2] <<- -10
        shinyjs::hide("actionmyplotshow")
        updateSliderTextInput(session,"sliderplotBinRange",
                              choices = reactive_values$slider_breaks$mylabels,
                              selected = c(first(reactive_values$slider_breaks$mylabels),
                                           last(reactive_values$slider_breaks$mylabels))
                              
        )
        updateSelectInput(session,"selectplotBinNorm",
                          choices = c("NA","min",
                                      reactive_values$slider_breaks$mylabels)
        )
        reactive_values$droplinesandlabels <- 1
        removeCssClass(selector = "a[data-value='grouptab']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='filenorm']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='genelists']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='sorttool']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='ratiotool']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='clustertool']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='groupiestool']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='cdftool']", class = "inactiveLink")
        removeCssClass(selector = "a[data-value='DataTableTool']", class = "inactiveLink")
      }
      if(all(LIST_DATA$meta_data$group == "self")){
        shinyjs::hide("hideplotgroup")
      } else {
        shinyjs::show("hideplotgroup")
      }
    }
    # sort/filter tab ----
    if (input$leftSideTabs == "sorttool"){
      ol <- input$sortGeneList
      if(!is.null(ol)){
        if (!ol %in% names(LIST_DATA$gene_file)) {
          ol <- "Complete"
        }
      }
      updateNumericInput(session, "peakfilternum",value = 0)
      output$rangeHelptext <- renderUI({helpText("select sample(s)")})
      updatePickerInput(session, "sortGeneList",
                        choices = names(LIST_DATA$gene_file),
                        selected = ol,
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", names(LIST_DATA$gene_file))
                        ))
      updatePickerInput(session, "sortSamples",
                        choices = c(distinct(LIST_DATA$meta_data, set)$set),
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$meta_data, set)$set)
                        )
      )
      output$valueboxsort <- renderValueBox({
        valueBox(0,
                 "Gene List Filter",
                 icon = icon("list"),
                 color = "green")
      })
      shinyjs::hide("hidesortplots1")
      shinyjs::hide("hidesortplots2")
      
    }
    # file group norm tab ----
    if (input$leftSideTabs == "grouptab") {
      updatePickerInput(
        session,
        "pickergroupsample",
        choices = distinct(LIST_DATA$meta_data, set)$set,
        choicesOpt = list(
          content = gsub("(.{55})", "\\1<br>", distinct(LIST_DATA$meta_data, set)$set)
        )
      )
      output$valueboxnormgroup <- renderValueBox({
        valueBox("0%",
                 "Done",
                 icon = icon("gears"),
                 color = "yellow")
      })
      gts <- LIST_DATA$meta_data %>% 
        filter(gene_list == "Complete") %>% 
        mutate(group = if_else(gsub("\n","", group) != set,as.character(as.integer(as.factor(group))),"self" )) %>% 
        select(set,group) %>%
        rename(File = set)
      dt <- datatable(
        gts,
        colnames = names(gts),
        rownames = FALSE,
        filter = "none",
        class = 'cell-border stripe compact',
        options = list(
          scrollX = TRUE,
          scrollY = TRUE,
          autoWidth = TRUE,
          width = 20,
          sDom  = '<"top">lrt<"bottom">ip',
          info = FALSE,
          paging = FALSE,
          lengthChange = FALSE,
          columnDefs = list(
            list(className = 'dt-center ', targets = "_all"),
            list(
              targets = 0,
              render = JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 25 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 27) + '...</span>' : data;",
                "}"
              )
            )
          )
        )
      )
      output$loadedfilestable2 <- DT::renderDataTable(dt)
    }
    # file norm tab ----
    if (input$leftSideTabs == "filenorm") {
      updatePickerInput(
        session,
        "pickernumerator",
        choices = distinct(LIST_DATA$meta_data, set)$set,
        choicesOpt = list(style = paste("color", dplyr::select(
          dplyr::filter(LIST_DATA$meta_data,
                        gene_list == names(LIST_DATA$gene_file)[1]),
          mycol)$mycol, sep = ":"),
          content = gsub("(.{55})", "\\1<br>", distinct(LIST_DATA$meta_data, set)$set)
        )
      )
      updatePickerInput(
        session,
        "pickerdenominator",
        choices = c(distinct(LIST_DATA$meta_data, set)$set,"Multiply by -1"),
        choicesOpt = list(style = paste("color", dplyr::select(
          dplyr::filter(LIST_DATA$meta_data,
                        gene_list == names(LIST_DATA$gene_file)[1]),
          mycol)$mycol, sep = ":"),
          content = gsub("(.{55})", "\\1<br>", c(distinct(LIST_DATA$meta_data, set)$set,"Multiply by -1"))
        )
      )
      output$valueboxnormfile <- renderValueBox({
        valueBox("0%",
                 "Done",
                 icon = icon("gears"),
                 color = "yellow")
      })
    }
    # genelists tab ----
    if (input$leftSideTabs == "genelists") {
      shinyjs::hide('actiongenelistsdatatable')
      updatePickerInput(
        session,
        "pickergenelists",
        choices = names(LIST_DATA$gene_file),
        choicesOpt = list(
          content = gsub("(.{55})", "\\1<br>", names(LIST_DATA$gene_file))
        )
      )
      output$valueboxgene1 <- renderValueBox({
        valueBox(0,
                 "Gene List innerjoin",
                 icon = icon("list"),
                 color = "green")
      })
      output$valueboxgene2 <- renderValueBox({
        valueBox(0,
                 "Gene List Total",
                 icon = icon("list"),
                 color = "yellow")
      })
      output$valueboxgene3 <- renderValueBox({
        valueBox(0,
                 "Gene List antijoin",
                 icon = icon("list"),
                 color = "red")
      })
    }
    # ratio switch tab ----
    if (input$leftSideTabs == "ratiotool"){
      if(!is.null(input$selectratiofile)){
        if(input$selectratiofile == "Load data file"){
          updateSliderTextInput(
            session,
            "sliderRatioBinNorm",
            choices = c("NA",
                        reactive_values$slider_breaks$mylabels),
            selected = "NA"
          )
        }
      } 
      updateSelectInput(
        session,
        "selectratiofile",
        choices = names(LIST_DATA$gene_file)
      )
      updatePickerInput(session, "pickerratio1file",
                        choices = c(distinct(LIST_DATA$meta_data, set)$set),
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$meta_data, set)$set)
                        ))
      updatePickerInput(session, "pickerratio2file",
                        choices = c("None", distinct(LIST_DATA$meta_data, set)$set),
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", c("None", distinct(LIST_DATA$meta_data, set)$set))
                        ))
      output$valueboxratio1 <- renderValueBox({
        valueBox(0,
                 "Ratio Up file1",
                 icon = icon("list"),
                 color = "green")
      })
      output$valueboxratio2 <- renderValueBox({
        valueBox(0,
                 "Ratio down file1",
                 icon = icon("list"),
                 color = "blue")
      })
      output$valueboxratio3 <- renderValueBox({
        valueBox(0,
                 "Ratio No Diff",
                 icon = icon("list"),
                 color = "yellow")
      })
    }
    # cluster switch tab ----
    if (input$leftSideTabs == "clustertool"){
      gl <- input$clusterSamples
      if(!is.null(gl)){ 
        if(gl[1] == "select sample"){
          shinyjs::hide("hideclusterplots1")
          shinyjs::hide("hideclustertable")
          gl <- c(distinct(LIST_DATA$meta_data, set)$set)[1]
        }
      }
      ol <- input$clusterGeneList
      if(!is.null(ol)){
        if (!ol %in% names(LIST_DATA$gene_file)) {
          ol <- "Complete"
        }
      }
      
      updatePickerInput(session, "clusterGeneList",
                        choices = names(LIST_DATA$gene_file),
                        selected = ol,
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", names(LIST_DATA$gene_file))
                        ))
      updatePickerInput(session, "clusterSamples",
                        choices = c(distinct(LIST_DATA$meta_data, set)$set),
                        selected = gl,
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$meta_data, set)$set)
                        )
      )
      output$valueboxsort <- renderValueBox({
        valueBox(0,
                 "Gene List Filter",
                 icon = icon("list"),
                 color = "green")
      })
    }
    
    # Groups switch tab ----
    if (input$leftSideTabs == "groupiestool"){
      gl <- input$groupiesSamples
      if(!is.null(gl)){
        if(gl[1] == "select sample"){
          shinyjs::hide("hidegroupiesplots1")
          shinyjs::hide("hidegroupiestable")
          shinyjs::hide("hidegroupiesplots2")
          gl <- c(distinct(LIST_DATA$meta_data, set)$set)[1]
        }
      }
      ol <- input$groupiesGeneList
      if(!is.null(ol)){
        if (!ol %in% names(LIST_DATA$gene_file)) {
          ol <- "Complete"
        }
      }
      
      updatePickerInput(session, "groupiesGeneList",
                        choices = names(LIST_DATA$gene_file),
                        selected = ol,
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", names(LIST_DATA$gene_file))
                        ))
      updatePickerInput(session, "groupiesSamples",
                        choices = c(distinct(LIST_DATA$meta_data, set)$set),
                        selected = gl,
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$meta_data, set)$set)
                        )
      )
      output$valueboxsort <- renderValueBox({
        valueBox(0,
                 "Gene List Filter",
                 icon = icon("list"),
                 color = "green")
      })
    }
    
    # CDF switch tab ----
    if(input$leftSideTabs == "cdftool"){
      shinyjs::hide('plotcdf')
      shinyjs::hide('plotcdfscatter')
      shinyjs::disable('actioncdfcolor')
      pickercdf <- list()
      for (i in names(LIST_DATA$gene_file)) {
        pickercdf[[i]] <-
          list(div(
            style = "margin-bottom: -10px;",
            pickerInput(
              inputId = paste0("-cdfspace3-", 
                               gsub(" ", "-cdfspace2-", 
                                    gsub("\n", "-cdfspace1-", i))),
              label = i,
              width = "99%",
              choices = distinct(LIST_DATA$meta_data, set)$set,
              multiple = T,
              options = list(
                `actions-box` = TRUE,
                `selected-text-format` = "count > 0"
              ),
              choicesOpt = list(style = paste("color", dplyr::select(
                dplyr::filter(LIST_DATA$meta_data,
                              gene_list == names(LIST_DATA$gene_file)[1]),
                mycol)$mycol, sep = ":"),
                content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$meta_data, set)$set))
            )
          ))
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Filter"))){
        output$DynamicCDFPicker_sort <- renderUI({
          pickercdf[str_detect(names(LIST_DATA$gene_file),"^Filter")]
        })
        shinyjs::show("showpickersort_cdf")
      } else {
        shinyjs::hide("showpickersort_cdf")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Gene_List_"))){
        output$DynamicCDFPicker_comparisons <- renderUI({
          pickercdf[str_detect(names(LIST_DATA$gene_file),"^Gene_List_")]
        })
        shinyjs::show("showpickercomparisons_cdf")
      } else {
        shinyjs::hide("showpickercomparisons_cdf")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Ratio_"))){
        output$DynamicCDFPicker_ratio <- renderUI({
          pickercdf[str_detect(names(LIST_DATA$gene_file),"^Ratio_")]
        })
        shinyjs::show("showpickerratio_cdf")
      } else {
        shinyjs::hide("showpickerratio_cdf")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Cluster_"))){
        output$DynamicCDFPicker_clusters <- renderUI({
          pickercdf[str_detect(names(LIST_DATA$gene_file),"^Cluster_")]
        })
        shinyjs::show("showpickercluster_cdf")
      } else {
        shinyjs::hide("showpickercluster_cdf")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Groups_"))){
        output$DynamicCDFPicker_groupies <- renderUI({
          pickercdf[str_detect(names(LIST_DATA$gene_file),"^Groups_")]
        })
        shinyjs::show("showpickergroupies_cdf")
      } else {
        shinyjs::hide("showpickergroupies_cdf")
      }
      output$DynamicCDFPicker_main <- renderUI({
        pickercdf[!str_detect(names(LIST_DATA$gene_file),"^Filter|^Gene_List_|^Ratio_|^Cluster_|^Groups_|^CDF")]
      })
      if (sum(grepl("CDF ", names(LIST_DATA$gene_file))) == 0) {
        output$plotcdf <- renderPlot({
          NULL
        })
        output$plotcdfscatter <- renderPlot({
          NULL
        })
        shinyjs::hide('plotcdf')
        shinyjs::hide('plotcdfscatter')
        shinyjs::hide('actioncdfdatatable')
        my_count <- 0
      } else {
        my_count <-
          n_distinct(LIST_DATA$gene_file[[grep("CDF ", names(LIST_DATA$gene_file))]]$full$gene, na.rm = T)
      }
    }
    # DataTable tab ----
    if(input$leftSideTabs == "DataTableTool"){
      
      ol <- input$pickerDT
      if(!is.null(ol)){
        if (!ol %in% names(LIST_DATA$gene_file)) {
          ol <- "Complete"
        }
      }
      updatePickerInput(session, "pickerDT",
                        choices = names(LIST_DATA$gene_file),
                        selected = ol,
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", names(LIST_DATA$gene_file))
                        ))
    }
    # QC tab ----
    if(input$leftSideTabs == "qcOptions"){
      shinyjs::hide("hidespinersQC")
      updatePickerInput(session, "QCsample",
                        choices = c(distinct(LIST_DATA$meta_data, set)$set),
                        selected = c(distinct(LIST_DATA$meta_data, set)$set)[1],
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", c(distinct(LIST_DATA$meta_data, set)$set))
                        ))
    }
  })
  
  # action reset button lines and labels ----
  observeEvent(input$resetlineslabels, ignoreInit = TRUE, ignoreNULL = T, {
    # print("action reset lines and labels")
    mynames <- LIST_DATA$meta_data_plot$tss_tes
    # if tss or tes location make sure there is text
    updateTextInput(session, "numerictssname", value = mynames[1])
    updateTextInput(session, "numerictesname", value = mynames[2])
    updateNumericInput(session,"numericlabelspaceing", value = input$numericlabelspaceing+1)
    updateNumericInput(session,"numericlabelspaceing", value = LIST_DATA$meta_data_plot$binning2[8])
    
  })
  
  # action button update lines and labels ----
  observeEvent(input$actionlineslabels, ignoreInit = TRUE, ignoreNULL = T, {
    # print("action lines and labels")
    my_pos <-
      suppressWarnings(as.numeric(unlist(
        strsplit(input$landlposition, split = "\\s+")
      )))
    my_label <- unlist(strsplit(input$landlnames, split = "\\s+"))
    if (length(my_pos) == 0) {
      my_label <- "none"
      my_pos <- LIST_DATA$meta_data_plot$x_plot_range[2] * 2
    }
    
    LIST_DATA$meta_data_plot$tss_tes <<- c(input$numerictssname,input$numerictesname)
    
    reactive_values$slider_breaks <- LinesLabelsSet(LIST_DATA$meta_data_plot$binning,
                                                    LIST_DATA$meta_data_plot$landmarks,
                                                    LIST_DATA$meta_data_plot$x_plot_range[2],
                                                    input$numerictssname,
                                                    input$numerictesname,
                                                    slider = T)
    reactive_values$setsliders <- SlidersSetsInfo(reactive_values$slider_breaks, 0)
    reactive_values$setsliders <- SlidersSetsInfo(reactive_values$slider_breaks, LIST_DATA$meta_data_plot$binning[1])
    
    reactive_values$slider_breaks$myselect  <- c(first(reactive_values$slider_breaks$mylabels),
                                                 last(reactive_values$slider_breaks$mylabels))
    updateSliderTextInput(session,"sliderplotBinRange",
                          choices = reactive_values$slider_breaks$mylabels,
                          selected = c(first(reactive_values$slider_breaks$mylabels),
                                       last(reactive_values$slider_breaks$mylabels))
                          
    )
    
    updateSelectInput(session,"selectplotBinNorm",
                      choices = c("NA","min",
                                  reactive_values$slider_breaks$mylabels),
                      selected = "NA"
    )
    
    reactive_values$Lines_Labels_List <-
      LinesLabelsPlot(
        LIST_DATA$meta_data_plot$binning,
        input$selectbody1color,
        input$selectbody1line,
        input$selectbody2color,
        input$selectbody2line,
        input$selecttsscolor,
        input$selecttssline,
        input$selecttescolor,
        input$selecttesline,
        my_label,
        my_pos,
        input$selectvlinesize,
        input$selectlinesize,
        input$selectfontsizex,
        input$selectfontsizey,
        input$selectlegendsize,
        input$selectlegendnewline,
        input$selectlegendnewlinespace,
        input$selectalpha,
        reactive_values$slider_breaks$lineloc
      )
    removeModal()
  })
  
  # update sliders ----
  observeEvent(reactive_values$setsliders, ignoreInit = TRUE, {
    # print("update sliders")
    updateSliderTextInput(
      session,
      "slidersortbinrange",
      choices = reactive_values$slider_breaks$mylabels,
      selected = reactive_values$setsliders[1:2]
      
    )
    updateSliderTextInput(
      session,
      "sliderbinratio1",
      choices = reactive_values$slider_breaks$mylabels,
      selected = reactive_values$setsliders[1:2]
    )
    updateSliderTextInput(
      session,
      "sliderbinratio2",
      choices = c("NA", reactive_values$slider_breaks$mylabels),
      selected = reactive_values$setsliders[3:4]
    )
    updateSliderTextInput(
      session,
      "sliderbincluster",
      choices = reactive_values$slider_breaks$mylabels,
      selected = reactive_values$setsliders[1:2]
    )
    updateSliderTextInput(
      session,
      "sliderbingroupies",
      choices = reactive_values$slider_breaks$mylabels,
      selected = reactive_values$setsliders[1:2]
    )
    updateSliderTextInput(
      session,
      "sliderbincdf1",
      choices = reactive_values$slider_breaks$mylabels,
      selected = reactive_values$setsliders[1:2]
    )
    updateSliderTextInput(
      session,
      "sliderbincdf2",
      choices = reactive_values$slider_breaks$mylabels,
      selected = reactive_values$setsliders[3:4]
    )
  })
  
  # keep sizes real numbers lines and labels ----
  observeEvent(
    c(
      input$selectvlinesize,
      input$selectlinesize,
      input$selectfontsizex,
      input$selectfontsizey,
      input$selectlegendsize,
      input$selectalpha,
      input$selectlegendnewline,
      input$selectlegendnewlinespace
    ),
    ignoreInit = TRUE,
    {
      if(LIST_DATA$STATE[2] > 0){
        # print("keep bin positions in bounds > 0")
        mynum <- c(2, 2.5, 13, 13, 10, 0.8, 20, 5)
        myset <- c(
          input$selectvlinesize,
          input$selectlinesize,
          input$selectfontsizex,
          input$selectfontsizey,
          input$selectlegendsize,
          input$selectalpha,
          input$selectlegendnewline,
          input$selectlegendnewlinespace
        )
        # keep bin positions in bounds > 0
        for (i in seq_along(myset)) {
          if (is.na(myset[i]) | myset[i] < 0) {
            myset[i] <- mynum[i]
            updateNumericInput(session, "selectvlinesize", value = myset[1])
            updateNumericInput(session, "selectlinesize", value = myset[2])
            updateNumericInput(session, "selectfontsizex", value = myset[3])
            updateNumericInput(session, "selectfontsizey", value = myset[4])
            updateNumericInput(session, "selectlegendsize", value = myset[5])
            updateNumericInput(session, "selectalpha", value = myset[6])
            updateNumericInput(session, "selectlegendnewline", value = myset[7])
            updateNumericInput(session, "selectlegendnewlinespace", value = myset[8])
          }
        }
      }
    })
  
  # Update lines and labels box's ----
  observeEvent(
    c(input$numerictssname,
      input$numerictesname,
      input$numericlabelspaceing
    ),
    ignoreInit = TRUE,
    ignoreNULL = TRUE,
    {
      # print("observe line and labels")
      myset <- input$numericlabelspaceing
      if (!is.na(myset) & (myset >= LIST_DATA$meta_data_plot$binning[2] &
                           as.double(myset)%%as.double(LIST_DATA$meta_data_plot$binning[2]) == 0) | myset == 1)  {
        shinyjs::enable("actionlineslabels")
        updateActionButton(session, "actionlineslabels", label = "SET and Plot")
        
        LIST_DATA$meta_data_plot$binning[8] <<- myset
        LIST_DATA$meta_data_plot$landmarks <<- LinesLableLandmarks(LIST_DATA$meta_data_plot$binning)
        Lines_Labels_List <- LinesLabelsSet(LIST_DATA$meta_data_plot$binning,
                                            LIST_DATA$meta_data_plot$landmarks,
                                            LIST_DATA$meta_data_plot$x_plot_range[2],
                                            input$numerictssname,
                                            input$numerictesname)
        
        # set label and position numbers
        updateTextInput(session,
                        "landlnames",
                        value = paste(Lines_Labels_List$mylabels, collapse = " "))
        updateTextInput(session,
                        "landlposition",
                        value = paste(Lines_Labels_List$mybrakes , collapse = " "))
        if(LIST_DATA$STATE[2] > 0) {
          LIST_DATA$STATE[2] <<- 1
        }
      } else {
        updateActionButton(session, "actionlineslabels", label = "bin size must be multiple of and not > other values")
        shinyjs::disable("actionlineslabels")
      }
    })
  
  # dropttest ----
  observeEvent(input$dropttest, ignoreInit = T, ignoreNULL = T, {
    if (is.null(reactive_values$ttest)){
      ttesttype <- "by files"
    } else {
      ttesttype <- reactive_values$ttest_values[1]
    }
    if(n_distinct(LIST_DATA$meta_data$gene_list, na.rm = T) > 1){
      reactive_values$ttest <- c("none","by files", "by lists")
    } else {
      reactive_values$ttest <- c("none","by files")
    }
    if(!is.null(LIST_DATA$ttest)){
      ttest_list <- c("select sample", distinct(LIST_DATA$ttest,set)$set)
    } else {
      ttest_list <- "select sample"
    }
    showModal(modalDialog(
      title = "Information message",
      " Update Nickname and color of samples",
      size = "xl",
      easyClose = F,
      footer = tagList(
        box(
          collapsed = F,
          collapsible = F,
          width = 12,
          status = "navy",
          solidHeader = T,
          title = "t.test settings",
          column(3,
                 pickerInput(inputId = "switchttest",
                             label = "Plot t.test",
                             choices = reactive_values$ttest,
                             selected = ttesttype
                 )
          ),
          column(3,
                 pickerInput(inputId = "switchttesttype",
                             label = "pick test",
                             choices = c("t.test","ks.test", "wilcox.test"),
                             selected = reactive_values$ttest_values[2])
          ),
          column(3,
                 pickerInput(inputId = "selectttestalt",
                             label = "alternative",
                             choices = c("two.sided", "less", "greater"),
                             selected = reactive_values$ttest_values[3])
          ),
          column(3,
                 pickerInput(inputId = "selectttestpaired",
                             label = "Paired",
                             choices = c("TRUE", "FALSE"),
                             selected = reactive_values$ttest_values[4])
          ),
          column(3,
                 pickerInput(inputId = "selectttestexact",
                             label = "exact",
                             choices = c("NULL", "TRUE", "FALSE"),
                             selected = reactive_values$ttest_values[5])
          ),
          column(3,
                 pickerInput(inputId = "selectttestlog",
                             label = "log p.value",
                             choices = c("none","-log", "-log10"),
                             selected = reactive_values$ttest_values[6])
          ),
          column(2,
                 pickerInput("padjust",
                             label = "p.adjust?",
                             choices = c("NO", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                         "fdr", "none"),
                             selected = reactive_values$ttest_values[7])
          )),
        box(
          collapsed = F,
          collapsible = F,
          width = 12,
          status = "navy",
          solidHeader = T,
          title = "t.test plot settings",
          column(
            3,
            numericInput("numericYRangeLowpval", label = "p.value Y min:", value = min(reactive_values$ttest_options[1:2]))
          ),
          column(
            3,
            numericInput("numericYRangeHighpval", label = "p.value Y max:", value = max(reactive_values$ttest_options[1:2]))
          ),
          column(3,
                 sliderInput(
                   "sliderplotOccupancy",
                   label = "p.value plot occupancy",
                   min = 1,
                   max = 3,
                   step = 0.5,
                   value = reactive_values$ttest_options[3])
          ),
          column(3,
                 numericInput("hlinettest",
                              "horizontal line p.val 0.05",
                              value = reactive_values$ttest_options[5],
                              min = -50,
                              max = 50,
                              step = .5)
          )
        ),
        box(
          collapsed = F,
          collapsible = F,
          width = 12,
          status = "navy",
          solidHeader = T,
          title = "t.test sample colors",
          pickerInput(inputId = "selectttestitem",
                      label = "Select to modify color",
                      choices = ttest_list,
                      selected = reactive_values$ttest_options[4]),
          colourInput("selectcolorttest", "Select color")
        ),
        actionButton("actionttest","Apply"),
        helpText("does not currently deal with grouped data using 'groups and single' option")
      )
    ))
  })
  
  # t.test gets line colors ----
  observeEvent(input$selectttestitem, ignoreInit = T,{
    if(!is.null(LIST_DATA$ttest)){
      if(input$selectttestitem != reactive_values$ttest_options[4]){
        # print("t.test gets line colors")
        mycol <- LIST_DATA$ttest %>% dplyr::filter(set == input$selectttestitem) %>% distinct(mycol)
        updateColourInput(session, "selectcolorttest", value = paste(mycol))
      }
    }
  })
  
  # t.test updates line colors ----
  observeEvent(input$selectcolorttest,ignoreInit = T,{
    if(!is.null(LIST_DATA$ttest)){
      # print("t.test updates line colors")
      LIST_DATA$ttest <<- LIST_DATA$ttest %>% dplyr::mutate(.,mycol = ifelse(set == input$selectttestitem, input$selectcolorttest, mycol))
      updatePickerInput(session, "selectttestitem", reactive_values$ttest_options[4])
    }
  })
  
  # t.test action button ----
  observeEvent(input$actionttest, ignoreInit = T, {
    ttest_values <- c(input$switchttest, input$switchttesttype, input$selectttestalt, 
                      input$selectttestpaired, input$selectttestexact, input$selectttestlog, 
                      input$padjust)
    reactive_values$ttest_options <- c(input$numericYRangeLowpval, input$numericYRangeHighpval, input$sliderplotOccupancy, 
                                       input$selectttestitem, input$hlinettest)
    
    if (LIST_DATA$STATE[2] != 2){
      # print("t.test button")
      list_data_frame <- Active_list_data(LIST_DATA,input$mygroup, input$checkboxfull, 
                                          input$selectlegendnewline, input$selectlegendnewlinespace)
      if (!is_empty(list_data_frame)) {
        if (sum(ttest_values == reactive_values$ttest_values) != 7){
          reactive_values$ttest_values <- ttest_values
          if (input$switchttest != "none"){
            if(input$switchttest == "by lists" & n_distinct(list_data_frame$gene_list, na.rm = T) == 1){
              showModal(modalDialog(
                title = "Information message",
                paste("Only 1 gene list active, replot with more than 1"),
                size = "s",
                easyClose = TRUE
              ))
              return()
            } else if (input$switchttest == "by files" & n_distinct(list_data_frame$set, na.rm = T) == 1) {
              showModal(modalDialog(
                title = "Information message",
                paste("Only 1 sample active, replot with more than 1"),
                size = "s",
                easyClose = TRUE
              ))
              return()
            }
            withProgress(message = 'Calculation in progress',
                         detail = 'This may take a while...',
                         value = 0,
                         {
                           LIST_DATA$ttest <<- ApplyTtest(list_data_frame,
                                                          input$switchttest,
                                                          input$selectttestlog,
                                                          input$switchttesttype,
                                                          input$padjust,
                                                          input$selectttestalt,
                                                          input$selectttestexact,
                                                          input$selectttestpaired,
                                                          input$mygroup,
                                                          input$selectlegendnewline, 
                                                          input$selectlegendnewlinespace)
                           mm <- YaxisValuetTest(LIST_DATA$ttest, input$hlinettest, input$selectttestlog )
                         })
          } else {
            LIST_DATA$ttest <<- NULL
            mm <- 0
          }
          reactive_values$ttest_options[1:2] <- mm
        }
        reactive_values$Plot_Options <- NULL
        reactive_values$Plot_Options <-
          MakePlotOptionFrame(LIST_DATA$meta_data)
      } else {
        LIST_DATA$STATE[2] <<- 2
        text = paste("Nothing selected to plot.\n")
        reactive_values$Plot_controler <- ggplot() +
          annotate(
            "text",
            x = 4,
            y = 25,
            size = 8,
            label = text
          ) +
          theme_void()
      }
    }
    removeModal()
  })
  
  # renders plots ----
  output$plot <- renderPlot({
    reactive_values$Plot_controler
  })
  # output$plot <- renderPlot({
  #   reactive_values$Plot_controler_Violin
  # })
  output$plot1sort <- renderPlot({
    reactive_values$Plot_controler_sort_min
  })
  output$plot2sort <- renderPlot({
    reactive_values$Plot_controler_sort_max
  })
  output$plot1cluster <- renderPlot({
    reactive_values$Plot_controler_cluster
  })
  output$plot1groupies <- renderPlot({
    reactive_values$Plot_controler_groupies
  })
  output$plot2groupies <- renderPlot({
    reactive_values$Plot_controler_dgroupies
  })
  output$plotratio <- renderPlot({
    reactive_values$Plot_controler_ratio
  })
  
  # reactive picker watcher (watches everything but grabs only gene lists with white space fix) ----
  observeEvent(reactiveValuesToList(input)[gsub(" ", "-bensspace2-", gsub("\n", "-bensspace1-", names(LIST_DATA$gene_file)))],
               ignoreNULL = FALSE,
               ignoreInit = TRUE,
               {
                 reactive_values$picker <-
                   reactiveValuesToList(input)[gsub(" ", "-bensspace2-", gsub("\n", "-bensspace1-", names(LIST_DATA$gene_file)))]
               })
  
  # records check box on/off ----
  observeEvent(reactive_values$picker,
               ignoreNULL = FALSE,
               ignoreInit = TRUE,
               {
                 # needed for controlling flow of first time auto plot
                 if (LIST_DATA$STATE[1] > .9) {
                   # print("checkbox on/off")
                   ttt <- reactive_values$picker
                   checkboxonoff <- list()
                   for(i in names(ttt)){
                     onoff_name <-
                       gsub("-bensspace2-", " ", gsub("-bensspace1-", "\n", i))
                     if(!is_empty(ttt[[i]])){
                       checkboxonoff <- bind_rows(checkboxonoff, tibble(gene_list = onoff_name,
                                                                        onoff = ttt[[i]], set = ttt[[i]]))
                     } else if(!all(is.na(names(ttt)))){
                       checkboxonoff <- bind_rows(checkboxonoff, tibble(gene_list = onoff_name,
                                                                        onoff = NA, set = NA))
                     }
                   }  
                   LIST_DATA$meta_data <<-
                     CheckBoxOnOff(checkboxonoff,
                                   LIST_DATA$meta_data)
                   if (LIST_DATA$STATE[2] > 0) {
                     shinyjs::show("actionmyplotshow")
                     LIST_DATA$STATE[2] <<- 2
                     LIST_DATA$ttest <<- NULL
                     reactive_values$Apply_Math <- NULL
                   } else if(LIST_DATA$STATE[2] == -10){
                     shinyjs::hide("actionmyplotshow")
                     LIST_DATA$STATE[2] <<- 1
                   }
                   # keeps plot button showing up unnecessarily 
                 } else if(LIST_DATA$STATE[2] > 0) {
                   LIST_DATA$STATE[1] <<- as.numeric(LIST_DATA$STATE[1]) + .25
                 }
               })
  
  # Dynamic Gene Pickers update ----
  observeEvent(
    reactive_values$Picker_controler,
    ignoreNULL = FALSE,
    ignoreInit = TRUE,
    {
      # print("Dynamic pickers update")
      
      pickerlist <- list()
      for (i in names(LIST_DATA$gene_file)) {
        pickerlist[[i]] <-
          list(div(
            style = "margin-bottom: -10px;",
            pickerInput(
              inputId = gsub(" ", "-bensspace2-", gsub("\n", "-bensspace1-", i)),
              label = i,
              width = "99%",
              choices = distinct(LIST_DATA$meta_data, set)$set,
              selected =  dplyr::select(dplyr::filter(LIST_DATA$meta_data, 
                                                      gene_list == i & onoff != 0), 
                                        onoff)$onoff,
              multiple = T,
              options = list(
                `actions-box` = TRUE,
                `selected-text-format` = "count > 0"),
              choicesOpt = list(style = paste("color", 
                                              dplyr::select(dplyr::filter(LIST_DATA$meta_data, 
                                                                          gene_list == i), mycol)$mycol,
                                              sep = ":"),
                                content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$meta_data, set)$set))
            )
          ))
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Filter"))){
        output$DynamicGenePicker_sort <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Filter")]
        })
        shinyjs::show("showpickersort")
      } else {
        shinyjs::hide("showpickersort")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Gene_List_"))){
        output$DynamicGenePicker_comparisons <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Gene_List_")]
        })
        shinyjs::show("showpickercomparisons")
      } else {
        shinyjs::hide("showpickercomparisons")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Ratio_"))){
        output$DynamicGenePicker_ratio <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Ratio_")]
        })
        shinyjs::show("showpickerratio")
      } else {
        shinyjs::hide("showpickerratio")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Cluster_"))){
        output$DynamicGenePicker_clusters <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Cluster_")]
        })
        shinyjs::show("showpickercluster")
      } else {
        shinyjs::hide("showpickercluster")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Groups_"))){
        output$DynamicGenePicker_groupies <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Groups_")]
        })
        shinyjs::show("showpickergroupies")
      } else {
        shinyjs::hide("showpickergroupies")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^CDF"))){
        output$DynamicGenePicker_cdf <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^CDF")]
        })
        shinyjs::show("showpickercdf")
      } else if(!any(str_detect(names(LIST_DATA$gene_file),"^CDF"))){
        shinyjs::hide("showpickercdf")
      }
      output$DynamicGenePicker_main <- renderUI({
        pickerlist[!str_detect(names(LIST_DATA$gene_file),"^Filter|^Gene_List_|^Ratio_|^Cluster_|^Groups_|^CDF")]
      })
      
    })
  
  

  # filter sum tool action ----
  observeEvent(input$actionsorttool, {
    # print("sort tool")
    if (input$slidersortpercent < 50 &
        input$selectsorttop == "Middle%") {
      updateSliderInput(session, "slidersortpercent", value = 50)
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   shinyjs::show("hidesortplots1")
                   reactive_values$Plot_controler_sort_min <- ggplot()
                   LD <- FilterTop(
                     LIST_DATA,
                     input$sortGeneList,
                     input$sortSamples,
                     floor(reactive_values$slider_breaks$mybrakes[
                       reactive_values$slider_breaks$mylabels %in% input$slidersortbinrange]),
                     input$slidersortbinrange,
                     input$slidersortpercent,
                     input$selectsorttop
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      LD <- LIST_DATA
      mylist <- last(grep("^Filter", names(LIST_DATA$gene_file)))
      LD$meta_data <- LD$meta_data %>%
        dplyr::mutate(onoff=if_else(gene_list == names(LD$gene_file)[mylist] &
                                      set %in% input$sortSamples, set, "0"))
      list_data_frame <- Active_list_data(LD, group="none", input$checkboxfull, 
                                          input$selectlegendnewline, input$selectlegendnewlinespace)
      if (!is_empty(list_data_frame)) {
        withProgress(message = 'Calculation in progress',
                     detail = 'This may take a while...',
                     value = 0,
                     {
                       Apply_Cluster_Math <-
                         ApplyMath(
                           list_data_frame,
                           "mean",
                           "none",
                           0,
                           0
                         )
                       reactive_values$Plot_controler_sort_min <- ggplot()
                       reactive_values$Plot_controler_sort_max <- ggplot()
                       gp1 <-
                         ggplot(Apply_Cluster_Math ,aes(as.numeric(bin),value,color=set)) +
                         geom_line() +
                         ylab("Mean bin value") +
                         theme(legend.position="bottom",
                               legend.title = element_blank(),
                               axis.title.x=element_blank())
                       print(gp1)
                       reactive_values$Plot_controler_sort_min <- gp1
                       shinyjs::show("hidesortplots1")
                       shinyjs::hide("hidesortplots2")
                     })
      }
      if (any(grep("^Filter", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxsort <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$full$gene, na.rm = T),
            "Gene List Filter",
            icon = icon("list"),
            color = "green"
          )
        })
      } else {
        output$valueboxsort <- renderValueBox({
          valueBox(0,
                   "Gene List Filter",
                   icon = icon("list"),
                   color = "green")
        })
      }
    } else {
      output$valueboxsort <- renderValueBox({
        valueBox(0,
                 "Gene List Filter",
                 icon = icon("list"),
                 color = "green")
      })
    }
    # updating select and keeping track if sort on sort
    ol <- input$sortGeneList
    if(!is.null(ol)){
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- "Complete"
      } 
    }
    updatePickerInput(session, "sortGeneList",
                      choices = names(LIST_DATA$gene_file),
                      selected = ol,
                      choicesOpt = list(
                        content = gsub("(.{35})", "\\1<br>", names(LIST_DATA$gene_file))
                      ))
    if(LIST_DATA$STATE[1] !=0 ){
      LIST_DATA$STATE[1] <<- 0.75
    }
  })
  
  
  
  observeEvent(input$actionaveragetool, {
    # print("average tool")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   shinyjs::show("hidesortplots1")
                   reactive_values$Plot_controler_sort_min <- ggplot()
                   LD <- FilterAverage(
                     LIST_DATA,
                     input$sortGeneList,
                     input$sortSamples,
                     floor(reactive_values$slider_breaks$mybrakes[
                       reactive_values$slider_breaks$mylabels %in% input$slidersortbinrange]),
                     input$slidersortbinrange,
                     input$selectaveragemath
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      LD <- LIST_DATA
      mylist <- grep("^Filter_all_bins", names(LIST_DATA$gene_file),value = T)
      LD$meta_data <- LD$meta_data %>%
        dplyr::mutate(onoff=if_else(gene_list %in% mylist &
                                      set %in% input$sortSamples, set, "0"))
      list_data_frame <- Active_list_data(LD, group="none", input$checkboxfull, 
                                          input$selectlegendnewline, input$selectlegendnewlinespace)
      if (!is_empty(list_data_frame)) {
        withProgress(message = 'Calculation in progress',
                     detail = 'This may take a while...',
                     value = 0,
                     {
                       Apply_Cluster_Math <-
                         ApplyMath(
                           list_data_frame,
                           "mean",
                           "none",
                           0,
                           0
                         )
                       reactive_values$Plot_controler_sort_min <- ggplot()
                       reactive_values$Plot_controler_sort_max <- ggplot()
                       gp1 <-
                         ggplot(Apply_Cluster_Math ,aes(as.numeric(bin),value,color=set)) +
                         geom_line() +
                         ylab("Mean bin value") +
                         theme(legend.position="bottom",
                               legend.title = element_blank(),
                               axis.title.x=element_blank())
                       print(gp1)
                       reactive_values$Plot_controler_sort_min <- gp1
                       shinyjs::show("hidesortplots1")
                       shinyjs::hide("hidesortplots2")
                     })
      }
      if (any(grep("^Filter all bins", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxsort <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[last(grep("^Filter all bins", names(LIST_DATA$gene_file)))]]$full$gene, na.rm = T),
            "All bins Filter",
            icon = icon("list"),
            color = "green"
          )
        })
      } else {
        output$valueboxsort <- renderValueBox({
          valueBox(0,
                   "All bins Filter",
                   icon = icon("list"),
                   color = "green")
        })
      }
    } else {
      output$valueboxsort <- renderValueBox({
        valueBox(0,
                 "All bins Filter",
                 icon = icon("list"),
                 color = "green")
      })
    }
    # updating select and keeping track if sort on sort
    ol <- input$sortGeneList
    if(!is.null(ol)){
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- "Complete"
      } 
    }
    updatePickerInput(session, "sortGeneList",
                      choices = names(LIST_DATA$gene_file),
                      selected = ol,
                      choicesOpt = list(
                        content = gsub("(.{35})", "\\1<br>", names(LIST_DATA$gene_file))
                      ))
    if(LIST_DATA$STATE[1] !=0 ){
      LIST_DATA$STATE[1] <<- 0.75
    }
  })
  
  # filter % numeric controller ----
  observeEvent(c(input$numericsortmin,input$numericsortmax, input$selectsortper), ignoreInit = TRUE, ignoreNULL = TRUE, {
    if (!is.numeric(input$numericsortmin)) {
      updateNumericInput(session, "numericsortmin", value = 1)
      return()
    }
    if (!is.numeric(input$numericsortmax)) {
      updateNumericInput(session, "numericsortmax", value = 99.5)
      return()
    }
    if (input$numericsortmin < 0 | input$numericsortmin > 100) {
      updateNumericInput(session, "numericsortmin", value = 1)
    }
    if (input$numericsortmax < 0 | input$numericsortmax > 100) {
      updateNumericInput(session, "numericsortmax", value = 99.5)
    }
    if(input$selectsortper == "between%"){
      if (input$numericsortmin > input$numericsortmax) {
        updateNumericInput(session, "numericsortmin", value = 1)
      }
      if (input$numericsortmax < input$numericsortmin) {
        updateNumericInput(session, "numericsortmax", value = 99.5)
      }
    }
  })
  
  # filter min max between % tool action ----
  observeEvent(input$actionsortper, ignoreInit = TRUE, {
    # print("sort % tool")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a bit ...',
                 value = 0,
                 {
                   shinyjs::show("hidesortplots1")
                   reactive_values$Plot_controler_sort_min <- ggplot()
                   sortmin <- FilterPer(LIST_DATA, 
                                        input$sortGeneList,
                                        input$sortSamples,
                                        floor(reactive_values$slider_breaks$mybrakes[
                                          reactive_values$slider_breaks$mylabels %in% input$slidersortbinrange]),
                                        c(input$numericsortmin,input$numericsortmax),
                                        input$selectsortper,
                                        input$slidersortbinrange)
                 })
    
    if(!is_empty(sortmin$sortplot)){
      LIST_DATA <<- sortmin
      if(input$selectsortper != "max%" ){  
        gp1 <-
          ggplot(sortmin$sortplot %>% dplyr::filter(!is.na(my_p_1)) ,aes(as.numeric(bin),my_p_1,color=set)) + 
          geom_line() + 
          ylab("Min % value") +
          theme(legend.position="bottom", 
                legend.title = element_blank(),
                axis.title.x=element_blank())
      } else{
        gp1 <- ggplot()
      }
      if (input$selectsortper != "min%" ){
        gp2 <-
          ggplot(sortmin$sortplot %>% dplyr::filter(!is.na(my_p_2)),aes(as.numeric(bin),my_p_2,color=set)) + 
          geom_line() + 
          ylab("Max % value") +
          theme(legend.position="bottom", 
                legend.title = element_blank(),
                axis.title.x=element_blank())
      } else {
        gp2 <- ggplot()
      }
      if (input$selectsortper == "min%" ){
        print(gp1)
        shinyjs::show("hidesortplots1")
        shinyjs::hide("hidesortplots2")
      } else if(input$selectsortper == "max%" ){
        print(gp2)
        shinyjs::show("hidesortplots2")
        shinyjs::hide("hidesortplots1")
      } else{
        print(gp1 + gp2 + plot_layout(ncol = 1))
        shinyjs::show("hidesortplots1")
        shinyjs::show("hidesortplots2")
      }
      reactive_values$Plot_controler_sort_min <- gp1
      reactive_values$Plot_controler_sort_max <- gp2
      output$valueboxsort <- renderValueBox({
        valueBox(
          n_distinct(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$full$gene, na.rm = T),
          "Gene List Filter",
          icon = icon("list"),
          color = "green"
        )
      })
      ol <- input$sortGeneList
      if(!is.null(ol)){
        if (!ol %in% names(LIST_DATA$gene_file)) {
          ol <- "Complete"
        } 
      }
      updatePickerInput(session, "sortGeneList",
                        choices = names(LIST_DATA$gene_file),
                        selected = ol,
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", names(LIST_DATA$gene_file))
                        ))
      ol <- last(grep("^Filter", names(LIST_DATA$gene_file), value = TRUE))
      updateNumericInput(session, "numericsortmin",
                         value = distinct(LIST_DATA$gene_file[[ol]]$full, min)$min)
      if(LIST_DATA$STATE[1] !=0 ){
        LIST_DATA$STATE[1] <<- 0.75
      }
      LIST_DATA$gene_file[[ol]]$full <<- LIST_DATA$gene_file[[ol]]$full %>% select(gene)
    } else {
      output$valueboxsort <- renderValueBox({
        valueBox(0,
                 "Gene List Filter",
                 icon = icon("list"),
                 color = "green")
      })
      reactive_values$Plot_controler_sort_min <-
        ggplot()
      reactive_values$Plot_controler_sort_max <-
        ggplot()
    }
  })
  
  # filter peak bin range info and suggestion ----
  observeEvent(c(input$sortSamples, input$slidersortbinrange), ignoreInit = TRUE, ignoreNULL = TRUE, {
    start_end_bin <- floor(reactive_values$slider_breaks$mybrakes[
      reactive_values$slider_breaks$mylabels %in% input$slidersortbinrange])
    if(length(start_end_bin)==1){
      start_end_bin <- c(start_end_bin,start_end_bin)
    }
    
    df <- LIST_DATA$table_file %>% 
      dplyr::filter(set %in% input$sortSamples) %>% 
      dplyr::filter(bin %in% start_end_bin[1]:start_end_bin[2])
    if(length(df$score)==0){
      rr <- c(0,0)
    } else {
      rr <- round(range(abs(df$score),na.rm = T),digits = 4)
    }
    
    updateNumericInput(session, "peakfilternum", value = max(rr,na.rm = T))
    output$rangeHelptext <- renderUI({helpText(paste("min",min(rr,na.rm = T),"max",max(rr,na.rm = T)))})
  })
  
  # filter peak tool action ----
  observeEvent(input$actionsortpeak, ignoreInit = TRUE, {
    # print("sort Peak")
    shinyjs::show("hidesortplots1")
    reactive_values$Plot_controler_sort_min <- ggplot()
    shinyjs::hide("hidesortplots2")
    if (is.null(input$sortSamples)) {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return()
    }
    if(!is.numeric(input$peakfilternum)){
      updateNumericInput(session, "peakfilternum", value = 1)
      showModal(modalDialog(
        title = "Information message",
        paste("please set score to filter on"),
        size = "s",
        easyClose = TRUE
      ))
      return()
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a bit ...',
                 value = 0,
                 {
                   sortmin <- FilterPeak(LIST_DATA, 
                                         input$sortGeneList,
                                         input$sortSamples,
                                         floor(reactive_values$slider_breaks$mybrakes[
                                           reactive_values$slider_breaks$mylabels %in% input$slidersortbinrange]),
                                         input$selectsortpeak,
                                         input$slidersortbinrange,
                                         input$peakfilternum)
                 })
    
    if(!is.null(sortmin)){
      LIST_DATA <<- sortmin
      # pull info for preview plot
      mylist <- last(grep("^Filter", names(sortmin$gene_file),value = T))
      sortmin$meta_data <- sortmin$meta_data %>%
        dplyr::mutate(onoff=if_else(gene_list == mylist &
                                      set %in% input$sortSamples, set, "0"))
      list_data_frame <- Active_list_data(sortmin, group="none", input$checkboxfull, 
                                          input$selectlegendnewline, input$selectlegendnewlinespace)
      if (!is_empty(list_data_frame)) {
        reactive_values$Plot_controler_sort_min <- ggplot()
        withProgress(message = 'Calculation in progress',
                     detail = 'This may take a while...',
                     value = 0,
                     {
                       Apply_Math <-
                         ApplyMath(
                           list_data_frame,
                           "mean",
                           "none",
                           0,
                           0
                         )
                     })
        gp1 <-
          ggplot(Apply_Math ,aes(as.numeric(bin),value,color=set)) +
          geom_line() +
          ylab("Filtered") +
          theme(legend.position="bottom",
                legend.title = element_blank(),
                axis.title.x=element_blank())
        print(gp1)
        reactive_values$Plot_controler_sort_min <- gp1
        shinyjs::show("hidesortplots1")
        output$valueboxsort <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$full$gene, na.rm = T),
            "Gene List Filter",
            icon = icon("list"),
            color = "green"
          )
        })
        ol <- input$sortGeneList
        if(!is.null(ol)){
          if (!ol %in% names(LIST_DATA$gene_file)) {
            ol <- "Complete"
          } 
        }
        updatePickerInput(session, "sortGeneList",
                          choices = names(LIST_DATA$gene_file),
                          selected = ol,
                          choicesOpt = list(
                            content = gsub("(.{35})", "\\1<br>", names(LIST_DATA$gene_file))
                          ))
        if(LIST_DATA$STATE[1] != 0 ){
          LIST_DATA$STATE[1] <<- 0.75
        }
      }
    } else {
      output$valueboxsort <- renderValueBox({
        valueBox(0,
                 "Gene List Filter",
                 icon = icon("list"),
                 color = "green")
      })
      showModal(modalDialog(
        title = "Information message",
        paste("no bins in non-peak range that are bigger than max in peak range"),
        size = "s",
        easyClose = TRUE
      ))
    }
  })
  
  # filters size and separation ----
  observeEvent(input$actionSizeSep,{
    genesep <- floor(input$geneSeparation)
    my_name <- "Complete_filtered"
    mysub <- "Complete_filtered:"
    if(!is.numeric(genesep) || genesep < 0){
      updateNumericInput(session, "geneSizeMin", value = 0)
      genesep <- 0
    } else {
      mysub <- paste(mysub, "gene sep >", genesep)
    }
    genesizemin <- floor(input$geneSizeMin)
    if(!is.numeric(genesizemin) || genesizemin <= 0){
      updateNumericInput(session, "geneSizeMin", value = 0)
      genesizemin <- 0
    } else {
      mysub <- paste(mysub, "gene length >", genesizemin)
    }
    genesizemax <- floor(input$geneSizeMax)
    if(!is.numeric(genesizemax) || genesizemax <= 0){
      updateNumericInput(session, "geneSizeMMax", value = 0)
      genesizemax <- 0
    }  else {
      mysub <- paste(mysub, "gene length <", genesizemax)
    }
    LIST_DATA$gene_file[[my_name]] <<- NULL
    LIST_DATA$meta_data <<- LIST_DATA$meta_data %>% filter(gene_list != my_name)
    
    LD <- FilterSepSize(distinct(LIST_DATA$table_file,gene,chrom,start,end,strand),
                        genesep,
                        genesizemin,
                        genesizemax,
                        input$checkboxStranded)
    if(n_distinct(LD$gene) > 0) {
      # adds full n count to nickname
      
      # preps meta data
      gene_info <- LIST_DATA$meta_data %>% 
        dplyr::filter(gene_list == "Complete") %>% 
        dplyr::mutate(gene_list = my_name, 
                      count = paste("n =", n_distinct(LD$gene, na.rm = T)),
                      sub = mysub, 
                      onoff = "0",
                      plot_legend = " ")
      # saves data in list of lists
      LIST_DATA$gene_file[[my_name]]$full <<- distinct(LD)
      LIST_DATA$gene_file[[my_name]]$info <<- tibble(loaded_info =
                                                       paste("Loaded gene list from file",
                                                             my_name,
                                                             Sys.Date()),
                                                     save_name = gsub(" ", "_", str_squish(paste(my_name,
                                                                                                 Sys.Date(), sep ="_"))),
                                                     col_info = "loaded file"
      )
      LIST_DATA$meta_data <<- bind_rows(LIST_DATA$meta_data, gene_info)
      output$valueboxsort <- renderValueBox({
        valueBox(
          n_distinct(LIST_DATA$gene_file[[my_name]]$full$gene, na.rm = T),
          "Gene List Filter",
          icon = icon("list"),
          color = "green"
        )
      })
    } else {
      output$valueboxsort <- renderValueBox({
        valueBox(
          0,
          "Gene List Filter",
          icon = icon("list"),
          color = "green"
        )
      })
    }
    
    
    ol <- input$sortGeneList
    if(!is.null(ol)){
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- "Complete"
      } 
    }
    updatePickerInput(session, "sortGeneList",
                      choices = names(LIST_DATA$gene_file),
                      selected = ol,
                      choicesOpt = list(
                        content = gsub("(.{35})", "\\1<br>", names(LIST_DATA$gene_file))
                      ))
    if(LIST_DATA$STATE[1] != 0 ){
      LIST_DATA$STATE[1] <<- 0.75
    }
    
  })
  
  # observe norm gene pickers ----
  observeEvent(c(input$pickernumerator, input$adddata,
                 input$pickerdenominator), {
                   if (input$pickernumerator != "") {
                     if(input$pickerdenominator != "Multiply by -1"){
                       updateTextInput(session, "textnromname",value = paste(input$pickernumerator, input$adddata,
                                                                             input$pickerdenominator))
                     } else {
                       updateTextInput(session, "textnromname",value = paste0(input$pickernumerator))
                     }
                     output$valueboxnormfile <- renderValueBox({
                       valueBox("0%",
                                "Done",
                                icon = icon("gears"),
                                color = "yellow")
                     })
                   }
                 })
  
  # create group file ----
  observeEvent(input$actionnormgroup, ignoreInit = TRUE, {
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- MakeGroupFile(
                     LIST_DATA,
                     input$pickergroupmath
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      
      updatePickerInput(
        session,
        "pickergroupsample",
        choices = distinct(LIST_DATA$meta_data, set)$set,
        choicesOpt = list(
          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$meta_data, set)$set)
        )
      )
      gts <- LIST_DATA$meta_data %>% 
        filter(gene_list == "Complete") %>% 
        mutate(group = if_else(gsub("\n","", group) != set,as.character(as.integer(as.factor(group))),"self" )) %>% 
        select(set,group) %>%
        rename(File = set)
      dt <- datatable(
        gts,
        colnames = names(gts),
        rownames = FALSE,
        filter = "none",
        class = 'cell-border stripe compact',
        options = list(
          scrollX = TRUE,
          scrollY = TRUE,
          autoWidth = TRUE,
          width = 20,
          sDom  = '<"top">lrt<"bottom">ip',
          info = FALSE,
          paging = FALSE,
          lengthChange = FALSE,
          columnDefs = list(
            list(className = 'dt-center ', targets = "_all"),
            list(
              targets = 0,
              render = JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 35 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 37) + '...</span>' : data;",
                "}"
              )
            )
          )
        )
      )
      output$loadedfilestable2 <- DT::renderDataTable(dt)
      updateTextInput(session, "textnromname", value = "")
      output$valueboxnormgroup <- renderValueBox({
        valueBox(
          "Done",
          paste("Compleat n =", n_distinct(LIST_DATA$gene_file[[1]]$full$gene, na.rm = T)),
          icon = icon("thumbs-up", lib = "glyphicon"),
          color = "green"
        )
      })
      ff <- distinct(LIST_DATA$meta_data, set)$set
      updateSelectInput(session,
                        "selectdataoption",
                        choices = ff)
    } else {
      #no new data file created
      output$valueboxnormgroup <- renderValueBox({
        valueBox("0%",
                 "Done",
                 icon = icon("gears"),
                 color = "red")
      })
    }
  })
  
  # create norm file ----
  observeEvent(input$actionnorm, ignoreInit = TRUE, {
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- MakeNormFile(
                     LIST_DATA,
                     input$pickernumerator,
                     input$pickerdenominator,
                     input$radiogenebygene,
                     input$checkboxnormzero,
                     input$checkboxnormzero2,
                     input$adddata,
                     input$textnromname
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      updatePickerInput(session,
                        "pickernumerator", selected = "",
                        choices = distinct(LIST_DATA$meta_data, set)$set,
                        choicesOpt = list(style = paste("color", 
                                                        dplyr::select(
                                                          dplyr::filter(LIST_DATA$meta_data,
                                                                        gene_list == names(
                                                                          LIST_DATA$gene_file)[1]), 
                                                          mycol)$mycol, 
                                                        sep = ":"),
                                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$meta_data, set)$set)
                        ))
      updatePickerInput(session,
                        "pickerdenominator", selected = "",
                        choices = c(distinct(LIST_DATA$meta_data, set)$set,"Multiply by -1"),
                        choicesOpt = list(style = paste("color", 
                                                        dplyr::select(
                                                          dplyr::filter(LIST_DATA$meta_data,
                                                                        gene_list == names(
                                                                          LIST_DATA$gene_file)[1]),
                                                          mycol)$mycol, 
                                                        sep = ":"),
                                          content = gsub("(.{35})", "\\1<br>", c(distinct(LIST_DATA$meta_data, set)$set,"Multiply by -1"))
                        ))
      updatePickerInput(
        session,
        "pickergroupsample",
        choices = distinct(LIST_DATA$meta_data, set)$set,
        choicesOpt = list(
          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$meta_data, set)$set)
        )
      )
      gts <- LIST_DATA$meta_data %>% 
        filter(gene_list == "Complete") %>% 
        mutate(group = if_else(gsub("\n","", group) != set,as.character(as.integer(as.factor(group))),"self" )) %>% 
        select(set,group) %>%
        rename(File = set)
      dt <- datatable(
        gts,
        colnames = names(gts),
        rownames = FALSE,
        filter = "none",
        class = 'cell-border stripe compact',
        options = list(
          scrollX = TRUE,
          scrollY = TRUE,
          autoWidth = TRUE,
          width = 20,
          sDom  = '<"top">lrt<"bottom">ip',
          info = FALSE,
          paging = FALSE,
          lengthChange = FALSE,
          columnDefs = list(
            list(className = 'dt-center ', targets = "_all"),
            list(
              targets = 0,
              render = JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 35 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 37) + '...</span>' : data;",
                "}"
              )
            )
          )
        )
      )
      output$loadedfilestable2 <- DT::renderDataTable(dt)
      updateTextInput(session, "textnromname", value = "")
      output$valueboxnormfile <- renderValueBox({
        valueBox(
          "Done",
          paste("Compleat n =", n_distinct(LIST_DATA$gene_file[[1]]$full$gene, na.rm = T)),
          icon = icon("thumbs-up", lib = "glyphicon"),
          color = "green"
        )
      })
      ff <- distinct(LIST_DATA$meta_data, set)$set
      updateSelectInput(session,
                        "selectdataoption",
                        choices = ff)
    } else {
      #no new data file created
      output$valueboxnormfile <- renderValueBox({
        valueBox("0%",
                 "Done",
                 icon = icon("gears"),
                 color = "red")
      })
    }
  })
  
  # observe norm gene pickers ----
  observeEvent(input$pickergroupsample, ignoreNULL = T,ignoreInit = T,{
    updateTextInput(session, "textgroupname",
                    value = input$pickergroupsample)
  })
  
  # create groups file ----
  observeEvent(input$actiongroup, ignoreInit = TRUE, {
    if (!is.null(input$pickergroupsample)) {
      LIST_DATA$meta_data <<- LIST_DATA$meta_data %>% 
        mutate(group = if_else(set %in% input$pickergroupsample,
                               as_tibble(insert_line_breaks(gsub(",",":",input$textgroupname), 
                                                            input$selectlegendnewline, 
                                                            input$selectlegendnewlinespace))$value, group))
      updatePickerInput(session,
                        "pickergroupsample", selected = "",
                        options = list(title = "Select at least 2 files"))
      updateTextInput(session, "textgroupname", value = "")
      gts <- LIST_DATA$meta_data %>% 
        filter(gene_list == "Complete") %>% 
        mutate(group = if_else(gsub("\n","", group) != set,as.character(as.integer(as.factor(group))),"self" )) %>% 
        select(set,group) %>%
        rename(File = set)
      dt <- datatable(
        gts,
        colnames = names(gts),
        rownames = FALSE,
        filter = "none",
        class = 'cell-border stripe compact',
        options = list(
          scrollX = TRUE,
          scrollY = TRUE,
          autoWidth = TRUE,
          width = 20,
          sDom  = '<"top">lrt<"bottom">ip',
          info = FALSE,
          paging = FALSE,
          lengthChange = FALSE,
          columnDefs = list(
            list(className = 'dt-center ', targets = "_all"),
            list(
              targets = 0,
              render = JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 35 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 37) + '...</span>' : data;",
                "}"
              )
            )
          )
        )
      )
      output$loadedfilestable2 <- DT::renderDataTable(dt)
      shinyjs::show("hideplotgroup")
    }
  })
  
  # Gene action ----
  observeEvent(input$actiongenelists, {
    # print("gene lists action")
    shinyjs::hide('actiongenelistsdatatable')
    shinyjs::hide('genelists1table')
    shinyjs::hide('genelists2table')
    shinyjs::hide('genelists3table')
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- IntersectGeneLists(LIST_DATA,
                                            input$pickergenelists)
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      
      ol <- input$pickergenelists
      if(!is.null(ol)){
        if (!any(ol %in% names(LIST_DATA$gene_file))) {
          ol <- grep("Gene_List_", names(LIST_DATA$gene_file), value = TRUE)
        }
      }
      updatePickerInput(
        session,
        "pickergenelists",
        choices = names(LIST_DATA$gene_file),
        selected = ol,
        choicesOpt = list(
          content = gsub("(.{55})", "\\1<br>", names(LIST_DATA$gene_file))
        )
      )
      shinyjs::show('actiongenelistsdatatable')
      if (any(grep("Gene_List_innerjoin\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxgene1 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Gene_List_innerjoin\nn =",
                                                 names(LIST_DATA$gene_file))]]$full$gene, na.rm = T),
            "Gene List innerjoin",
            icon = icon("list"),
            color = "green"
          )
        })
      } else{
        output$valueboxgene1 <- renderValueBox({
          valueBox(0,
                   "Gene List innerjoin",
                   icon = icon("list"),
                   color = "green")
        })
      }
      if (any(grep("Gene_List_Total\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxgene2 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Gene_List_Total\nn =", names(LIST_DATA$gene_file))]]$full$gene, na.rm = T),
            "Gene List Total",
            icon = icon("list"),
            color = "yellow"
          )
        })
      } else{
        output$valueboxgene2 <- renderValueBox({
          valueBox(0,
                   "Gene List Total",
                   icon = icon("list"),
                   color = "yellow")
        })
      }
      if (any(grep("Gene_List_antijoin\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxgene3 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Gene_List_antijoin\nn =",
                                                 names(LIST_DATA$gene_file))]]$full$gene, na.rm = T),
            "Gene List antijoin",
            icon = icon("list"),
            color = "red"
          )
        })
      } else{
        output$valueboxgene3 <- renderValueBox({
          valueBox(0,
                   "Gene List antijoin",
                   icon = icon("list"),
                   color = "red")
        })
      }
    } else {
      output$valueboxgene1 <- renderValueBox({
        valueBox(0,
                 "Gene List innerjoin",
                 icon = icon("list"),
                 color = "green")
      })
      output$valueboxgene2 <- renderValueBox({
        valueBox(0,
                 "Gene List Total",
                 icon = icon("list"),
                 color = "yellow")
      })
      output$valueboxgene3 <- renderValueBox({
        valueBox(0,
                 "Gene List antijoin",
                 icon = icon("list"),
                 color = "red")
      })
      return()
    }
  })
  
  # Gene lists DT show gene list ----
  observeEvent(input$actiongenelistsdatatable, ignoreInit = TRUE, {
    # print("generate gene lists table")
    shinyjs::hide('actiongenelistsdatatable')
    if (any(grep("Gene_List_innerjoin\nn =", names(LIST_DATA$gene_file)) >
            0)) {
      newnames1 <-
        gsub("\n",
             " ",
             grep(
               "Gene_List_innerjoin\nn =",
               names(LIST_DATA$gene_file),
               value = TRUE
             ))
      mytab <- "innerjoined Gene Lists"
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$genelists1table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep("Gene_List_innerjoin\nn =",
                                                     names(LIST_DATA$gene_file))]]$full,
                           rownames = FALSE,
                           colnames = newnames1,
                           class = 'cell-border stripe compact',
                           filter = "none",
                           caption = LIST_DATA$gene_file[[grep("Gene_List_innerjoin\nn =",
                                                               names(LIST_DATA$gene_file))]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
                             columnDefs = list(
                               list(className = 'dt-center ', targets = "_all"),
                               list(
                                 targets = 0,
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data.length > 44 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 39) + '...</span>' : data;",
                                   "}"
                                 )
                               )
                             )
                           )
                         )
                       )
                   })
      shinyjs::show('genelists1table')
    } else {
      output$genelists1table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = "Gene_List_antijoin n = 0",
            options = list(searching = FALSE)
          )
        )
      mytab <- "Total Gene Lists"
    }
    if (any(grep("Gene_List_Total\nn =", names(LIST_DATA$gene_file)) >
            0)) {
      newnames2 <-
        gsub("\n",
             " ",
             grep(
               "Gene_List_Total\nn =",
               names(LIST_DATA$gene_file),
               value = TRUE
             ))
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$genelists2table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep("Gene_List_Total\nn =",
                                                     names(LIST_DATA$gene_file))]]$full,
                           rownames = FALSE,
                           colnames = newnames2,
                           class = 'cell-border stripe compact',
                           filter = "none",
                           caption = LIST_DATA$gene_file[[grep("Gene_List_Total\nn =",
                                                               names(LIST_DATA$gene_file))]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
                             columnDefs = list(
                               list(className = 'dt-center ', targets = "_all"),
                               list(
                                 targets = 0,
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data.length > 44 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 39) + '...</span>' : data;",
                                   "}"
                                 )
                               )
                             )
                           )
                         )
                       )
                   })
      shinyjs::show('genelists2table')
    } else {
      output$genelists2table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = "Gene_List_antijoin n = 0",
            options = list(searching = FALSE)
          )
        )
      if (mytab == "Total Gene Lists") {
        mytab <- "antijoin Gene Lists"
      }
    }
    if (any(grep("Gene_List_antijoin\nn =", names(LIST_DATA$gene_file)) >
            0)) {
      newnames3 <-
        gsub("\n",
             " ",
             grep(
               "Gene_List_antijoin\nn =",
               names(LIST_DATA$gene_file),
               value = T
             ))
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$genelists3table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep("Gene_List_antijoin\nn =",
                                                     names(LIST_DATA$gene_file))]]$full,
                           rownames = FALSE,
                           colnames = newnames3,
                           class = 'cell-border stripe compact',
                           filter = "none",
                           caption = LIST_DATA$gene_file[[grep("Gene_List_antijoin\nn =",
                                                               names(LIST_DATA$gene_file))]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
                             columnDefs = list(
                               list(className = 'dt-center ', targets = "_all"),
                               list(
                                 targets = 0,
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data.length > 44 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 39) + '...</span>' : data;",
                                   "}"
                                 )
                               )
                             )
                           )
                         )
                       )
                   })
      shinyjs::show('genelists3table')
    } else {
      output$genelists3table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = "Gene_List_antijoin n = 0",
            options = list(searching = FALSE)
          )
        )
      if (mytab == "antijoin Gene Lists") {
        mytab <- "Total Gene Lists"
      }
    }
    updateTabItems(session, "geneliststooltab", mytab)
  })
  
  # Cluster tool action ----
  observeEvent(input$actionclustertool, ignoreInit = TRUE, {
    # print("cluster tool action")
    shinyjs::hide('plot1cluster')
    if (n_distinct(LIST_DATA$gene_file[[input$clusterGeneList]]$full, na.rm = T) < as.numeric(input$selectclusternumber) |
        is.null(input$clusterSamples)) {
      showModal(modalDialog(
        title = "Information message",
        paste("Can't make more clusters than number of genes"),
        size = "s",
        easyClose = TRUE
      ))
      return()
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <-
                     FindClusters(
                       LIST_DATA,
                       input$clusterGeneList,
                       input$clusterSamples,
                       floor(reactive_values$slider_breaks$mybrakes[
                         reactive_values$slider_breaks$mylabels %in% input$sliderbincluster]),
                       input$clustpattern
                     )
                 })
    reactive_values$clustergroups <- NULL
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      if(LIST_DATA$STATE[1] != 0){
        LIST_DATA$STATE[1] <<- .5
      }
      reactive_values$clustergroups <- str_detect(names(LIST_DATA$gene_file),"^Cluster_")
    } else {
      return()
    }
  })
  
  
  # Plots Clusters based on number selected ----
  observeEvent(c(input$selectclusternumber, reactive_values$clustergroups, input$clusterRF),
               ignoreInit = TRUE, ignoreNULL = TRUE,
               {
                 # print("cluster tool number")
                 if (is.null(reactive_values$clustergroups)) {
                   return()
                 }
                 shinyjs::hide('hideclusterplots1')
                 shinyjs::hide("hideclustertable")
                 withProgress(message = 'Calculation in progress',
                              detail = 'This may take a while...',
                              value = 0,
                              {
                                LD <-
                                  ClusterNumList(
                                    LIST_DATA,
                                    input$clusterGeneList,
                                    input$clusterSamples,
                                    input$sliderbincluster,
                                    input$selectclusternumber
                                  )
                              })
                 if (!is_empty(LD$table_file)) {
                   LIST_DATA <<- LD
                   ol <- input$selectclusterfile
                   shinyjs::show('hideclusterplots1')
                   shinyjs::show("hideclustertable")
                   shinyjs::show('plot1cluster')
                   updateSelectInput(
                     session,
                     "selectclusterfile",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   LD$meta_data <- LD$meta_data %>%
                     dplyr::mutate(onoff=if_else(str_detect(gene_list,"^Cluster_") &
                                                   set == input$clusterSamples, set, "0"))
                   withProgress(message = 'Calculation in progress',
                                detail = 'This may take a while...',
                                value = 0,
                                {
                                  list_data_frame <- Active_list_data(LD, group="none", input$checkboxfull, 
                                                                      input$selectlegendnewline, input$selectlegendnewlinespace)
                                  if (!is_empty(list_data_frame)) {
                                    if(input$clusterRF){
                                      clusterRF <- "relative frequency"
                                      ylabRF <- "relative frequency"
                                    } else {
                                      clusterRF <- "none"
                                      ylabRF <- "mean of bin counts"
                                    }
                                    Apply_Cluster_Math <- ApplyMath(
                                      list_data_frame,
                                      relative_frequency=clusterRF
                                    )
                                  }
                                  reactive_values$Plot_controler_cluster <- ggplot()
                                  gp1 <-
                                    ggplot(Apply_Cluster_Math ,aes(as.numeric(bin),value,color=gene_list)) +
                                    geom_line(linewidth=1) +
                                    ylab(ylabRF) +
                                    theme(legend.position="bottom",
                                          legend.title = element_blank(),
                                          axis.title.x=element_blank())
                                  print(gp1)
                                  reactive_values$Plot_controler_cluster <- gp1
                                })
                   gts <- list_data_frame %>% 
                     distinct(gene_list) %>% 
                     separate(gene_list,c("Cluster","number_of_genes"),sep = "\nn = ",extra = "drop")
                   dt <- datatable(
                     gts,
                     colnames = names(gts),
                     rownames = FALSE,
                     filter = "none",
                     class = 'cell-border stripe compact',
                     options = list(
                       scrollX = TRUE,
                       scrollY = TRUE,
                       autoWidth = TRUE,
                       width = 20,
                       sDom  = '<"top">lrt<"bottom">ip',
                       info = FALSE,
                       paging = FALSE,
                       lengthChange = FALSE,
                       columnDefs = list(
                         list(className = 'dt-center ', targets = "_all"),
                         list(
                           targets = 0,
                           render = JS(
                             "function(data, type, row, meta) {",
                             "return type === 'display' && data.length > 25 ?",
                             "'<span title=\"' + data + '\">' + data.substr(0, 27) + '...</span>' : data;",
                             "}"
                           )
                         )
                       )
                     )
                   )
                   output$clustertable <- DT::renderDataTable(dt)
                 } else {
                   return()
                 }
               })
  
  # Cluster reset controller ----
  observeEvent(c(input$clusterGeneList,
                 input$clusterSamples), ignoreNULL = TRUE, ignoreInit = TRUE, {
                   reactive_values$clustergroups <- NULL
                   shinyjs::hide("hideclusterplots1")
                   shinyjs::hide("hideclustertable")
                 })
  
  # groupies tool action ----
  observeEvent(input$actiongroupiestool, ignoreInit = TRUE, {
    # print("groupies tool action")
    shinyjs::hide('plot1groupies')
    shinyjs::hide('plot2groupies')
    if (n_distinct(LIST_DATA$gene_file[[input$groupiesGeneList]]$full, na.rm = T) < as.numeric(input$selectgroupiesnumber) |
        is.null(input$groupiesSamples)) {
      showModal(modalDialog(
        title = "Information message",
        paste("Can't make more groups than number of genes"),
        size = "s",
        easyClose = TRUE
      ))
      return()
    }
    if(n_distinct(input$sliderbingroupies) < 2){
      start_end_bin <- c(reactive_values$slider_breaks$mylabels[1], last(reactive_values$slider_breaks$mylabels))
    } else {
      start_end_bin <- floor(reactive_values$slider_breaks$mybrakes[
        reactive_values$slider_breaks$mylabels %in% input$sliderbingroupies])
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <-
                     FindGroups(
                       LIST_DATA,
                       input$groupiesGeneList,
                       input$groupiesSamples,
                       start_end_bin
                     )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      if(LIST_DATA$STATE[1] != 0){
        LIST_DATA$STATE[1] <<- .5
      }
      reactive_values$groupiesgroups <- str_detect(names(LIST_DATA$gene_file),"^Groups_")
    } else {
      return()
    }
  })
  
  # Plots groupies based on number selected ----
  observeEvent(c(input$selectgroupiesnumber, reactive_values$groupiesgroups),
               ignoreInit = TRUE, ignoreNULL = TRUE,
               {
                 # print("groupies tool number")
                 if (is.null(reactive_values$groupiesgroups)) {
                   return()
                 }
                 shinyjs::hide('hidegroupiesplots1')
                 shinyjs::hide("hidegroupiestable")
                 shinyjs::hide('hidegroupiesplots2')
                 withProgress(message = 'Calculation in progress',
                              detail = 'This may take a while...',
                              value = 0,
                              {
                                LD <-
                                  GroupsNumList(
                                    LIST_DATA,
                                    input$groupiesGeneList,
                                    input$groupiesSamples,
                                    input$sliderbingroupies,
                                    input$selectgroupiesnumber
                                  )
                              })
                 if (!is_empty(LD$table_file)) {
                   LIST_DATA <<- LD
                   ol <- input$selectgroupiesfile
                   shinyjs::show('hidegroupiesplots1')
                   shinyjs::show("hidegroupiestable")
                   shinyjs::show('plot1groupies')
                   updateSelectInput(
                     session,
                     "selectgroupiesfile",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   LD$meta_data <- LD$meta_data %>%
                     dplyr::mutate(onoff=if_else(str_detect(gene_list,"^Groups_") &
                                                   set == input$groupiesSamples, set, "0"))
                   withProgress(message = 'Calculation in progress',
                                detail = 'This may take a while...',
                                value = 0,
                                {
                                  list_data_frame <- Active_list_data(LD, group="none", input$checkboxfull, 
                                                                      input$selectlegendnewline, input$selectlegendnewlinespace)
                                  if (!is_empty(list_data_frame)) {
                                    
                                    Apply_groupies_Math <- ApplyMath(
                                      list_data_frame,
                                      "mean",
                                      "none",
                                      0,
                                      0
                                    )
                                  }
                                  reactive_values$Plot_controler_groupies <- ggplot()
                                  reactive_values$Plot_controler_dgroupies <- ggplot()
                                  gp1 <-
                                    ggplot(Apply_groupies_Math ,aes(as.numeric(bin),value,color=gene_list)) +
                                    geom_line() +
                                    ylab("Mean bin value") +
                                    theme(legend.position="bottom",
                                          legend.title = element_blank(),
                                          axis.title.x=element_blank())
                                  print(gp1)
                                  reactive_values$Plot_controler_groupies <- gp1
                                })
                   gts <- list_data_frame %>%
                     distinct(gene_list) %>%
                     separate(gene_list,c("Groups","number_of_genes"),sep = "\nn = ",extra = "drop")
                   dt <- datatable(
                     gts,
                     colnames = names(gts),
                     rownames = FALSE,
                     filter = "none",
                     class = 'cell-border stripe compact',
                     options = list(
                       scrollX = TRUE,
                       scrollY = TRUE,
                       autoWidth = TRUE,
                       width = 20,
                       sDom  = '<"top">lrt<"bottom">ip',
                       info = FALSE,
                       paging = FALSE,
                       lengthChange = FALSE,
                       columnDefs = list(
                         list(className = 'dt-center ', targets = "_all"),
                         list(
                           targets = 0,
                           render = JS(
                             "function(data, type, row, meta) {",
                             "return type === 'display' && data.length > 25 ?",
                             "'<span title=\"' + data + '\">' + data.substr(0, 27) + '...</span>' : data;",
                             "}"
                           )
                         )
                       )
                     )
                   )
                   output$groupiestable <- DT::renderDataTable(dt)
                 } else {
                   return()
                 }
               })
  
  # groupies reset controller ----
  observeEvent(c(input$groupiesGeneList,
                 input$groupiesSamples), ignoreNULL = TRUE, ignoreInit = TRUE, {
                   reactive_values$groupiesgroups <- NULL
                   shinyjs::hide("hidegroupiesplots1")
                   shinyjs::hide("hidegroupiestable")
                   shinyjs::hide("hidegroupiesplots2")
                 })
  
  # Ratio tool action ----
  observeEvent(input$actionratiotool, ignoreInit = TRUE, {
    # print("ratio tool action")
    shinyjs::hide('ratio1table')
    shinyjs::hide('ratio2table')
    shinyjs::hide('ratio3table')
    numericratio <- input$numericratio
    if (is.numeric(input$numericratio)) {
      if (input$numericratio < 0) {
        updateNumericInput(session, "numericratio", value = 2)
        numericratio <- 2
      }
    } else {
      updateNumericInput(session, "numericratio", value = 2)
      numericratio <- 2
    }
    if (any(is.na(input$sliderbinratio2))) {
      updateSliderTextInput(session,
                            "sliderbinratio2",
                            selected = c("NA","NA"))
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <-
                     CompareRatios(
                       LIST_DATA,
                       input$selectratiofile,
                       input$pickerratio1file,
                       input$pickerratio2file,
                       floor(reactive_values$slider_breaks$mybrakes[
                         reactive_values$slider_breaks$mylabels %in% input$sliderbinratio1]),
                       input$sliderbinratio1,
                       floor(reactive_values$slider_breaks$mybrakes[
                         reactive_values$slider_breaks$mylabels %in% input$sliderbinratio2]),
                       input$sliderbinratio2,
                       numericratio,
                       input$checkratiozero,
                       floor(reactive_values$slider_breaks$mybrakes[
                         reactive_values$slider_breaks$mylabels %in% input$sliderRatioBinNorm]),
                       input$sliderRatioBinNorm
                     )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      if(LIST_DATA$STATE[1] != 0){
        LIST_DATA$STATE[2] <<- 0.5
      }
      
      ol <- input$selectratiofile
      updateSelectInput(
        session,
        "selectratiofile",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
      if (any(grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxratio1 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_file))]]$full$gene, na.rm = T),
            "Ratio Up file1",
            icon = icon("list"),
            color = "green"
          )
        })
      } else{
        output$valueboxratio1 <- renderValueBox({
          valueBox(0,
                   "Ratio Up file1",
                   icon = icon("list"),
                   color = "green")
        })
      }
      if (any(grep("Ratio_Down_file1\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxratio2 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Ratio_Down_file1\nn =",
                                                 names(LIST_DATA$gene_file))]]$full$gene, na.rm = T),
            "Ratio Down file1",
            icon = icon("list"),
            color = "blue"
          )
        })
      } else{
        output$valueboxratio2 <- renderValueBox({
          valueBox(0,
                   "Ratio Down file1",
                   icon = icon("list"),
                   color = "blue")
        })
      }
      if (any(grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxratio3 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_file))]]$full$gene, na.rm = T),
            "Ratio No Diff",
            icon = icon("list"),
            color = "yellow"
          )
        })
      } else{
        output$valueboxratio3 <- renderValueBox({
          valueBox(0,
                   "Ratio No Diff",
                   icon = icon("list"),
                   color = "yellow")
        })
      }
      if(!is.null(LIST_DATA$boxRatio)){
        updateNumericInput(session, "textboxmaxratio",
                           value = 0)
        updateNumericInput(session, "textboxminratio",
                           value = 0)
        my_range <- range(LIST_DATA$boxRatio$Ratio,na.rm = T) 
        
        updateNumericInput(session, "textboxmaxratio",
                           value = my_range[2])
        updateNumericInput(session, "textboxminratio",
                           value = my_range[1])
      } 
    } else {
      output$valueboxratio1 <- renderValueBox({
        valueBox(0,
                 "Ratio Up file1",
                 icon = icon("list"),
                 color = "green")
      })
      output$valueboxratio2 <- renderValueBox({
        valueBox(0,
                 "Ratio Down file1",
                 icon = icon("list"),
                 color = "blue")
      })
      output$valueboxratio3 <- renderValueBox({
        valueBox(0,
                 "Ratio No Diff",
                 icon = icon("list"),
                 color = "yellow")
      })
      return()
    }
  })
  
  # update Ratio violinplot ----
  observeEvent(c(input$textboxmaxratio, input$textboxminratio, input$checkboxviolinlog),ignoreInit = TRUE, ignoreNULL = TRUE,{
    if(!is.null(LIST_DATA$boxRatio)){
      if(input$checkboxviolinlog){
        shinyjs::disable("textboxminratio")
        shinyjs::disable("textboxmaxratio")
        gb <- ggviolin(LIST_DATA$boxRatio, x= "set", y = "Ratio", fill="set",
                       color="set",add = "boxplot", yscale = "log2",
                       add.params = list(fill = "white"),
                       xlab = "",ylab = "log2(Ratio)") + 
          theme(legend.position = 'none') 
      } else {
        shinyjs::enable("textboxminratio")
        shinyjs::enable("textboxmaxratio")
        my_range <- c(floor(input$textboxminratio), ceiling(input$textboxmaxratio)) 
        gb <- ggviolin(LIST_DATA$boxRatio, x= "set", y = "Ratio", fill="set",
                       color="set",add = "boxplot", add.params = list(fill = "white"),xlab = "") + 
          theme(legend.position = 'none') +
          coord_cartesian(ylim = my_range)
      
      }
      print(gb)
      reactive_values$Plot_controler_ratio <- gb 
    } 
    
  })
  
  # CDF tool action ----
  observeEvent(input$actioncdftool, ignoreInit = TRUE, {
    # print("CDF tool action")
    shinyjs::hide('plotcdf')
    shinyjs::hide('plotcdfscatter')
    if (any(between(
      floor(reactive_values$slider_breaks$mybrakes[
        reactive_values$slider_breaks$mylabels %in% input$sliderbincdf1]),
      floor(reactive_values$slider_breaks$mybrakes[
        reactive_values$slider_breaks$mylabels %in% input$sliderbincdf2[1]]),
      floor(reactive_values$slider_breaks$mybrakes[
        reactive_values$slider_breaks$mylabels %in% input$sliderbincdf2[2]])
    )) |
    any(between(
      floor(reactive_values$slider_breaks$mybrakes[
        reactive_values$slider_breaks$mylabels %in% input$sliderbincdf2]),
      floor(reactive_values$slider_breaks$mybrakes[
        reactive_values$slider_breaks$mylabels %in% input$sliderbincdf1[1]]),
      floor(reactive_values$slider_breaks$mybrakes[
        reactive_values$slider_breaks$mylabels %in% input$sliderbincdf1[2]])
    ))) {
      showModal(modalDialog(
        title = "Information message",
        paste("Bins regions should not overlab, \nBins set to 1/3 2/3"),
        size = "s",
        easyClose = TRUE
      ))
      updateSliderTextInput(
        session,
        "sliderbincdf1",
        selected = c(
          LIST_DATA$x_plot_range[1],
          floor(LIST_DATA$x_plot_range[2] / 4)
        )%>% replace(.,.==0,1))
      updateSliderTextInput(
        session,
        "sliderbincdf2",
        selected = c(
          ceiling(LIST_DATA$x_plot_range[2] / 4) + 1,
          LIST_DATA$x_plot_range[2]
        ))
    }
    ttt <-
      reactiveValuesToList(input)[paste0("-cdfspace3-", 
                                         gsub(" ", "-cdfspace2-", 
                                              gsub("\n", "-cdfspace1-", names(LIST_DATA$gene_file))))]
    checkboxonoff <- list()
    for (i in names(ttt)) {
      for (tt in ttt[i]) {
        selectgenelistonoff <-
          gsub("-cdfspace3-", "",
               gsub("-cdfspace2-", " ", 
                    gsub("-cdfspace1-", "\n", i)))
        checkboxonoff[[selectgenelistonoff]] <-
          c(checkboxonoff[[selectgenelistonoff]], tt)
      }
    }
    if (is_empty(checkboxonoff)) {
      return()
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <-
                     CumulativeDistribution(
                       LIST_DATA,
                       checkboxonoff,
                       floor(reactive_values$slider_breaks$mybrakes[
                         reactive_values$slider_breaks$mylabels %in% input$sliderbincdf1]),
                       input$sliderbincdf1,
                       floor(reactive_values$slider_breaks$mybrakes[
                         reactive_values$slider_breaks$mylabels %in% input$sliderbincdf2]),
                       input$sliderbincdf2, 
                       input$selectlegendnewline, input$selectlegendnewlinespace
                     )
                 })
    LD <<- LD
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      if(LIST_DATA$STATE[1] != 0){
        LIST_DATA$STATE[1] <<- 0.5
      }
      shinyjs::show('actioncdfdatatable')
      shinyjs::show('plotcdf')
      shinyjs::show('plotcdfscatter')
      newname <-
        grep("CDF ", names(LIST_DATA$gene_file), value = TRUE)
      rr <- range(LIST_DATA$gene_file[[newname]]$full$value,na.rm = T)
      # force trigger
      updateNumericInput(session,
                         "numericcdfmin",
                         value = 0)
      updateNumericInput(session,
                         "numericcdfmax",
                         value = 0)
      updateNumericInput(session,
                         "numericcdfmin",
                         value = floor(rr[1]))
      updateNumericInput(session,
                         "numericcdfmax",
                         value = ceiling(rr[2]))
    } else {
      shinyjs::disable('actioncdfcolor')
      return()
    }
  })
  
  # CDF x plot range ----
  observeEvent(c(input$numericcdfmin, input$numericcdfmax, 
                 reactive_values$df_options), ignoreInit = TRUE, {
                   # print("cdf plot observe range")
                   newname <-
                     grep("CDF ", names(LIST_DATA$gene_file), value = TRUE)
                   if(is_empty(newname)){
                     return()
                   }
                   df_options <-
                     LIST_DATA$meta_data %>%
                     dplyr::filter(gene_list ==  newname) %>%
                     dplyr::mutate(set = paste(
                       count,
                       sub(" - ","\n", plot_legend),sep = "\n")
                     )
                   df <- LIST_DATA$gene_file[[newname]]$full %>%
                     full_join(.,df_options %>% select(set,plot_legend),by="plot_legend") %>%
                     dplyr::mutate(set=set.y) %>% select(-set.x,-set.y)
                   
                   use_header <- pull(distinct(df_options, myheader))
                   if (n_groups(group_by(df_options, set)) == 2 &
                       n_distinct(df$gene, na.rm = T) > 1) {
                     tt_name <- pull(distinct(df_options, set))
                     tt <-
                       suppressWarnings(ks.test(pull(dplyr::filter(
                         df, set == tt_name[1]
                       ), value),
                       pull(dplyr::filter(
                         df, set == tt_name[2]
                       ), value)))
                     if (tt[[2]] == 0) {
                       use_header <- paste(use_header, "  p-value < 2.2e-16 ")
                     } else {
                       use_header <-
                         paste(use_header, paste("  p-value = ", format(tt[[2]], scientific = TRUE)))
                     }
                   }
                   mycdfss <- ggscatter(df %>% arrange(value), y = "value",x="gene",color = "plot_legend",
                             alpha=.8, palette = df_options$mycol,
                             ylim = as.numeric(c(input$numericcdfmin, input$numericcdfmax))) + 
                     rremove("x.text") + rremove("legend")
                   print(mycdfss)
                   mycdf <- GGplotC(df, df_options, use_header,as.numeric(c(input$numericcdfmin, input$numericcdfmax)))
                   output$plotcdf <- renderPlot({
                     mycdf
                   })
                   output$plotcdfscatter <- renderPlot({
                     mycdfss
                   })
                   shinyjs::enable('actioncdfcolor')
                 })
  
  # show DataTable ----
  observeEvent(input$pickerDT, ignoreNULL = TRUE, ignoreInit = TRUE,{
    if (length(LIST_DATA$gene_file) > 0) {
      gg <- LIST_DATA$gene_file[[input$pickerDT]]$full 
      dtg <- datatable(
        gg,
        colnames = names(gg),
        rownames = FALSE,
        filter = "top",
        class = 'cell-border stripe compact',
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          scrollY = TRUE,
          autoWidth = TRUE,
          width = 5,
          sDom  = '<"top">lrt<"bottom">ip',
          info = FALSE,
          paging = TRUE,
          lengthChange = FALSE,
          columnDefs = list(
            list(className = 'dt-center ', targets = "_all"),
            list(
              targets = 0,
              render = JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 64 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 59) + '...</span>' : data;",
                "}"
              )
            ),
            list(targets = -0, width = 5)
          )
        )
      )
      output$showgenelist <- DT::renderDataTable(dtg)
    }
  })
  
  
  # capture filtered DataTable ----
  observeEvent(input$showgenelist_rows_all, ignoreNULL = TRUE, ignoreInit = TRUE,{
    listname <- str_split_fixed(input$pickerDT,"\nn = ",2)[1]
    if(!str_detect(listname,"^filtered_")){
      listname <- paste0("filtered_",listname)
    } 
   
    if (length(input$showgenelist_rows_all) > 0 & 
        length(input$showgenelist_rows_all) != length(LIST_DATA$gene_file[[input$pickerDT]]$full$gene)) {
      
      if(any(str_detect(names(LIST_DATA$gene_file),listname))){
        # remove old filtered gene list
        old_names <- names(LIST_DATA$gene_file)[str_detect(names(LIST_DATA$gene_file),listname)]
        for(i in old_names){
          LIST_DATA$gene_file[[i]] <<- NULL
          LIST_DATA$meta_data <<- dplyr::filter(LIST_DATA$meta_data,
                                               gene_list != i)
        }
      } 
      listname <- paste0(listname,"\nn = ",length(input$showgenelist_rows_all))
      genelist <- LIST_DATA$gene_file[["Complete"]]$full[input$showgenelist_rows_all,]
      # record for info
      LIST_DATA$gene_file[[listname]]$full <<- genelist
      LIST_DATA$gene_file[[listname]]$info <<- tibble(loaded_info =
                                                         paste(listname,
                                                               Sys.Date()),
                                                       save_name = gsub(" ", "_", paste(listname, Sys.Date(), sep = "_")),
                                                       col_info = "gene"
      )
      LIST_DATA$meta_data <<- 
        distinct(bind_rows(LIST_DATA$meta_data,
                           LIST_DATA$meta_data %>% 
                             dplyr::filter(gene_list == names(LIST_DATA$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = listname,
                                           sub =  listname, 
                                           count = paste0("n = ",length(input$showgenelist_rows_all)),
                                           onoff = "0",
                                           plot_legend = " ")))
    }
  })
  # QC action ----
  observeEvent(input$buttonFilterzero, ignoreInit = TRUE, {
    # print("QC % tool") 
    out_list <- LIST_DATA$table_file %>% 
      dplyr::filter(set == input$QCsample)
    if(input$QCpickerplot == "% 0's per bin"){
      sortmin <- out_list %>% 
        group_by(set,bin) %>% 
        summarise(value = sum(score == 0)/n_distinct(gene)*100, .groups = "drop")
      ylab = "% zero"
    } else if(input$QCpickerplot == "quadrille"){
      qu <- out_list %>%
        group_by(gene) %>%
        summarise(set = sum(score, na.rm = TRUE),.groups="drop") %>% 
        dplyr::mutate(., set = ntile(set, 5))
      sortmin <- out_list %>% dplyr::select(-set) %>% inner_join(.,qu,by="gene") %>% 
        group_by(set,bin) %>% 
        summarise(value = sum(score), .groups = "drop") %>% 
        dplyr::mutate(set=paste("quadrille",set))
      ylab = "value"
    } else {
      if(input$QCpickerplot == "low range percentile"){
        my_per <- c(0.01, 0.025, 0.05, 0.075, 0.1)
      } else {
        my_per <- c(0.1, 0.25, 0.5, 0.75, 0.9, 1)
      }
      p_funs <- map(my_per, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
        set_names(my_per)
      sortmin <- out_list %>%
        group_by(set,bin) %>% summarize_at(vars(score), p_funs) %>% ungroup() %>% 
        dplyr::mutate(across(where(is.double),as.character)) %>% 
        gather(.,key = set,value = "value",-bin,-set) %>% 
        dplyr::mutate(value=as.double(value))
      ylab = "value"
    }
    if(!is.null(sortmin)){
      gp1 <-
        ggplot(sortmin, aes(as.numeric(bin), value, color=set)) + 
        geom_line() + 
        ylab(ylab) +
        theme(legend.position="bottom", 
              legend.title = element_blank(),
              axis.title.x=element_blank())
    } else{
      gp1 <- ggplot()
    }
    shinyjs::show("hidespinersQC")
    output$plotQC <- renderPlot({gp1})
    print(gp1)
  })
  
}

# UI -----
ui <- dashboardPage(
  skin = "purple-light",
  options = list(sidebarExpandOnHover = F),
  header = dashboardHeader(userOutput("user")),
  sidebar = dashboardSidebar(
    id = "sidebar",
    minified = TRUE,
    collapsed = TRUE,
    tags$head(
      tags$style(HTML("
        .inactiveLink {
          pointer-events: none;
          color: gray !important; 
          cursor: not-allowed;
        }
    
        .shiny-notification {
          position: fixed;
          top: 50%;
          left: 50%;
          transform: translate(-50%, -50%);
          color: purple;
          font-size: 20px;
          font-style: italic;
        }
        /* Outline all tabs */
        .nav-pills > li > a {
          color: white !important;       /* keep text white */
          border: 1px solid #ffffff;     /* white outline */
          border-radius: 6px;            /* rounded corners */
          margin: 1px;                   /* spacing between tabs */
        }
    
        /* Active tab styling */
        .nav-pills > li.active > a {
          color: white !important;
          background-color: purple-light !important; /* blue fill for active tab */
          border: 1px solid #ffffff;            /* keep outline */
          font-weight: bold;
        }
    
      "))
    )
    ,
    # disables tabs on start
    sidebarMenu(
      id = "leftSideTabs",
      menuItem("Load Data", tabName = "loaddata", icon = icon("file-import")),
      menuItem("Plot", tabName = "mainplot", icon = icon("chart-area")),
      menuItem("QC/Options", tabName = "qcOptions", icon = icon("clipboard-check")),
      menuItem("Group data", tabName = "grouptab", icon = icon("sitemap")),
      menuItem("Norm data", tabName = "filenorm", icon = icon("scale-balanced")),
      menuItem("Compare Lists", tabName = "genelists", icon = icon("grip-lines")),
      menuItem("Filter Tool", tabName = "sorttool", icon = icon("filter")),
      menuItem("Ratio Tool", tabName = "ratiotool", icon = icon("divide")),
      menuItem("Cluster Tools", tabName = "clustertool", icon = icon("object-group")),
      menuItem("Groups Tools", tabName = "groupiestool", icon = icon("object-group")),
      menuItem("CDF Tools", tabName = "cdftool", icon = icon("ruler-combined")),
      menuItem("Data Table Veiw", tabName = "DataTableTool", icon = icon("table"))
    )
  ),
  body = dashboardBody(
    useShinyjs(),
    tabItems(
      # load data tab ----
      tabItem(tabName = "loaddata",
              tabsetPanel(
                tabPanel(
                  "LOAD",
                  box(
                    status = "navy",
                    solidHeader = TRUE,
                    title = "Load .matrix.gz/URL.txt/url.tsv file",
                    width = 6,
                    style = "height: 150px;",
                    align = "center",
                    fileInput(
                      "filetable",
                      width = "75%",
                      label = "",
                      accept = NULL,
                      multiple = FALSE
                    ),
                    helpText("load table file(s)"),
                    br(),
                    hidden(div(id = "hidespiners", shinycssloaders::withSpinner(DT::dataTableOutput('loadedfilestable'), type = 4))),
                    DT::dataTableOutput('loadedfilestotaltable')
                  )
                ),
                tabPanel("SAVE",
                         hidden(div(
                           id = "startoff2",
                           box(
                             status = "navy",
                             solidHeader = TRUE,
                             title = "save .table",
                             width = 6,
                             style = "height: 150px;",
                             align = "center",
                             selectInput("selectsave", "", choices = "select list", ),
                             downloadButton("downloadGeneList", "Save List"),
                             helpText("save windowed bedGraph file(s)")
                           )
                         )))
              ),
              hidden(div(
                id = "startoff",
                box(
                  title = "Load Gene list, .tsv/.txt/.bed",
                  width = 6,
                  style = "height: 150px;" ,
                  solidHeader = TRUE,
                  status = "navy",
                  align = "center",
                  fileInput(
                    "filegene1",
                    width = "75%",
                    label = "",
                    accept = c('text/plain', 'application/gzip', 'application/x-gzip', '.txt', '.bed.gz', '.bed'),
                    multiple = FALSE
                  ),
                  helpText("load gene list"),
                  checkboxInput(inputId = "checkboxgenematch",label = "non-exact gene match"),
                  br(),
                  DT::dataTableOutput('loadedgenetable')
                )
              ))),
      tabItem(
        # mainplot ----
        tabName = "mainplot",
        fluidRow(box(
          width = 12, 
          status = "navy",
          solidHeader = TRUE,
          title = "",
          tags$style(".list {color:#00FF00}"),
          dropdownMenu = boxDropdown(
            icon = icon("list",class = "list"),
            boxDropdownItem("Update sample color/name", id = "dropcolor", icon = icon("palette")),
            dropdownDivider(),
            boxDropdownItem("Lines and Labels", id = "droplinesandlabels", icon = icon("chart-bar")),
            dropdownDivider(),
            boxDropdownItem("t-Test", id = "dropttest", icon = icon("chart-line"))
          ),
          sidebar = boxSidebar(
            id = "sidebarmath",
            width = 45,
            tags$style(".calculator {color:#FF0000}"),
            icon = icon("calculator", class = "calculator"),
            
            actionBttn(
              inputId = "actionMathUpDatePlot",
              label = "Update Line Plot",
              style = "unite",
              color = "primary",
              size = "md",
              block = TRUE,
              icon = icon("chart-line")
            ),
            
            br(),
            # Accordion organization
            tabsetPanel(
              id = "sidebar_tabs",
              type = "pills",
              
              # Tab 1: Basic Settings
              tabPanel(
                title = tagList(
                  tags$span(icon("chart-line"), style = "color:white"),
                  tags$span(" Main", style = "color:white")
                ),
                value = "basic",
                hr(),
                
                fluidRow(
                  column(
                    5,
                    selectInput("myMath",
                                label = "Math Function:",
                                choices = c("mean", "sum", "median", "var"),
                                selected = "mean"
                    )
                  ),
                  column(
                    6,
                    selectInput(
                      "selectplotnrom",
                      label = "Y Normalization:",
                      choices = c("none", "relative frequency", "rel gene frequency"),
                      selected = "none",
                      selectize = FALSE
                    )
                  )
                ),
                
                fluidRow(
                  column(
                    5,
                    selectInput(
                      "selectplotBinNorm",
                      label = "Bin Normalization:",
                      choices = c("NA"),
                      selected = "NA"
                    )
                  ),
                  column(
                    6,
                    awesomeRadio("checkboxbin",
                                 label = "Norm bin method:",
                                 choices = c("divide", "subtract"),
                                 selected = "divide",
                                 inline = TRUE)
                  )
                ),
                
                hidden(div(
                  id = "hideplotgroup",
                  selectInput("mygroup",
                              label = "Plot Group:",
                              choices = c("none", "groups"),
                              selected = "none"
                  )
                ))
              ),
              
              # Tab 2: Transform
              tabPanel(
                title = tagList(
                  tags$span(icon("wave-square"), style = "color:white"),
                  tags$span(" Smoothing & Transform", style = "color:white")
                ),
                value = "transform",
                br(),
                h4("Data Transformation", style = "margin-top: 0;"),
                hr(),
                
                awesomeCheckbox("checkboxsmooth", 
                                label = "Enable smoothing",
                                value = FALSE),
                
                conditionalPanel(
                  condition = "input.checkboxsmooth",
                  numericInput("numericsmooth", 
                               label = "Smoothing span (0-1):", 
                               value = 0.2,
                               min = 0,
                               max = 1,
                               step = 0.05)
                ),
                
                hr(),
                
                h5("Data Transformations:"),
                
                fluidRow(
                  column(
                    6,
                    awesomeCheckbox("checkboxlog2", 
                                    label = "Log2 transform",
                                    value = FALSE),
                    awesomeCheckbox("checkboxabs", 
                                    label = "Absolute value",
                                    value = FALSE)
                  ),
                  column(
                    6,
                    awesomeCheckbox("checkboxauc", 
                                    label = "Show AUC",
                                    value = FALSE),
                    awesomeCheckbox("checkboxfull", 
                                    label = "Include zeros",
                                    value = FALSE)
                  )
                )
              ),
              
              # Tab 3: Regions
              tabPanel(
                title = tagList(
                  tags$span(icon("ruler"), style = "color:white"),
                  tags$span(" Axis & Range", style = "color:white")
                ),
                value = "regions",
                br(),
                h4("Y-axis Range:", style = "margin-top: -5px; margin-bottom: -10px;"),
                hr(),
                
                fluidRow(div(style = "margin-top: -10px; margin-bottom: -10px;",
                  column(
                    4,
                    numericInput("numericYRangeLow", 
                                 label = "Y min:", 
                                 value = 0)
                  ),
                  column(
                    4,
                    numericInput("numericYRangeHigh", 
                                 label = "Y max:", 
                                 value = 0)
                  )
                )
                ),
                div(style = "margin-top: -10px; margin-bottom: -10px;",
                hr(),
                ),
                h5("X-axis Range (Bins):"),
                sliderTextInput(
                  "sliderplotBinRange",
                  label = NULL,
                  grid = TRUE,
                  choices = c("100", "100"),
                  selected = c("100", "100")
                ),
                
                tags$small("Adjust the bin range to zoom in/out on specific regions")
              ),
              
              # Tab 4: Statistics
              tabPanel(
                title = tagList(
                  tags$span(icon("play-circle"), style = "color:white"),
                  tags$span(" Box/Violin Plot", style = "color:white")
                ),
                value = "stats",
                br(),
                
                div(style = "margin-top: -10px;",
                hr()
                ),
                actionBttn(
                  inputId = "actionViolinPlot",
                  label = "Create Box/Violin Plot",
                  style = "unite",
                  color = "warning",
                  size = "md",
                  block = TRUE,
                  icon = icon("chart-bar")
                ),
                
                fluidRow(
                  column(
                    4,
                    div(style = "margin-top: 20px;",
                    awesomeCheckbox("checkboxlog2Violin", 
                                    label = "Log2 transform",
                                    value = FALSE)
                  )),
                  column(
                    4,
                    selectInput("VplotType",
                                label = "plot type:",
                                choices = c("violin", "boxplot", "both"),
                                selected = "boxplot")
                    ),
                  column(
                    4,
                    numericInput("VmergeBin", 
                                 label = "Aggregate bins:", 
                                 value = 10)
                  )
                ),
                fluidRow(div(style = "margin-top: -10px; margin-bottom: -10px;",
                             column(
                               4,
                               numericInput("numericYRangeLowViolin", 
                                            label = "Y min:", 
                                            value = 0)
                             ),
                             column(
                               4,
                               numericInput("numericYRangeHighViolin", 
                                            label = "Y max:", 
                                            value = 0)
                             )
                )
                ),
                
                
                hr(),
                
                tags$small(
                  icon("info-circle"),
                  " Click 'Update Line Plot' to refresh the main visualization or 'Create Box/Violin Plot' to view distribution plots."
                )
              )
            )
          ),
          shinycssloaders::withSpinner(plotOutput("plot"), type = 4),
          div(
            id = "actionmyplotshow",
            style = "position: absolute; z-index: 1; left: 45%; top: 50%;",
            actionButton(
              "actionmyplot",
              "Update Plot",
              icon = icon("chart-area"),
              style = "color: #fff; background-color: #337ab7; border-color: #2e6da4;"
            )
          )
        )),
        fluidRow(box(
          width = 12, 
          status = "primary",
          solidHeader = TRUE,
          title = "", 
          box(title = "Main",
              width = 6,
              status = "navy",
              solidHeader = T,
              collapsible = T,
              collapsed = F,
              uiOutput("DynamicGenePicker_main")
          ),
          hidden(
            div(
              id = "showpickersort",
              box(
                title = "Filter (max 4 lists)",
                width = 6,
                status = "navy",
                solidHeader = T,
                collapsible = T,
                collapsed = F,
                uiOutput("DynamicGenePicker_sort")
              )
            )),
          hidden(
            div(
              id = "showpickercomparisons",
              box(
                title = "Gene comparisons",
                width = 6,
                status = "navy",
                solidHeader = T,
                collapsible = T,
                collapsed = F,
                uiOutput("DynamicGenePicker_comparisons")
              )
            )),
          hidden(
            div(
              id = "showpickerratio",
              box(
                title = "Ratio",
                width = 6,
                status = "navy",
                solidHeader = T,
                collapsible = T,
                collapsed = F,
                uiOutput("DynamicGenePicker_ratio")
              )
            )),
          hidden(
            div(
              id = "showpickercluster",
              box(
                title = "Clusters",
                width = 6,
                status = "navy",
                solidHeader = T,
                collapsible = T,
                collapsed = F,
                uiOutput("DynamicGenePicker_clusters")
              )
            )),
          hidden(
            div(
              id = "showpickergroupies",
              box(
                title = "Groups",
                width = 6,
                status = "navy",
                solidHeader = T,
                collapsible = T,
                collapsed = F,
                uiOutput("DynamicGenePicker_groupies")
              )
            )),
          hidden(
            div(
              id = "showpickercdf",
              box(
                title = "CDF",
                width = 6,
                status = "primary",
                solidHeader = T,
                collapsible = T,
                collapsed = F,
                uiOutput("DynamicGenePicker_cdf")
              )
            ))
        ))
      ),
      tabItem(
        # QC ----
        tabName = "qcOptions",
        box(
          width = 12,
          status = "primary",
          title = "QC Options",
          solidHeader = T,
          hidden(div(id = "hidespinersQC", shinycssloaders::withSpinner(plotOutput("plotQC"), type = 4)))
        ),
        box(width = 12,
            status = "primary",
            title = "QC Options",
            solidHeader = T,
            pickerInput("QCsample",
                        label = "select file",
                        width = "60%",
                        choices = "Load data file",
                        multiple = F,
                        options = list(title = "Select file")),
            pickerInput("QCpickerplot",
                        label = "select plot type",
                        width = "60%",
                        choices = c("low range percentile", "braod range percentile",
                                    "% 0's per bin","quadrille"),
                        selected = "quadrille",
                        multiple = F,
                        options = list(title = "Select file")),
            actionButton("buttonFilterzero",label = "Plot")
        )
      ),
      tabItem(
        # grouptab ----
        tabName = "grouptab",
        div(
          box(
            title = "Select files for setting groups",
            width = 12,
            status = "primary",
            solidHeader = T,
            collapsible = T,
            align = "center",
            div(style = "padding-left: 15%;",
                fluidRow(
                  pickerInput(
                    "pickergroupsample",
                    label = "Pick group samples",
                    width = "90%",
                    choices = "Load data file",
                    multiple = T,
                    options = list(title = "Select at least 2 files",
                                   `selected-text-format` = "count > 0"
                    )
                  ),
                  helpText(icon("exclamation-triangle"),"samples can be set in more that one group at a time"),
                  DT::dataTableOutput('loadedfilestable2')
                )),
            div(style = "padding-left: 15%;",
                fluidRow(
                  textInput("textgroupname", "group file name",
                            width = "90%",)),
                column(4, style = "padding-top: 4%;",
                       actionButton("actiongroup", label = "create group"),
                       helpText("will have the same color as the top sample in the group"))
            )
          ),
          box(
            title = "Select math for combinding group to one sample",
            width = 12,
            status = "primary",
            solidHeader = T,
            collapsible = T,
            div(style = "padding-left: 15%;",
                fluidRow(
                  pickerInput(
                    "pickergroupmath",
                    label = "Math",
                    width = "90%",
                    choices = c("mean", "sum", "median"),
                    selected = "mean"
                  )
                )),
            div(style = "padding-left: 15%;",
                fluidRow(
                  column(4, style = "padding-top: 4%;",
                         actionButton("actionnormgroup", label = "create file"))
                )),
            valueBoxOutput("valueboxgroupfile")
          )
        )
      ),
      tabItem(
        # filenorm ----
        tabName = "filenorm",
        div(
          box(
            title = "Select files for normalization",
            width = 12,
            status = "primary",
            solidHeader = T,
            collapsible = T,
            div(style = "padding-left: 15%;",
                fluidRow(
                  pickerInput(
                    "pickernumerator",
                    label = "numerator",
                    width = "90%",
                    choices = "Load data file",
                    multiple = F,
                    options = list(title = "Select first file")
                  )
                )),
            div(style = "padding-left: 15%;",
                fluidRow(
                  column(
                    3,
                    radioGroupButtons(
                      "adddata",
                      label = "",
                      status = "primary",
                      choices = c("/", "+", "-"),
                      selected = "/"
                    )
                  ),
                  column(4, style = "padding-top: 4%;",
                         actionButton("actionnorm", label = "create file"))
                )),
            div(style = "padding-left: 15%;",
                fluidRow(
                  pickerInput(
                    "pickerdenominator",
                    label = "denominator",
                    width = "90%",
                    choices = "Load data file",
                    multiple = F,
                    options = list(title = "Select second file")
                  )
                ),
                fluidRow(
                  textInput("textnromname", "Norm file name",
                            width = "90%",))
            ),
            awesomeRadio(
              "radiogenebygene",
              label = "",
              choices = c("bin by bin", "mean of bins by mean of bins"),
              selected = "bin by bin"
            ),
            awesomeCheckbox(
              "checkboxnormzero",
              label = "replace denom 0's with min/2",value = FALSE),
            awesomeCheckbox(
              "checkboxnormzero2",
              label = "replace all 0's with min/2",value = FALSE),
            valueBoxOutput("valueboxnormfile")
          )
        )
      ),
      tabItem(
        # genelists ----
        tabName = "genelists",
        div(
          id = "enablemaingenelists",
          box(
            title = "Gene Lists",
            status = "primary",
            solidHeader = T,
            width = 12,
            fluidRow(column(width = 4,
                            pickerInput(
                              inputId = "pickergenelists",
                              label = "Select Gene lists",
                              choices = "Load data file",
                              multiple = T,
                              options = list(`selected-text-format` = "count > 1")
                            )
            )
            ),
            actionButton("actiongenelists", "Compare Gene lists"),
            helpText("Shows innerjoined, Exlusive, and Total gene lists")
          ),
          box(
            title = "Gene List Tables",
            status = "primary",
            solidHeader = T,
            collapsible = T,
            width = 12,
            helpText("Needs at least 2 gene lists"),
            actionButton("actiongenelistsdatatable", "Show gene list"),
            tabBox(
              id = "geneliststooltab",
              width = 12,
              tabPanel(
                "innerjoined Gene Lists",
                helpText("All filtering applied to gene list usage elsewhere"),
                DT::dataTableOutput('genelists1table')
              ),
              tabPanel(
                "Total Gene Lists",
                helpText("All filtering applied to gene list usage elsewhere"),
                DT::dataTableOutput('genelists2table')
              ),
              tabPanel(
                "antijoin Gene Lists",
                helpText("All filtering applied to gene list usage elsewhere"),
                DT::dataTableOutput('genelists3table')
              )
            )
          ),
          fluidRow(
            valueBoxOutput("valueboxgene1"),
            valueBoxOutput("valueboxgene2"),
            valueBoxOutput("valueboxgene3")
          )
        )
      ),
      tabItem(
        # filter/sort tools ----
        tabName = "sorttool",
        div(
          id = "enablemainsort",
          box(
            width = 12,
            solidHeader = T,
            status = "primary",
            collapsible = T,
            collapsed = T,
            title = "Size and Separation filter",
            column(3,
                   numericInput(
                     "geneSeparation",
                     label = "gene separated by bp",
                     value = 0,
                     step = 100,
                     min = 0,
                     max = 1e7
                   ),
                   awesomeCheckbox("checkboxStranded","stranded?"),
                   helpText("0 = no filter")
            ),
            column(3,
                   numericInput(
                     "geneSizeMin",
                     label = "Min gene size bp",
                     value = 0,
                     step = 100,
                     min = 0,
                     max = 1e7
                   ),
                   helpText("0/empty = no filter")
            ),
            column(3,
                   numericInput(
                     "geneSizeMax",
                     label = "Max gene size bp",
                     value = 0,
                     step = 100,
                     min = 0,
                     max = 1e7
                   ),
                   helpText("0/empty = no filter")
            ),
            column(3,
                   helpText("filters on separation then size"),
                   actionButton("actionSizeSep", "filter")
            )
          ),
          box(
            width = 12, solidHeader = T,
            status = "primary",
            column(width = 6,
                   pickerInput("sortGeneList", label = "select list",
                               choices = (LIST_DATA$meta_data)),
                   div(
                     style = "margin-bottom: -20px;",
                     sliderTextInput(
                       "slidersortbinrange",
                       label = "Select Bin Range:",
                       grid = TRUE,
                       c("100","100"),
                       selected = c("100","100")
                     )
                   )
            ),
            column(width = 6,
                   pickerInput("sortSamples", label = "select sample(s)",
                               choices = "select sample(s)",selected = "select sample(s)",
                               multiple = TRUE,
                               options = list(
                                 `actions-box` = FALSE,
                                 `selected-text-format` = "count > 0"))
            )
          ),
          box(
            title = "Filter by sum rank",
            solidHeader = T,
            width = 6,
            status = "navy",
            collapsible = T,
            fluidRow(column(
              12,
              style = "margin-bottom: -20px;",
              sliderInput(
                "slidersortpercent",
                label = "% select:",
                post = "%",
                min = 1,
                max = 100,
                value = 75,
                step = 0.2
              )
            ),
            ),
            fluidRow(align="center",
                     column(
                       6,
                       pickerInput(
                         "selectsorttop",
                         "Filter to keep",
                         choices = c("Top%", "Middle%", "Bottom%"),
                         selected = "Middle%"
                       )
                     )),
            fluidRow(align="center",
                     actionButton("actionsorttool", "filter sum")
            )
          ),
          box(
            title = "Filter on average genes",
            solidHeader = T,
            width = 6,
            status = "navy",
            collapsible = T,
            fluidRow(align="center",
                     helpText(HTML("For selected bins creates gene lists<br>
                       &nbsp;&nbsp;&nbsp;&nbsp;1. All bins above<br>
                       &nbsp;&nbsp;&nbsp;&nbsp;2. All bins below<br>
                       &nbsp;&nbsp;&nbsp;3. Mixed"))),
            fluidRow(align="center",
                     column(
                       6,
                       style = "margin-bottom: 10px;",
                       pickerInput(
                         "selectaveragemath",
                         choices = c("mean", "median"),
                         selected = "mean"
                       )
                     ),
                     column(
                       6,
                       helpText("Tip: Use 1 sample and small number of bins")
                       )),
            fluidRow(align="center",
                     actionButton("actionaveragetool", "filter list")
            )
          ),
          box(
            title = "filter peaks",
            solidHeader = T,
            width = 6,
            status = "navy",
            collapsible = T,
            fluidRow(column(
              12,
              style = "margin-bottom: 10px;",
              numericInputIcon("peakfilternum",
                               "signal hight", 
                               value = "1",
                               step = ".25"
              ),
              uiOutput('rangeHelptext'),
            )),
            fluidRow(align="center",
                     column(6,
                            pickerInput(
                              "selectsortpeak",
                              "Filter out Option",
                              choices = c("peak","keep peak"),
                              selected = "peak"
                            )
                     )),
            fluidRow(align="center",
                     actionButton("actionsortpeak", "filter")
            )
          ),
          box(
            title = "filter by percentile distribution",
            solidHeader = T,
            width = 6,
            status = "navy",
            collapsible = T,
            fluidRow(column(
              6,
              numericInputIcon("numericsortmin",
                               "min", 
                               value = "1",
                               max = "100", min="1",
                               step = ".25",
                               icon = icon("percent")
              )
            ),
            column(
              6,
              numericInputIcon("numericsortmax",
                               "max", 
                               value = "99.5",
                               max = "100", min="1",
                               step = ".25",
                               icon = icon("percent")
              )
            )
            ),
            fluidRow(align="center",column(6,
                                           pickerInput(
                                             "selectsortper",
                                             "Filter Option",
                                             choices = c("min%", "between%", "max%"),
                                             selected = "min%"
                                           )
            )),
            fluidRow(align="center",
                     actionButton("actionsortper", "filter percentile")
            ),
            helpText("Hint: use 1 file to display a range of %'s")
          ),
          div(
            id = "hidesortplots1",
            box(headerBorder = F,
                width = 6,
                withSpinner(plotOutput("plot1sort",height = "200px"), type = 4)
            )
          ),
          div(
            id = "hidesortplots2",
            box(headerBorder = F,
                width = 6,
                withSpinner(plotOutput("plot2sort",height = "200px"), type = 4)
            )
          ),
          valueBoxOutput("valueboxsort")
        )
      ),
      tabItem(
        # ratio tool ----
        tabName = "ratiotool",
        div(
          id = "enablemainratio",
          box(title = "Ratio tool",
              status = "primary",
              solidHeader = T,
              width = 12,
              column(width = 6,
                     selectInput(
                       inputId = "selectratiofile",
                       label = "Select gene list to sort on",
                       choices = "Load data file",
                       width = "99%"
                     ),
                     actionButton("actionratiotool", "Get fold changes"),
                     awesomeCheckbox(
                       "checkratiozero",
                       label = "replace denom 0's with min/2",
                       value = FALSE
                     )
              ),
              column(width = 6,
                     pickerInput(
                       inputId = "pickerratio1file",
                       width = "99%",
                       label = "Select first file",
                       choices = "Load data file",
                       multiple = F,
                       options = list(title = "Select first file")
                     ),
                     pickerInput(
                       inputId = "pickerratio2file",
                       width = "99%",
                       label = "Select second file",
                       choices = "Load data file",
                       multiple = F,
                       options = list(title = "Select second file")
                     )
              )
          ),
          box(
            title = "Ratio tool",
            status = "primary",
            solidHeader = T,
            width = 12,
            fluidRow(
              column(
                2,
                numericInput(
                  "numericratio",
                  "Fold Change",
                  value = 2,
                  min = 0,
                  max = 10,
                  step = 0.1
                )
              ),
              column(
                5,
                sliderTextInput(
                  "sliderbinratio1",
                  label = "Select Bin Range:",
                  grid = TRUE,
                  choices = c("100","100"),
                  selected = c("100","100")
                )
              ),
              column(
                5,
                sliderTextInput(
                  "sliderbinratio2",
                  label = "Select Bin Range:",
                  grid = TRUE,
                  choices = c("100","100"),
                  selected = c("100","100")
                )
              )
            ),
            helpText("(file1[1]/file1[2])/(file2[1]/file2[2]) or file1[1]/file2[2]"),
            sliderTextInput(
              "sliderRatioBinNorm",
              label = "Select Bin To Norm first:",
              grid = TRUE,
              choices = c("NA"),
              selected = c("NA")
            )
          ),
          box(
            title = "Violin  Plot",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            withSpinner(plotOutput("plotratio"), type = 4),
            numericInput(inputId = 'textboxmaxratio',
                         "yaxis max",
                         value = 0,
                         min = 0,
                         max = 1000,
                         step = .5),
            numericInput(inputId = 'textboxminratio',
                         "yaxis min",
                         value = 0,
                         min = 0,
                         max = 1000,
                         step = .5),
            checkboxInput(inputId = 'checkboxviolinlog',
                          label = "log2",value = TRUE)
          ),
          fluidRow(
            valueBoxOutput("valueboxratio1"),
            valueBoxOutput("valueboxratio2"),
            valueBoxOutput("valueboxratio3")
          )
        )
      ),
      tabItem(
        # cluster tools ----
        tabName = "clustertool",
        div(
          id = "enablemaincluster",
          box(title = "Cluster tools",
              status = "primary",
              solidHeader = T,
              width = 6,
              pickerInput("clusterGeneList", label = "select list",
                          choices = (LIST_DATA$meta_data)),
              pickerInput("clusterSamples", label = "select sample",
                          choices = "select sample",selected = "select sample",
                          multiple = F
              )
              
          ),
          box(
            title = "Cluster tools",
            status = "primary",
            solidHeader = T,
            width = 6,
            id = "test",
            style = "margin-bottom: 15px;",
            fluidRow(column(
              4,
              selectInput(
                inputId = "selectclusternumber",
                label = "Select number of clusters",
                choices = c(10:2),
                selected = 4,
                width = "99%"
              )
            ),
            column(
              8,
              sliderTextInput(
                "sliderbincluster",
                label = "Select Bin Range:",
                grid = TRUE,
                choices = c("100","100"),
                selected = c("100","100")
              )
            )),
            column(width = 5,
                   pickerInput("clustpattern",label = "Cluster on",
                               choices = c("pattern","expression"),
                               selected = "expression",multiple = F)),
            column(width = 6,
                   checkboxInput("clusterRF","plot relative frequency",value = T)),
            column(width = 6,
                   actionButton("actionclustertool", "Get clusters"))
          ),
          div(
            id = "hideclusterplots1",
            box(headerBorder = F,
                width = 8,
                withSpinner(plotOutput("plot1cluster",height = "300px"), type = 4)
            )
          ),
          div(
            id = "hideclustertable",
            box(headerBorder = F,
                style = "padding: 0px 2px;",
                width = 4,
                DT::dataTableOutput('clustertable',height = "320px")     
            ))
        )
      ),
      tabItem(
        # groupies tools ----
        tabName = "groupiestool",
        div(
          id = "enablemaingroupies",
          box(title = "Groups tools",
              status = "primary",
              solidHeader = T,
              width = 6,
              pickerInput("groupiesGeneList", label = "select list",
                          choices = (LIST_DATA$meta_data)),
              pickerInput("groupiesSamples", label = "select sample",
                          choices = "select sample",selected = "select sample",
                          multiple = F
              )
              
          ),
          box(
            title = "groupies tools",
            status = "primary",
            solidHeader = T,
            width = 6,
            id = "groupies_test",
            style = "margin-bottom: 15px;",
            fluidRow(column(
              4,
              selectInput(
                inputId = "selectgroupiesnumber",
                label = "Select number of groupies",
                choices = c(10:2),
                selected = 4,
                width = "99%"
              )
            ),
            column(
              8,
              sliderTextInput(
                "sliderbingroupies",
                label = "Select Bin Range:",
                grid = TRUE,
                choices = c("100","100"),
                selected = c("100","100")
              )
            )),
            column(width = 6,
                   actionButton("actiongroupiestool", "Get groups"))
          ),
          div(
            id = "hidegroupiesplots1",
            box(headerBorder = F,
                width = 8,
                withSpinner(plotOutput("plot1groupies",height = "300px"), type = 4)
            )
          ),
          div(
            id = "hidegroupiestable",
            box(headerBorder = F,
                style = "padding: 0px 2px;",
                width = 4,
                DT::dataTableOutput('groupiestable',height = "320px")
            )),
          div(
            id = "hidegroupiesplots2",
            box(headerBorder = F,
                width = 12,
                withSpinner(plotOutput("plot2groupies",height = "200px"), type = 4)
            )
          )
        )
      ),
      tabItem(
        # cdf ----
        tabName = "cdftool",
        box(
          title = "CDF tool",
          status = "primary",
          solidHeader = T,
          collapsible = T,
          width = 12,
          box(title = "Main",
              width = 6,
              status = "navy",
              solidHeader = T,
              collapsible = T,
              collapsed = F,
              uiOutput("DynamicCDFPicker_main")
          ),
          hidden(
            div(
              id = "showpickersort_cdf",
              box(
                title = "Filter (max 4 lists)",
                width = 6,
                status = "navy",
                solidHeader = T,
                collapsible = T,
                collapsed = T,
                uiOutput("DynamicCDFPicker_sort")
              )
            )),
          hidden(
            div(
              id = "showpickercomparisons_cdf",
              box(
                title = "Gene comparisons",
                width = 6,
                status = "navy",
                solidHeader = T,
                collapsible = T,
                collapsed = T,
                uiOutput("DynamicCDFPicker_comparisons")
              )
            )),
          hidden(
            div(
              id = "showpickerratio_cdf",
              box(
                title = "Ratio",
                width = 6,
                status = "navy",
                solidHeader = T,
                collapsible = T,
                collapsed = T,
                uiOutput("DynamicCDFPicker_ratio")
              )
            )),
          hidden(
            div(
              id = "showpickercluster_cdf",
              box(
                title = "Clusters",
                width = 6,
                status = "navy",
                solidHeader = T,
                collapsible = T,
                collapsed = T,
                uiOutput("DynamicCDFPicker_clusters")
              )
            ))
          ,
          hidden(
            div(
              id = "showpickergroupies_cdf",
              box(
                title = "Groups",
                width = 6,
                status = "navy",
                solidHeader = T,
                collapsible = T,
                collapsed = T,
                uiOutput("DynamicCDFPicker_groupies")
              )
            ))
        ),
        box(
          title = "CDF Plot",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          column(
            width = 6,
            sliderTextInput(
              "sliderbincdf1",
              label = "Select numerator Bin Range:",
              grid = TRUE,
              choices = c("0","100"),
              selected = c("0","100")
            )),
          column(
            width = 6,
            sliderTextInput(
              "sliderbincdf2",
              label = "Select denominator Bin Range:",
              grid = TRUE,
              choices = c("100","100"),
              selected = c("100","100")
            )),
          fluidRow(column(
            width = 2,
            numericInput(
              "numericcdfmin",
              label = "Xaxis min",
              value = 0
            )),
            column(
              width = 2,
              numericInput(
                "numericcdfmax",
                label = "Xaxis max",
                value = 1
              ))),
          actionButton("actioncdftool", "Plot CDF"),
          actionButton("actioncdfcolor", "Set Plot colors")
        ),
        box(
          title = "CDF Plot",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          withSpinner(plotOutput("plotcdf"), type = 4)
        ),
        box(
          title = "Scatter Plot",
          status = "primary",
          solidHeader = TRUE,
          width = 12,
          collapsible = TRUE,
          withSpinner(plotOutput("plotcdfscatter"), type = 4)
        )
      ),
      tabItem(
        # DataTable ----
        tabName = "DataTableTool",
        box(
          status = "primary",
          solidHeader = T,
          collapsible = T,
          width = 12,
          title = "Data Table",
          pickerInput("pickerDT",
                      label = "select gene list",
                      choices = "select list",selected = "select list",
                      multiple = F),
          DT::dataTableOutput('showgenelist')
        )
      )
    )
  ),
  controlbar = dashboardControlbar(disable = TRUE),
  title = "DashboardPage"
)

# execute ----
shinyApp(ui = ui, server = server)