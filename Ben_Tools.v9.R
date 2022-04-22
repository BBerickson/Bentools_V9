# Created by Benjamin Erickson BBErickson@gmail.com

source("R_scripts/helpers.R", local = TRUE)

# run load needed packages using my_packages(x) ----
  # "factoextra" for dendogram plot fviz_dend()
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
    "DT",
    "patchwork",
    "zip",
    "ggpubr",
    "fastcluster",
    "factoextra"
  )
))

source("R_scripts/functions.R", local = TRUE)

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 500MB. ----
options(shiny.maxRequestSize = 500 * 1024 ^ 2)

LIST_DATA <<- list(
  table_file = NULL,
  # gene bin score set
  gene_file = NULL,
  # holds $Complete genes from files and $gene file(s)
  gene_info = NULL,
  # for holding meta data gene file(s) [c("gene_list", "count", "set", "color", plot?, "legend", "plot_set")]
  ttest = NULL,
  # t.test results $full is for numbers $gene_info for holding plotting options
  clust = NULL,
  # Cluster holder
  x_plot_range = c(0, 0),
  STATE = c(0, 0, 5) # flow control,
  # [1] 1 = at least one file has been loaded and lets reactive fill in info
  #
  # [2] 0 = first time switching tab auto plotting
  #     1 = hidden plot button, reactive for plot enabled
  #     2 = on/off reactive picker changed, shows plot button, reactive for plot disabled
  # [3] line and label type i.e. 5, 4, 3, 543
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
    Y_Axis_numbers = c(0,100),
    Lines_Labels_List = list(mybrakes="",mylabels="",
                             myline = tibble(use_virtical_line_color=c("green","red","black","black"),
                                             use_virtical_line_type=c("dotted","dotted","solid","solid")),
                             mysize = c(2.0, 2.5, 13.0, 13.0, 10.0),
                             myset = c(20, 40, 15, 45, 100, 5)),
    Picker_controler = NULL,
    mymath = c("mean", "none", "0", "FALSE", "FALSE", "1", "80", "none"),
    ttest = NULL,
    ttest_values = c("none", "wilcox.test", "two.sided", "FALSE", "FALSE", "-log10", "fdr"),
    ttest_options = c(0, 0, 1, "select sample", 0.05),
    clustergroups = NULL,
    cluster_control = NULL,
    setsliders = NULL
  )
  
  output$user <- renderUser({
    dashboardUser(
      name = "BenTools V9.a",
      image = "ben head.jpg",
      title = "Benjamin Erickson",
      subtitle = "BBErickson@gmail.com"
    )
  })
  # disables tabs on start ----
  addCssClass(selector = "a[data-value='mainplot']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='qcOptions']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='filenorm']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='genelists']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='sorttool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='ratiotool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='clustertool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='cdftool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='DataTableTool']", class = "inactiveLink")
  
  # loads data file(s) ----
  observeEvent(input$filetable, {
    print("load file")
    shinyjs::disable("startoff")
    shinyjs::disable("startoff2")
    shinyjs::show("hidespiners")
    output$loadedfilestotaltable <- DT::renderDataTable(datatable())
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                
    meta_data <- PrepMetaFile(input$filetable$datapath,
                              input$filetable$name)
    if (!is_empty(meta_data)) {
    for (i in seq_along(meta_data$filepath)) {
      setProgress(i/length(meta_data$filepath), 
                  detail = paste("Gathering info on",meta_data$nick[i]))
      if (meta_data$nick[i] %in% names(LIST_DATA$table_file)) {
        showModal(modalDialog(
          title = "Information message",
          paste(meta_data$nick[i], "has already been loaded"),
          size = "s",
          easyClose = TRUE
        ))
        next()
      }
      bin_colname <- tableTestbin(meta_data[i, ])
      if(is_empty(bin_colname)){
        meta_data <- meta_data %>% filter(filepath != meta_data$filepath[i])
        next()
      }
      meta_data$color[i] <- RgbToHex(meta_data$color[i])
      shinyjs::reset("filetable")
      setProgress(i/length(meta_data$filepath), 
                  detail = paste("Loading file",meta_data$nick[i]))
      LD <- LoadTableFile(meta_data[i, ], bin_colname)
      setProgress(i/length(meta_data$filepath), 
                  detail = paste("Testing compatiblity",meta_data$nick[i]))
        if (LIST_DATA$x_plot_range[2] == 0) {
          LIST_DATA$x_plot_range <<- range(LD$bin)
          LIST_DATA$STATE[3] <<- meta_data$type[i]
        } else if (max(LD$bin) != LIST_DATA$x_plot_range[2]) {
          showModal(
            modalDialog(
              title = "Information message",
              "Can't load file, different number of bins",
              size = "s",
              easyClose = TRUE
            )
          )
          next()
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
          next()
        }
    setProgress(i/length(meta_data$filepath), 
                detail = paste("Finishing up",meta_data$nick[i]))
          # make complete gene list
          LIST_DATA$gene_file$Complete$full <<-
            full_join(LD, LIST_DATA$gene_file$Complete$full, by = "gene") %>%
            distinct(gene)
          LIST_DATA$gene_info <<- LIST_DATA$gene_info %>%
            dplyr::mutate(count = if_else(gene_list == "Complete",
                                          paste("n =", n_distinct(LIST_DATA$gene_file$Complete$full, na.rm = T)), count))
        } else {
          LIST_DATA$gene_file$Complete$full <<- distinct(LD, gene)
        }
        if (LIST_DATA$STATE[2] == 0 &
            n_distinct(LIST_DATA$table_file$set) < 2) {
          oo <- meta_data$nick[i]
        } else {
          oo <- "0"
        }
        LIST_DATA$table_file <<- distinct(bind_rows(LIST_DATA$table_file, LD))
        LIST_DATA$gene_file$Complete$info <<- tibble(loaded_info = paste("all loaded genes",
                                                                             Sys.Date()))
        LIST_DATA$gene_info <<- distinct(bind_rows(LIST_DATA$gene_info,tibble(
          gene_list = "Complete",
          count = paste("n =", n_distinct(LIST_DATA$gene_file$Complete$full, na.rm = T)),
          set = meta_data$nick[i],
          mycol = meta_data$color[i],
          onoff = oo,
          sub = " ",
          plot_set = " ",
          group = meta_data$nick[i]
        )))
        ### check if grep of gene has occurred
        for(gg in seq_along(LIST_DATA$gene_file)[-1]){
          if("org_gene" %in% LIST_DATA$gene_file[[gg]]$full){
            setProgress(i/length(meta_data$filepath), 
                        detail = paste("re-matching gene list",meta_data$nick[i]))
            new_gene_match <- MatchGenes(LIST_DATA$gene_file[[1]]$full,
                                         LIST_DATA$gene_file[[gg]]$full %>% select(org_gene) %>%
                                           dplyr::rename(gene = org_gene))
            if (n_distinct(new_gene_match$gene, na.rm = T) != 0) {
              # fix name, fix info
              listname <- names(LIST_DATA$gene_file)[gg]
              LIST_DATA$gene_info <<- LIST_DATA$gene_info %>%
                dplyr::mutate(count=if_else(gene_list == listname,paste("n =",n_distinct(new_gene_match$gene, na.rm = T)),
                                            count))
            }
          }
          setProgress(i/length(meta_data$filepath), 
                      detail = paste("Updating gene lists",meta_data$nick[i]))
          # add missing data
          LIST_DATA$gene_info <<- LIST_DATA$gene_info %>%
            dplyr::filter(gene_list == names(LIST_DATA$gene_file)[gg]) %>%
            dplyr::mutate(set = meta_data$nick[i],
                          count = count[i],
                          mycol = meta_data$color[i],
                          onoff = "0",
                          sub = " ",
                          plot_set = " ") %>%
            bind_rows(LIST_DATA$gene_info, .)
        }
    }
    } else {
      return()
      }
  })
    if(is_empty(meta_data$filepath)){
      return()
    }
    # first time starting
    if (LIST_DATA$STATE[1] == 0) {
      shinyjs::show("startoff")
      shinyjs::show("startoff2")
      # sets lines and labels
      # tries to guess lines and labels type
      num_bins <- LIST_DATA$x_plot_range[2]
      if (num_bins == 80 & LIST_DATA$STATE[3] == '543') {
        reactive_values$setsliders <- c(1,80,14,18,19,45)
        LIST_DATA$STATE[3] <<- kLinesandlabels[1]
      } else if (num_bins == 80 & LIST_DATA$STATE[3] == '5') {
        reactive_values$setsliders <- c(1,80,20,60,0,0)
        LIST_DATA$STATE[3] <<- kLinesandlabels[2]
      } else if (num_bins == 2 & LIST_DATA$STATE[3] == 'PI') {
        reactive_values$setsliders <- c(1,2,1,1,2,2)
        LIST_DATA$STATE[3] <<- kLinesandlabels[3]
      } else if (num_bins == 205 & LIST_DATA$STATE[3] == '5') {
        reactive_values$setsliders <- c(1,205,1,15,0,0)
        LIST_DATA$STATE[3] <<- kLinesandlabels[4]
      } else if (LIST_DATA$STATE[3] == '5') {
        reactive_values$setsliders <- c(1,num_bins,1,
                                        num_bins/2,0,0)
      } else if (LIST_DATA$STATE[3] == '4') {
        reactive_values$setsliders <- c(1,num_bins,1,num_bins,0,0)
        LIST_DATA$STATE[3] <<- kLinesandlabels[6]
      } else if (LIST_DATA$STATE[3] == '3') {
        reactive_values$setsliders <- c(1,num_bins,
                                        num_bins/2,
                                        num_bins,0,0)
        LIST_DATA$STATE[3] <<- kLinesandlabels[7]
      } else {
        reactive_values$setsliders <- c(
        1, num_bins, 
        floor(num_bins / 5.5),
        floor(num_bins / 4.4),
        floor(num_bins / 4.4) + 1,
        floor(num_bins / 1.77))
        LIST_DATA$STATE[3] <<- kLinesandlabels[8]
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
    removeCssClass(selector = "a[data-value='filenorm']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='genelists']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='sorttool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='ratiotool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='clustertool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='cdftool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='DataTableTool']", class = "inactiveLink")
    
    gts <- LIST_DATA$table_file %>% group_by(set) %>%
      summarise(number_of_genes = n_distinct(gene, na.rm = T)) %>%
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
      gg <- LIST_DATA$gene_info %>% dplyr::filter(!str_detect(gene_list,"^Complete")) %>%
        select(., gene_list, count) %>% dplyr::rename(Usable = count) %>% 
        distinct()
      ggg <- NULL
      for (i in gg$gene_list) {
        ggg <-
          c(
            ggg,
            sapply(LIST_DATA$gene_file[i], "[[", "full") %>% bind_cols(.) %>% suppressMessages() %>% n_distinct(1)
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
    print("load gene file")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   # load info, update select boxes, switching works and changing info and plotting
                   LD <- LoadGeneFile(input$filegene1$datapath,
                                      input$filegene1$name,
                                      LIST_DATA)
                 })
    if (!is.null(LD)) {
      LIST_DATA <<- LD
    }
    shinyjs::reset("filegene1")
    
    gg <-
      LIST_DATA$gene_info %>% dplyr::filter(!str_detect(gene_list,"^Complete")) %>%
      select(., gene_list, count) %>% dplyr::rename(Usable = count) %>%
      distinct()
    ggg <- NULL
    for (i in gg$gene_list) {
      ggg <-
        c(
          ggg,
          sapply(LIST_DATA$gene_file[i], "[[", "full") %>% bind_cols(.) %>% suppressMessages() %>% n_distinct(1)
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
             ".txt")
    },
    content = function(file) {
      new_comments <-
        new_comments <- paste("#", Sys.Date(), "\n# File(s) used:")
      new_comments <-
        c(new_comments, paste("#", distinct(LIST_DATA$gene_info, set)$set))
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
          LIST_DATA$gene_file[[input$selectsave]]$full
      if(input$selectsave == "CDF Log2 PI Cumulative plot"){
        new_comments2 <- new_comments2 %>% select(-plot_set)
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
    print("plot button")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   list_data_frame <- Active_list_data(LIST_DATA,input$mygroup)
                   if (!is_empty(list_data_frame)) {
                     LIST_DATA$gene_info <<- rows_update(LIST_DATA$gene_info,
                                                         list_data_frame %>% 
                                                           distinct(set,gene_list,plot_set), 
                                                         by=c("set","gene_list"))
                     reactive_values$Apply_Math <- ApplyMath(
                       list_data_frame,
                       input$myMath,
                       input$selectplotnrom,
                       as.numeric(input$selectplotBinNorm),
                       input$mygroup
                     )
                     reactive_values$Y_Axis_numbers <-
                       YAxisValues(
                         reactive_values$Apply_Math,
                         input$sliderplotBinRange,
                         c(0,100),
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
                                                      input$mygroup)
                       } else {
                       LIST_DATA$ttest <<- NULL
                       }
                     }
                     reactive_values$Plot_Options <- NULL
                     reactive_values$Plot_Options <-
                       MakePlotOptionFrame(LIST_DATA$gene_info)
                     
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
    print("actionMathUpDatePlot")
    mymath <- c(input$myMath,
                input$selectplotnrom,
                input$selectplotBinNorm,
                input$checkboxsmooth,
                input$checkboxlog2,
                input$sliderplotBinRange,
                input$mygroup,
                input$checkboxauc
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
          MakePlotOptionFrame(LIST_DATA$gene_info)
      }
    }
    updateBoxSidebar(id = "sidebarmath")
  })
  
  # reactive Apply_Math, sets Y axis min max ----
  observeEvent(reactive_values$Y_Axis_numbers_set, ignoreInit = T, {
    print("updates reactive_values$Y_Axis_numbers")
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
  
  # reactive Plot_Options, triggers plot ----
  observeEvent(reactive_values$Plot_Options, ignoreInit = T, ignoreNULL = T, {
    print("plot")
    if(!is.null(LIST_DATA$ttest)){
      mm <- reactive_values$ttest_options[1:2]
      Plot_Options_ttest <- MakePlotOptionttest(LIST_DATA$ttest, as.numeric(mm),
                                                 input$selectttestlog,input$hlinettest,input$padjust,input$switchttesttype)
    }
    Y_Axis_Label <- YAxisLabel(input$myMath,
                               input$selectplotnrom,
                               as.numeric(input$selectplotBinNorm),
                               input$checkboxsmooth,
                               input$checkboxlog2)
    Lines_Labels_List <- reactive_values$Lines_Labels_List
    # make line smaller so group ribbon stands out more
    if(input$mygroup == "groups only"){
      Lines_Labels_List$mysize[2] <- Lines_Labels_List$mysize[2] / 2
    } 
    reactive_values$Plot_controler <-
      GGplotLineDot(
        reactive_values$Apply_Math,
        input$sliderplotBinRange,
        reactive_values$Plot_Options,
        reactive_values$Y_Axis_numbers,
        Lines_Labels_List,
        input$checkboxsmooth, Plot_Options_ttest,
        input$checkboxlog2,
        Y_Axis_Label,
        input$sliderplotOccupancy,
        input$checkboxauc
      )
    LIST_DATA$STATE[2] <<- 1
  })
  
  # checks that number of names == position ----
  observeEvent(c(input$landlnames, input$landlposition), ignoreInit = TRUE, {
    if(LIST_DATA$STATE[2] > 0){
      print("names == position")
      my_pos <-
        suppressWarnings(as.numeric(unlist(
          strsplit(input$landlposition, split = "\\s+")
        )))
      my_label <- unlist(strsplit(input$landlnames, split = "\\s+"))
      if (any(is.na(my_pos))) {
        my_pos <- my_pos[is.na(my_pos)]
        updateTextInput(session, "landlposition", value = my_pos)
      }
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
    print("dropcolor")
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
                colourInput("colourhex", "Select color HEX",value = distinct(LIST_DATA$gene_info,mycol)$mycol[1]),
                tags$hr(),
                textInput("textrgbtohex", "RGB", value = RgbToHex(x = distinct(LIST_DATA$gene_info,mycol)$mycol[1], convert = "rgb")),
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
                        choices = distinct(LIST_DATA$gene_info,gene_list)$gene_list,
                        selected = distinct(LIST_DATA$gene_info,gene_list)$gene_list[1],
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$gene_info,gene_list)$gene_list)
                        )),
            pickerInput("selectdataoption", "", choices = distinct(LIST_DATA$gene_info,set)$set, 
                        selected = distinct(LIST_DATA$gene_info,set)$set[1],
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$gene_info,set)$set)
                        )),
            tags$hr(style = "color: #2e6da4; background-color: #2e6da4; border-color: #2e6da4;"),
            textInput("textnickname", "Update Nickname",value = distinct(LIST_DATA$gene_info,set)$set[1]),
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
    print("actioncdfcolor")
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
                colourInput("colourhex", "Select color HEX",value = distinct(LIST_DATA$gene_info %>% 
                                                                               dplyr::filter(str_detect(gene_list,"^CDF")),mycol)$mycol[1]),
                tags$hr(),
                textInput("textrgbtohex", "RGB", value = RgbToHex(x = distinct(LIST_DATA$gene_info %>% 
                                                                                 dplyr::filter(str_detect(gene_list,"^CDF")),mycol)$mycol[1], convert = "rgb")),
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
            pickerInput("selectgenelistoptions", "", width = 300, choices = "CDF Log2 PI Cumulative plot",
                        selected = "CDF Log2 PI Cumulative plot"),
            pickerInput("selectdataoption", "", choices = distinct(LIST_DATA$gene_info,set)$set,
                        selected = distinct(LIST_DATA$gene_info,set)$set[1],
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$gene_info,set)$set)
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
                 my_sel <- LIST_DATA$gene_info %>% 
                   dplyr::filter(gene_list == input$selectgenelistoptions & 
                                   set == input$selectdataoption)
                 print("options update")
                 updateColourInput(session, "colourhex", value = paste(my_sel$mycol))
                 
                 updateTextInput(session,
                                 "textnickname",
                                 value = paste(my_sel$set))
               })
  
  # update color based on rgb text input ----
  observeEvent(input$actionmyrgb, {
    print("color rgb")
    updateColourInput(session, "colourhex", value = RgbToHex(input$textrgbtohex, convert = "hex"))
  })
  
  # save color selected and update plot ----
  observeEvent(input$colourhex, ignoreInit = TRUE, {
    print("update text color")
    updateTextInput(session,
                    "textrgbtohex",
                    value = RgbToHex(x = input$colourhex, convert = "rgb"))
    if (!is.null(names(LIST_DATA$gene_file))) {
      my_sel <- LIST_DATA$gene_info %>% 
        dplyr::filter(gene_list == input$selectgenelistoptions & 
                        set == input$selectdataoption)
      if (input$colourhex != my_sel$mycol) {
        print("color new")
        LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
          dplyr::mutate(mycol=if_else(gene_list == input$selectgenelistoptions & 
                                        set == input$selectdataoption,
                                      input$colourhex, mycol))
        reactive_values$Picker_controler <- 
          c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, mycol)$mycol)
        if(my_sel$onoff != 0){
          if (!is.null(reactive_values$Apply_Math) & input$leftSideTabs == "mainplot") {
            reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA$gene_info)
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
      print("new nickname")
      if (any(input$textnickname == distinct(LIST_DATA$gene_info, set)$set)) {
        updateTextInput(session,
                        "textnickname",
                        value = paste0(input$selectdataoption,"-",input$textnickname,
                                       "-dup"))
      }
      LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
        dplyr::mutate(set = if_else(set == input$selectdataoption,
                                    input$textnickname, set)) %>% 
        dplyr::mutate(onoff = if_else(onoff == input$selectdataoption,
                                      input$textnickname, onoff)) %>% 
        dplyr::mutate(plot_set = str_replace(LIST_DATA$gene_info$plot_set,
                                             paste0("^",input$selectdataoption),input$textnickname)) %>% 
        dplyr::mutate(group = str_replace(LIST_DATA$gene_info$group,
                                          input$selectdataoption,input$textnickname))
      
      LIST_DATA$table_file <<- LIST_DATA$table_file %>%
        dplyr::mutate(set = if_else(set == input$selectdataoption,
                                    input$textnickname, set))
      reactive_values$Picker_controler <- 
        c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
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
  observeEvent(c(input$droplinesandlabels, reactive_values$droplinesandlabels), ignoreInit = T, {
    print("droplinesandlabels")
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
          div(
            style = "padding-left:25px; display:inline-block; text-align:left !important;",
            selectInput(
              selectize = T,
              "selectlineslabels",
              width = "200px",
              label = "quick set lines and labels",
              choices = kLinesandlabels,
              selected = LIST_DATA$STATE[3]
            )
          ),
          column(12,
                 div(
                   style = "padding:2px; display:inline-block; text-align:center;",
                   numericInput(
                     "numerictss",
                     "TSS bin",
                     value = reactive_values$Lines_Labels_List$myset[3],
                     min = 0,
                     max = 100
                   )
                 ),
                 div(
                   style = "padding:2px; display:inline-block; text-align:center;",
                   textInput("numerictssname", value = "TSS", label = "TSS label",width = "60px"),
                 ),
                 div(
                   style = "padding:2px; display:inline-block;",
                   numericInput(
                     "numericbody1",
                     "5|4 bin",
                     value = reactive_values$Lines_Labels_List$myset[1],
                     min = 0,
                     max = 100
                   )
                 ),
                 div(
                   style = "padding:2px; display:inline-block; text-align:center;",
                   numericInput(
                     "numericbody2",
                     "4|3 bin",
                     value = reactive_values$Lines_Labels_List$myset[2],
                     min = 0,
                     max = 100
                   )
                 ),
                 div(
                   style = "padding:2px; display:inline-block; text-align:center;",
                   numericInput(
                     "numerictes",
                     "pA bin",
                     value = reactive_values$Lines_Labels_List$myset[4],
                     min = 0,
                     max = 100
                   )
                 ),
                 div(
                   style = "padding:2px; display:inline-block; text-align:center;",
                   textInput("numerictesname", value = "pA", label = " TES label",width = "60px"),
                 ),
                 div(
                   style = "padding:2px; display:inline-block; text-align:center;",
                   numericInput(
                     "numericbinsize",
                     "bp/bin",
                     value = reactive_values$Lines_Labels_List$myset[5],
                     min = 20,
                     max = 1000,
                     step = 5
                   )
                 ),
                 div(
                   style = "padding:2px 8px 2px 2px; display:inline-block; text-align:center;",
                   numericInput(
                     "numericlabelspaceing",
                     "every bin",
                     value = reactive_values$Lines_Labels_List$myset[6],
                     min = 0,
                     max = 100
                   )
                 )
          ),
          helpText("For 543 style 0 > TSS < 5|4 < 4|3 < pA < max bin"),
          div(
            textInput("landlnames", reactive_values$Lines_Labels_List$mylabels, label = "Yaxis labels"),
            textInput("landlposition", reactive_values$Lines_Labels_List$mybrakes, label = "Yaxis label position (numbers only)")
          ),
          helpText("select buttons for more options"),
          column(
            12,
            div(
              style = "padding-left: -5px; display:inline-block;",
              dropdownButton(
                tags$h3("Set TSS Options"),
                
                selectInput(
                  inputId = 'selecttsscolor',
                  label = 'TSS line and label color',
                  choices = c("red", "green", "blue", "brown", "black", "white"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_color[1]
                ),
                selectInput(
                  inputId = 'selecttssline',
                  label = 'TSS line type',
                  choices = c("dotted", "solid"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_type[1]
                ),
                icon = icon("sliders-h"),
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
                  choices = c("red", "green", "blue", "brown", "black", "white"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_color[3]
                ),
                selectInput(
                  inputId = 'selectbody1line',
                  label = '5|4 line type',
                  choices = c("dotted", "solid"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_type[3]
                ),
                icon = icon("sliders-h"),
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
                  choices = c("red", "green", "blue", "brown", "black", "white"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_color[4]
                ),
                selectInput(
                  inputId = 'selectbody2line',
                  label = '4|3 line type',
                  choices = c("dotted", "solid"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_type[4]
                ),
                icon = icon("sliders-h"),
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
                  choices = c("dotted", "solid"),
                  selected = reactive_values$Lines_Labels_List$myline$use_virtical_line_type[2]
                ),
                icon = icon("sliders-h"),
                status = "danger",
                tooltip = tooltipOptions(title = "TES Options")
              )
            ),
            div(
              style = "padding-left: 25px; display:inline-block;",
              dropdownButton(
                tags$h3("Set font Options"),
                
                numericInput(
                  inputId = 'selectvlinesize',
                  "Set vertcal line size",
                  value = reactive_values$Lines_Labels_List$mysize[1],
                  min = .5,
                  max = 10,
                  step = .5
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
                icon = icon("sliders-h"),
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
                  inputId = 'selectlegendsize',
                  "Set legend size",
                  value = reactive_values$Lines_Labels_List$mysize[5],
                  min = 1,
                  max = 20,
                  step = 1
                ),
                icon = icon("sliders-h"),
                status = "warning",
                tooltip = tooltipOptions(title = "Line Options")
              )
            )
          )
        ),
        actionButton("actionlineslabels", "SET and Plot"),
      )
    ))
  })
  
  # observe switching tabs ----
  observeEvent(input$leftSideTabs, ignoreInit = TRUE, {
    print("switch tab")
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
        c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
      if(LIST_DATA$STATE[1] == 0){
        LIST_DATA$STATE[1] <<- 1
        LIST_DATA$STATE[2] <<- -10
        shinyjs::hide("actionmyplotshow")
        updateSliderInput(session,"sliderplotBinRange",
                          min = reactive_values$setsliders[1],
                          max = reactive_values$setsliders[2],
                          value = LIST_DATA$x_plot_range)
        reactive_values$droplinesandlabels <- 1
      }
    }
    # sort/filter tab ----
    if (input$leftSideTabs == "sorttool"){
      if(!is.null(input$sortSamples)){ 
        if(input$sortSamples[1] == "select sample(s)"){
          updateSliderInput(
            session,
            "slidersortbinrange",
            min = reactive_values$setsliders[1],
            max = reactive_values$setsliders[2],
            value = reactive_values$setsliders[3:4]
          )
          updateSliderInput(
            session,
            "slidersortbinrangefilter",
            min = reactive_values$setsliders[1],
            max = reactive_values$setsliders[2],
            value = reactive_values$setsliders[c(5,2)]
          )
        }
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
      updatePickerInput(session, "sortSamples",
                        choices = c(distinct(LIST_DATA$gene_info, set)$set),
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$gene_info, set)$set)
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
    # file norm tab ----
    if (input$leftSideTabs == "filenorm") {
      updatePickerInput(
        session,
        "pickernumerator",
        choices = distinct(LIST_DATA$gene_info, set)$set,
        choicesOpt = list(style = paste("color", dplyr::select(
          dplyr::filter(LIST_DATA$gene_info,
                        gene_list == names(LIST_DATA$gene_file)[1]),
          mycol)$mycol, sep = ":"),
          content = gsub("(.{55})", "\\1<br>", distinct(LIST_DATA$gene_info, set)$set)
          )
      )
      updatePickerInput(
        session,
        "pickerdenominator",
        choices = c(distinct(LIST_DATA$gene_info, set)$set,"-1"),
        choicesOpt = list(style = paste("color", dplyr::select(
          dplyr::filter(LIST_DATA$gene_info,
                        gene_list == names(LIST_DATA$gene_file)[1]),
          mycol)$mycol, sep = ":"),
          content = gsub("(.{55})", "\\1<br>", c(distinct(LIST_DATA$gene_info, set)$set,"-1"))
          )
      )
      updatePickerInput(
        session,
        "pickergroupsample",
        choices = distinct(LIST_DATA$gene_info, set)$set,
        choicesOpt = list(
          content = gsub("(.{55})", "\\1<br>", distinct(LIST_DATA$gene_info, set)$set)
        )
      )
      output$valueboxnormfile <- renderValueBox({
        valueBox("0%",
                 "Done",
                 icon = icon("cogs"),
                 color = "yellow")
      })
      gts <- LIST_DATA$gene_info %>% 
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
          updateSliderInput(
            session,
            "sliderbinratio1",
            min = reactive_values$setsliders[1],
            max = reactive_values$setsliders[2],
            value = reactive_values$setsliders[3:4]
          )
          updateSliderInput(
            session,
            "sliderbinratio2",
            min = 0,
            max = reactive_values$setsliders[2],
            value = reactive_values$setsliders[5:6]
          )
          updateSliderInput(
            session,
            "sliderRatioBinNorm",
            min = 0,
            max = reactive_values$setsliders[2],
            value = 0
          )
        }
      } 
      updateSelectInput(
        session,
        "selectratiofile",
        choices = names(LIST_DATA$gene_file)
      )
      updatePickerInput(session, "pickerratio1file",
                        choices = c(distinct(LIST_DATA$gene_info, set)$set),
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$gene_info, set)$set)
                        ))
      updatePickerInput(session, "pickerratio2file",
                        choices = c("None", distinct(LIST_DATA$gene_info, set)$set),
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", c("None", distinct(LIST_DATA$gene_info, set)$set))
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
      if(!is.null(input$clusterSamples)){ 
        if(input$clusterSamples[1] == "select sample(s)"){
          shinyjs::disable("onoffdendrogram")
          shinyjs::hide("hideclusterplots1")
          shinyjs::hide("hideclustertable")
          shinyjs::hide("hideclusterplots2")
          updateSliderInput(
            session,
            "sliderbincluster",
            min = reactive_values$setsliders[1],
            max = reactive_values$setsliders[2],
            value = reactive_values$setsliders[3:4]
          )
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
                        choices = c(distinct(LIST_DATA$gene_info, set)$set),
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$gene_info, set)$set)
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
      if(input$sliderbincdf1[1] == 0){
      updateSliderInput(
        session,
        "sliderbincdf1",
        min = reactive_values$setsliders[1],
        max = reactive_values$setsliders[2],
        value = reactive_values$setsliders[3:4]
      )
      updateSliderInput(
        session,
        "sliderbincdf2",
        min = reactive_values$setsliders[1],
        max = reactive_values$setsliders[2],
        value = reactive_values$setsliders[5:6]
      )
      }
    shinyjs::hide('plotcdf')
    shinyjs::hide('plotcdfscatter')
    shinyjs::disable('actioncdfcolor')
    pickercdf <- list()
    for (i in names(LIST_DATA$gene_file)) {
      pickercdf[[i]] <-
        list(div(
          style = "margin-bottom: -10px;",
          pickerInput(
            inputId = gsub(" ", "-cdfspace2-", gsub("\n", "-cdfspace1-", i)),
            label = i,
            width = "99%",
            choices = distinct(LIST_DATA$gene_info, set)$set,
            multiple = T,
            options = list(
              `actions-box` = TRUE,
              `selected-text-format` = "count > 0"
            ),
            choicesOpt = list(style = paste("color", dplyr::select(
              dplyr::filter(LIST_DATA$gene_info,
                            gene_list == names(LIST_DATA$gene_file)[1]),
              mycol)$mycol, sep = ":"),
              content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$gene_info, set)$set))
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
    output$DynamicCDFPicker_main <- renderUI({
      pickercdf[!str_detect(names(LIST_DATA$gene_file),"^Filter|^Gene_List_|^Ratio_|^Cluster_|^CDF")]
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
                        choices = c(distinct(LIST_DATA$gene_info, set)$set),
                        selected = c(distinct(LIST_DATA$gene_info, set)$set)[1],
                        choicesOpt = list(
                          content = gsub("(.{35})", "\\1<br>", c(distinct(LIST_DATA$gene_info, set)$set))
                        ))
    }
    
  })
  
  # action button update lines and labels ----
  observeEvent(input$actionlineslabels, ignoreInit = TRUE, {
    print("action lines and labels")
    LIST_DATA$STATE[3] <<- input$selectlineslabels
    my_pos <-
      suppressWarnings(as.numeric(unlist(
        strsplit(input$landlposition, split = "\\s+")
      )))
    my_label <- unlist(strsplit(input$landlnames, split = "\\s+"))
    if (length(my_pos) == 0) {
      my_label <- "none"
      my_pos <- reactive_values$setsliders[2] * 2
    }
    
    # if tss or tes location make sure there is text
    if(nchar(trimws(input$numerictssname)) == 0 & input$numerictss > 0){
      updateTextInput(session, "numerictssname", value = "TSS")
    }
    if(nchar(trimws(input$numerictesname)) == 0 & input$numerictes > 0){
      updateTextInput(session, "numerictesname", value = "pA")
    }
    
    reactive_values$Lines_Labels_List <-
      LinesLabelsListPlot(
        input$numericbody1,
        input$selectbody1color,
        input$selectbody1line,
        input$numericbody2,
        input$selectbody2color,
        input$selectbody2line,
        input$numerictss,
        input$selecttsscolor,
        input$selecttssline,
        input$numerictes,
        input$selecttescolor,
        input$selecttesline,
        my_label,
        my_pos,
        input$selectvlinesize,
        input$selectlinesize,
        input$selectfontsizex,
        input$selectfontsizey,
        input$selectlegendsize,
        input$numericbinsize,
        input$numericlabelspaceing
      )
    removeModal()
  })
  
  # quick lines and labels preset change ----
  observeEvent(input$selectlineslabels, ignoreInit = TRUE, {
    if (input$selectlineslabels == "") {
      return()
    }
    print("quick Lines & Labels")
    myset <- LinesLabelsPreSet(input$selectlineslabels)
    updateNumericInput(session, "numericbody1", value = myset[1])
    updateNumericInput(session, "numericbody2", value = myset[2])
    updateNumericInput(session, "numerictss", value = myset[3])
    updateNumericInput(session, "numerictes", value = myset[4])
    updateNumericInput(session, "numericbinsize", value = myset[5])
    updateNumericInput(session, "numericlabelspaceing", value = myset[7])
    
  })
  
  # keep sizes real numbers lines and labels ----
  observeEvent(
    c(
      input$selectvlinesize,
      input$selectlinesize,
      input$selectfontsizex,
      input$selectfontsizey,
      input$selectlegendsize
    ),
    ignoreInit = TRUE,
    {
      if(LIST_DATA$STATE[2] > 0){
        print("keep bin positions in bounds > 0")
        mynum <- c(2, 2.5, 13, 13, 10)
        myset <- c(
          input$selectvlinesize,
          input$selectlinesize,
          input$selectfontsizex,
          input$selectfontsizey,
          input$selectlegendsize
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
          }
        }
      }
    })
  
  # Update lines and labels box's ----
  observeEvent(
    c(
      input$numericbody1,
      input$numericbody2,
      input$numerictss,
      input$numerictssname,
      input$numerictes,
      input$numerictesname,
      input$numericbinsize,
      input$numericlabelspaceing
    ),
    ignoreInit = TRUE,
    {
      if(LIST_DATA$STATE[2] > 0){
        print("observe line and labels")
        myset <- c(
          input$numericbody1,
          input$numericbody2,
          input$numerictss,
          input$numerictes,
          input$numericbinsize,
          input$numericlabelspaceing
        )
        # keep bin positions in bounds > 0
        for (i in seq_along(myset)) {
          if (is.na(myset[i]) | myset[i] < 0) {
            myset[i] <- 0
            updateNumericInput(session, "numericbody1", value = myset[1])
            updateNumericInput(session, "numericbody2", value = myset[2])
            updateNumericInput(session, "numerictss", value = myset[3])
            updateNumericInput(session, "numerictes", value = myset[4])
            updateNumericInput(session, "numericbinsize", value = myset[5])
            updateNumericInput(session, "numericlabelspaceing", value = myset[6])
          }
        }
        Lines_Labels_List <- LinesLabelsListset(myset[1],
                                                myset[2],
                                                myset[3],
                                                myset[4],
                                                myset[5],
                                                LIST_DATA$x_plot_range[2],
                                                myset[6],
                                                input$numerictssname,
                                                input$numerictesname)
        
        # set label and position numbers
        updateTextInput(session,
                        "landlnames",
                        value = paste(Lines_Labels_List$mylabels, collapse = " "))
        updateTextInput(session,
                        "landlposition",
                        value = paste(Lines_Labels_List$mybrakes , collapse = " "))
      } else {
        LIST_DATA$STATE[2] <<- 1
      }
    })
  
  # dropttest ----
  observeEvent(input$dropttest, ignoreInit = T, {
    if (is.null(reactive_values$ttest)){
      ttesttype <- "by files"
    } else {
      ttesttype <- reactive_values$ttest_values[1]
    }
    if(n_distinct(LIST_DATA$gene_info$gene_list, na.rm = T) > 1){
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
        print("t.test gets line colors")
        mycol <- LIST_DATA$ttest %>% dplyr::filter(set == input$selectttestitem) %>% distinct(mycol)
        updateColourInput(session, "selectcolorttest", value = paste(mycol))
      }
    }
  })
  
  # t.test updates line colors ----
  observeEvent(input$selectcolorttest,ignoreInit = T,{
    if(!is.null(LIST_DATA$ttest)){
      print("t.test updates line colors")
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
      print("t.test button")
      list_data_frame <- Active_list_data(LIST_DATA,input$mygroup)
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
                                                          input$mygroup)
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
          MakePlotOptionFrame(LIST_DATA$gene_info)
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
  output$plot1sort <- renderPlot({
    reactive_values$Plot_controler_sort_min
  })
  output$plot2sort <- renderPlot({
    reactive_values$Plot_controler_sort_max
  })
  output$plot1cluster <- renderPlot({
    reactive_values$Plot_controler_cluster
  })
  output$plot2cluster <- renderPlot({
    reactive_values$Plot_controler_dcluster
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
                   print("checkbox on/off")
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
                   LIST_DATA$gene_info <<-
                     CheckBoxOnOff(checkboxonoff,
                                   LIST_DATA$gene_info)
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
      print("Dynamic pickers update")
      
      pickerlist <- list()
      for (i in names(LIST_DATA$gene_file)) {
        pickerlist[[i]] <-
          list(div(
            style = "margin-bottom: -10px;",
            pickerInput(
              inputId = gsub(" ", "-bensspace2-", gsub("\n", "-bensspace1-", i)),
              label = i,
              width = "99%",
              choices = distinct(LIST_DATA$gene_info, set)$set,
              selected =  dplyr::select(dplyr::filter(LIST_DATA$gene_info, 
                                                      gene_list == i & onoff != 0), 
                                        onoff)$onoff,
              multiple = T,
              options = list(
                `actions-box` = TRUE,
                `selected-text-format` = "count > 0"),
              choicesOpt = list(style = paste("color", 
                                              dplyr::select(dplyr::filter(LIST_DATA$gene_info, 
                                                                          gene_list == i), mycol)$mycol,
                                              sep = ":"),
                                content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$gene_info, set)$set))
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
      if(any(str_detect(names(LIST_DATA$gene_file),"^CDF"))){
        output$DynamicGenePicker_cdf <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^CDF")]
        })
        shinyjs::show("showpickercdf")
      } else if(!any(str_detect(names(LIST_DATA$gene_file),"^CDF"))){
        shinyjs::hide("showpickercdf")
      }
      output$DynamicGenePicker_main <- renderUI({
        pickerlist[!str_detect(names(LIST_DATA$gene_file),"^Filter|^Gene_List_|^Ratio_|^Cluster_|^CDF")]
      })
      
    })
  
  # filter sum tool action ----
  observeEvent(input$actionsorttool, {
    print("sort tool")
    if (input$slidersortpercent < 50 &
        input$selectsorttop == "Middle%") {
      updateSliderInput(session, "slidersortpercent", value = 50)
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   reactive_values$Plot_controler_sort_min <- ggplot()
                   LD <- FilterTop(
                     LIST_DATA,
                     input$sortGeneList,
                     input$sortSamples,
                     input$slidersortbinrange[1],
                     input$slidersortbinrange[2],
                     input$slidersortpercent,
                     input$selectsorttop
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      LD <- LIST_DATA
      mylist <- last(grep("^Filter", names(LIST_DATA$gene_file)))
      LD$gene_info <- LD$gene_info %>%
        dplyr::mutate(onoff=if_else(gene_list == names(LD$gene_file)[mylist] &
                                      set %in% input$sortSamples, set, "0"))
      list_data_frame <- Active_list_data(LD)
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
    print("sort % tool")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a bit ...',
                 value = 0,
                 {
    reactive_values$Plot_controler_sort_min <- ggplot()
    sortmin <- FilterPer(LIST_DATA, 
                         input$sortGeneList,
                         input$sortSamples,
                         input$slidersortbinrange,
                         c(input$numericsortmin,input$numericsortmax),
                         input$selectsortper)
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
  
  # filter peak tool action ----
  observeEvent(input$actionsortpeak, ignoreInit = TRUE, {
    print("sort Peak")
    shinyjs::hide("hidesortplots1")
    shinyjs::hide("hidesortplots2")
    if (any(between(input$slidersortbinrange,
                    input$slidersortbinrangefilter[1],
                    input$slidersortbinrangefilter[2]))) {
      showModal(modalDialog(
        title = "Information message",
        paste("Bins regions should not overlap, \nBins set to default"),
        size = "s",
        easyClose = TRUE
      ))
      updateSliderInput(
        session,
        "slidersortbinrange",
        min = reactive_values$setsliders[1],
        max = reactive_values$setsliders[2],
        value = reactive_values$setsliders[3:4]
      )
      updateSliderInput(
        session,
        "slidersortbinrangefilter",
        min = reactive_values$setsliders[1],
        max = reactive_values$setsliders[2],
        value = reactive_values$setsliders[c(5,2)]
      )
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a bit ...',
                 value = 0,
                 {
    sortmin <- FilterPeak(LIST_DATA, 
                         input$sortGeneList,
                         input$sortSamples,
                         input$slidersortbinrange,
                         input$slidersortbinrangefilter,
                         input$selectsortpeak)
                 })
    
    if(!is.null(sortmin)){
      LIST_DATA <<- sortmin
      # pull info for preview plot
      mylist <- last(grep("^Filter", names(sortmin$gene_file),value = T))
      sortmin$gene_info <- sortmin$gene_info %>%
        dplyr::mutate(onoff=if_else(gene_list == mylist &
                                      set %in% input$sortSamples, set, "0"))
      list_data_frame <- Active_list_data(sortmin)
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
  
  # observe norm gene pickers ----
  observeEvent(c(input$pickernumerator, input$adddata,
                 input$pickerdenominator), {
                   if (input$pickernumerator != "") {
                     updateTextInput(session, "textnromname",value = paste(input$pickernumerator, input$adddata,
                                                                           input$pickerdenominator))
                     output$valueboxnormfile <- renderValueBox({
                       valueBox("0%",
                                "Done",
                                icon = icon("cogs"),
                                color = "yellow")
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
                     input$adddata,
                     input$textnromname
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      updatePickerInput(session,
                        "pickernumerator", selected = "",
                        choices = distinct(LIST_DATA$gene_info, set)$set,
                        choicesOpt = list(style = paste("color", 
                                                        dplyr::select(
                                                          dplyr::filter(LIST_DATA$gene_info,
                                                                        gene_list == names(
                                                                          LIST_DATA$gene_file)[1]), 
                                                          mycol)$mycol, 
                                                        sep = ":"),
                                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$gene_info, set)$set)
                                          ))
      updatePickerInput(session,
                        "pickerdenominator", selected = "",
                        choices = distinct(LIST_DATA$gene_info, set)$set,
                        choicesOpt = list(style = paste("color", 
                                                        dplyr::select(
                                                          dplyr::filter(LIST_DATA$gene_info,
                                                                        gene_list == names(
                                                                          LIST_DATA$gene_file)[1]),
                                                          mycol)$mycol, 
                                                        sep = ":"),
                                          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$gene_info, set)$set)
                                          ))
      updatePickerInput(
        session,
        "pickergroupsample",
        choices = distinct(LIST_DATA$gene_info, set)$set,
        choicesOpt = list(
          content = gsub("(.{35})", "\\1<br>", distinct(LIST_DATA$gene_info, set)$set)
        )
      )
      gts <- LIST_DATA$gene_info %>% 
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
      ff <- distinct(LIST_DATA$gene_info, set)$set
      updateSelectInput(session,
                        "selectdataoption",
                        choices = ff)
    } else {
      #no new data file created
      output$valueboxnormfile <- renderValueBox({
        valueBox("0%",
                 "Done",
                 icon = icon("cogs"),
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
    LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
      mutate(group = if_else(set %in% input$pickergroupsample,
                             gsub("(.{20})", "\\1\n",
                                  gsub(",",":",input$textgroupname)), group))
      updatePickerInput(session,
                        "pickergroupsample", selected = "",
                        options = list(title = "Select at least 2 files"))
      updateTextInput(session, "textgroupname", value = "")
      gts <- LIST_DATA$gene_info %>% 
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
    print("gene lists action")
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
    print("generiate gene lists table")
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
    print("cluster tool action")
    shinyjs::hide('plot1cluster')
    shinyjs::hide('plot2cluster')
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
                       input$sliderbincluster[1],
                       input$sliderbincluster[2]
                     )
                 })
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
  observeEvent(c(input$selectclusternumber, reactive_values$clustergroups),
               ignoreInit = TRUE, ignoreNULL = TRUE,
               {
                 print("cluster tool number")
                 if (is.null(reactive_values$clustergroups)) {
                   return()
                 }
                 shinyjs::hide('hideclusterplots1')
                 shinyjs::hide("hideclustertable")
                 shinyjs::hide('hideclusterplots2')
                 withProgress(message = 'Calculation in progress',
                              detail = 'This may take a while...',
                              value = 0,
                              {
                                LD <-
                                  ClusterNumList(
                                    LIST_DATA,
                                    input$clusterGeneList,
                                    input$clusterSamples,
                                    input$sliderbincluster[1],
                                    input$sliderbincluster[2],
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
                   LD$gene_info <- LD$gene_info %>%
                     dplyr::mutate(onoff=if_else(str_detect(gene_list,"^Cluster_") &
                                                   set == input$clusterSamples, set, "0"))
                   withProgress(message = 'Calculation in progress',
                                detail = 'This may take a while...',
                                value = 0,
                                {
                                  list_data_frame <- Active_list_data(LD)
                                  if (!is_empty(list_data_frame)) {
                                    
                                    Apply_Cluster_Math <- ApplyMath(
                                      list_data_frame,
                                      "mean",
                                      "none",
                                      0
                                    )
                                  }
                                  reactive_values$Plot_controler_cluster <- ggplot()
                                  reactive_values$Plot_controler_dcluster <- ggplot()
                                  gp1 <-
                                    ggplot(Apply_Cluster_Math ,aes(as.numeric(bin),value,color=gene_list)) +
                                    geom_line() +
                                    ylab("Mean bin value") +
                                    theme(legend.position="bottom",
                                          legend.title = element_blank(),
                                          axis.title.x=element_blank())
                                  print(gp1)
                                  reactive_values$Plot_controler_cluster <- gp1
                                })
                   shinyjs::enable("onoffdendrogram")
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
  
  # plots deprogram of Clusters ----
  observeEvent(input$actiondclustertool, ignoreInit = TRUE,{
    shinyjs::show('hideclusterplots2')
    shinyjs::show('plot2cluster')
    gp2 <-
      suppressWarnings(fviz_dend(LIST_DATA$clust$cm, k = as.numeric(input$selectclusternumber),
                show_labels = F))
    print(gp2)
    reactive_values$Plot_controler_dcluster <- gp2
  })
  
  # Cluster reset controller ----
  observeEvent(c(input$clusterGeneList,
                 input$clusterSamples), ignoreNULL = TRUE, ignoreInit = TRUE, {
    reactive_values$clustergroups <- NULL
    shinyjs::disable("onoffdendrogram")
    shinyjs::hide("hideclusterplots1")
    shinyjs::hide("hideclustertable")
    shinyjs::hide("hideclusterplots2")
  })
  
  # Ratio tool action ----
  observeEvent(input$actionratiotool, ignoreInit = TRUE, {
    print("ratio tool action")
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
    if (input$sliderbinratio2[1] == 0 & input$sliderbinratio2[2] > 0) {
      updateSliderInput(session,
                        "sliderbinratio2",
                        value = c(0,0))
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
                       input$sliderbinratio1[1],
                       input$sliderbinratio1[2],
                       input$sliderbinratio2[1],
                       input$sliderbinratio2[2],
                       numericratio,
                       input$checkratiozero,
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
        my_range <- range(LIST_DATA$boxRatio$Ratio,na.rm = T) 
        updateNumericInput(session, "textboxmaxratio",
                           value = my_range[2])
        updateNumericInput(session, "textboxminratio",
                           value = my_range[1])
      } else {
        updateNumericInput(session, "textboxmaxratio",
                           value = 0)
        updateNumericInput(session, "textboxminratio",
                           value = 0)
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
  observeEvent(c(input$textboxmaxratio, input$textboxminratio),ignoreInit = TRUE, ignoreNULL = TRUE,{
                   if(!is.null(LIST_DATA$boxRatio)){
                     my_range <- c(floor(input$textboxminratio), ceiling(input$textboxmaxratio)) 
                     gb <- ggviolin(LIST_DATA$boxRatio, x= "set", y = "Ratio", fill="set",
                                       color="set",add = "boxplot", add.params = list(fill = "white"),xlab = "") + 
                         theme(legend.position = 'none') +
                       coord_cartesian(ylim = my_range)
                     print(gb)
                     reactive_values$Plot_controler_ratio <- gb 
                   } 
                   
                 })
  
  # CDF tool action ----
  observeEvent(input$actioncdftool, ignoreInit = TRUE, {
    print("CDF tool action")
    shinyjs::hide('plotcdf')
    shinyjs::hide('plotcdfscatter')
    if (any(between(
      input$sliderbincdf1,
      input$sliderbincdf2[1],
      input$sliderbincdf2[2]
    )) |
    any(between(
      input$sliderbincdf2,
      input$sliderbincdf1[1],
      input$sliderbincdf1[2]
    ))) {
      showModal(modalDialog(
        title = "Information message",
        paste("Bins regions should not overlab, \nBins set to 1/3 2/3"),
        size = "s",
        easyClose = TRUE
      ))
      updateSliderInput(session,
                        "sliderbincdf1",
                        value = c(
                          LIST_DATA$x_plot_range[1],
                          floor(LIST_DATA$x_plot_range[2] / 4)
                        ))
      updateSliderInput(session,
                        "sliderbincdf2",
                        value = c(
                          ceiling(LIST_DATA$x_plot_range[2] / 4) + 1,
                          LIST_DATA$x_plot_range[2]
                        ))
    }
    ttt <-
      reactiveValuesToList(input)[gsub(" ", "-cdfspace2-", gsub("\n", "-cdfspace1-", names(LIST_DATA$gene_file)))]
    checkboxonoff <- list()
    for (i in names(ttt)) {
      for (tt in ttt[i]) {
        selectgenelistonoff <-
          gsub("-cdfspace2-", " ", gsub("-cdfspace1-", "\n", i))
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
                       input$sliderbincdf1[1],
                       input$sliderbincdf1[2],
                       input$sliderbincdf2[1],
                       input$sliderbincdf2[2]
                     )
                 })
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
      rr <- range(LIST_DATA$gene_file[[newname]]$full$value)
      updateSliderInput(session, "sliderrangecdf",
                        min = floor(rr[1]),
                        max = ceiling(rr[2]),
                        value = c(0,0))
      updateSliderInput(session, "sliderrangecdf",
                        min = floor(rr[1]),
                        max = ceiling(rr[2]),
                        value = c(floor(rr[1]), ceiling(rr[2])))
    } else {
      shinyjs::disable('actioncdfcolor')
      return()
    }
  })
  
  # CDF x plot range ----
  observeEvent(c(input$sliderrangecdf, reactive_values$df_options), ignoreInit = TRUE, {
    print("cdf plot observe range")
    newname <-
      grep("CDF ", names(LIST_DATA$gene_file), value = TRUE)
    if(is_empty(newname)){
      return()
    }
    df_options <-
      LIST_DATA$gene_info %>%
      dplyr::filter(gene_list ==  newname) %>%
      dplyr::mutate(set = paste(
        count,
        sub(" - ","\n", plot_set),sep = "\n")
      )
    df <- LIST_DATA$gene_file[[newname]]$full %>%
      full_join(.,df_options %>% select(set,plot_set),by="plot_set") %>%
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
    mycdf <- GGplotC(df, df_options, use_header,as.numeric(input$sliderrangecdf))
    output$plotcdf <- renderPlot({
      mycdf
    })
    output$plotcdfscatter <- renderPlot({
      ggscatter(df %>% arrange(gene,value), y = "value",x="gene",color = "plot_set",
                alpha=.8, palette = df_options$mycol) + 
        rremove("x.text") + rremove("legend")
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
        filter = "none",
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
  
  # QC action ----
  observeEvent(input$buttonFilterzero, ignoreInit = TRUE, {
    print("QC % tool") 
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
      tags$style(".inactiveLink {
                            pointer-events: none;
                           cursor: default;
                           }",
                 ".shiny-notification {
                         position: fixed;
                         top: 50%;
                         left: 50%;
                         transform: translate(-50%, -50%);
                         color: purple;
                         font-size: 20px;
                         font-style: italic;}"
                 )
    ),
    # disables tabs on start
    sidebarMenu(
      id = "leftSideTabs",
      menuItem("Load Data", tabName = "loaddata", icon = icon("file-import")),
      menuItem("Plot", tabName = "mainplot", icon = icon("chart-area")),
      menuItem("QC/Options", tabName = "qcOptions", icon = icon("clipboard-check")),
      menuItem("Norm data", tabName = "filenorm", icon = icon("copy")),
      menuItem("Compare Lists", tabName = "genelists", icon = icon("grip-lines")),
      menuItem("Filter Tool", tabName = "sorttool", icon = icon("filter")),
      menuItem("Ratio Tool", tabName = "ratiotool", icon = icon("percentage")),
      menuItem("Cluster Tools", tabName = "clustertool", icon = icon("object-group")),
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
                    title = "Load .table/URL.txt file",
                    width = 6,
                    style = "height: 150px;",
                    align = "center",
                    fileInput(
                      "filetable",
                      width = "75%",
                      label = "",
                      accept = c('.table','.url.txt'),
                      multiple = FALSE
                    ),
                    helpText("load windowed bedGraph file(s)"),
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
                  title = "Load Gene list, .txt/.bed",
                  width = 6,
                  style = "height: 150px;" ,
                  solidHeader = TRUE,
                  status = "navy",
                  align = "center",
                  fileInput(
                    "filegene1",
                    width = "75%",
                    label = "",
                    accept = c('.txt', 'bed'),
                    multiple = FALSE
                  ),
                  helpText("load gene list"),
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
            width = 50,
            tags$style(".calculator {color:#FF0000}"),
            icon = icon("calculator", class = "calculator"),
            fluidRow(
              column(
                4,
                selectInput("myMath",
                            label = "Math",
                            choices = c("mean", "sum", "median", "var"),
                            selected = "mean"
                )),
              column(
                4,
                selectInput(
                  "selectplotnrom",
                  label = "Y Normalization",
                  choices = c("none", "relative frequency", "rel gene frequency"),
                  selected = "none"
                )),
              column(
                3,
                selectInput(
                  "selectplotBinNorm",
                  label = "Bin Norm:",
                  choices = c(0:80),
                  selected = 0
                )),
              column(
                3,
                awesomeCheckbox("checkboxsmooth", label = "smooth")),
              column(
                2,
                awesomeCheckbox("checkboxlog2", label = "log2"),
              awesomeCheckbox("checkboxauc", label = "AUC")),
              column(
                3,
                numericInput("numericYRangeLow", label = "Plot Y min:", value = 0)
              ),
              column(
                3,
                numericInput("numericYRangeHigh", label = "Plot Y max:", value = 0)
              ),
              hidden(div(id = "hideplotgroup",
                column(
                  4,
                  selectInput("mygroup",
                              label = "plot group",
                              choices = c("none", "groups only", "groups and single"),
                              selected = "none"
                  ))
              )),
              column(
                11,
                sliderInput(
                  "sliderplotBinRange",
                  label = "Plot Bin Range:",
                  min = 0,
                  max = 80,
                  value = c(0, 80)
                )),
              column(
                6,offset = 4,
                actionBttn(
                  inputId = "actionMathUpDatePlot",
                  label = "Update Plot",
                  style = "unite",
                  color = "default",
                  size = "sm"
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
            label = "replace 0 with min/2",value = FALSE),
          valueBoxOutput("valueboxnormfile")
        ),
        box(
          title = "Select files for group mean",
          width = 12,
          status = "primary",
          solidHeader = T,
          collapsible = T,
          collapsed = T,
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
            width = 12, solidHeader = T,
            status = "primary",
            column(width = 6,
                   pickerInput("sortGeneList", label = "select list",
                               choices = (LIST_DATA$gene_info)),
                   div(
                     style = "margin-bottom: -20px;",
                     sliderInput(
                       "slidersortbinrange",
                       label = "Select Bin Range:",
                       min = 1,
                       max = 80,
                       value = c(1, 80)
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
                step = 0.1
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
            title = "filter peaks",
            solidHeader = T,
            width = 6,
            status = "navy",
            collapsible = T,
            fluidRow(column(
              12,
              style = "margin-bottom: -20px;",
              sliderInput(
                "slidersortbinrangefilter",
                label = "Select Bin Range of non-peak:",
                min = 1,
                max = 80,
                value = c(1, 80)
              )
            )
            ),
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
                label = "replace 0 with min/2",
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
                sliderInput(
                  "sliderbinratio1",
                  label = "Select numerator Bin Range:",
                  min = 1,
                  max = 80,
                  value = c(1, 1)
                )
              ),
              column(
                5,
                sliderInput(
                  "sliderbinratio2",
                  label = "Select denominator Bin Range:",
                  min = 0,
                  max = 80,
                  value = c(0, 80)
                )
              )
            ),
            sliderInput(
              "sliderRatioBinNorm",
              label = "Bin Norm:",
              min = 0,
              max = 80,
              value = 0
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
                         step = .5)
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
                                 choices = (LIST_DATA$gene_info)),
         pickerInput("clusterSamples", label = "select sample(s)",
                             choices = "select sample(s)",selected = "select sample(s)",
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
             sliderInput(
               "sliderbincluster",
               label = "Select Bin Range:",
               min = 1,
               max = 80,
               value = c(1, 80)
             )
           )),
           column(width = 6,
           actionButton("actionclustertool", "Get clusters")),
           div(
             id = "onoffdendrogram",
             column(width = 6,
           actionButton("actiondclustertool", "plot dendrogram")
           ))
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
                )),
         div(
           id = "hideclusterplots2",
           box(headerBorder = F,
               width = 12,
               withSpinner(plotOutput("plot2cluster",height = "200px"), type = 4)
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
                  collapsed = F,
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
                  collapsed = F,
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
                  collapsed = F,
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
                  collapsed = F,
                  uiOutput("DynamicCDFPicker_clusters")
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
              sliderInput(
                "sliderbincdf1",
                label = "Select numerator Bin Range:",
                min = 0,
                max = 80,
                value = c(0, 0)
              )),
              column(
                width = 6,
              sliderInput(
                "sliderbincdf2",
                label = "Select denominator Bin Range:",
                min = 0,
                max = 80,
                value = c(0, 80)
              )),
            sliderInput(
              "sliderrangecdf",
              label = "Select plot Range:",
              min = 0,
              max = 1,
              value = c(0, 1)
            ),
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
  controlbar = dashboardControlbar(),
  title = "DashboardPage"
)

# execute ----
shinyApp(ui = ui, server = server)