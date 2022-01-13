# Created by Benjamin Erickson BBErickson@gmail.com

source("R_scripts/helpers.R", local = TRUE)

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
    "DT",
    "patchwork",
    "zip",
    "ggpubr",
    "fastcluster"
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
  # t.test results $use is for numbers $gene_info for holding plotting options
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
  
  # reactive values ----
  reactive_values <- reactiveValues(
    Apply_Math = NULL,
    Plot_Options = NULL,
    Y_Axis_numbers = c(0,100),
    Lines_Labels_List = list(mybrakes="",mylabels=""),
    Picker_controler = NULL,
    mymath = c("mean", "none", "0", "FALSE", "FALSE", "0", "80"),
    ttest = NULL,
    ttest_values = c("none", "wilcox.test", "two.sided", "FALSE", "FALSE", "-log10", "fdr"),
    ttest_options = c(0, 0, 1, "select sample", 0.05)
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
  
  # sidebar observe
  
  # loads data file(s) ----
  observeEvent(input$filetable, {
    print("load file")
    shinyjs::disable("startoff")
    shinyjs::disable("startoff2")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- LoadTableFile(input$filetable$datapath,
                                       input$filetable$name,
                                       LIST_DATA)
                 })
    
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
    } else {
      showModal(modalDialog(
        title = "Information message",
        paste("No files loaded"),
        size = "s",
        easyClose = TRUE
      ))
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
        LIST_DATA$STATE[3] <<- kLinesandlabels[1]
      } else if (num_bins == 80 & LIST_DATA$STATE[3] == '5') {
        LIST_DATA$STATE[3] <<- kLinesandlabels[2]
      } else if (num_bins <= 60 & LIST_DATA$STATE[3] == '543') {
        LIST_DATA$STATE[3] <<- kLinesandlabels[3]
      } else if (num_bins == 205 & LIST_DATA$STATE[3] == '5') {
        LIST_DATA$STATE[3] <<- kLinesandlabels[4]
      } else if (LIST_DATA$STATE[3] == '5') {
        LIST_DATA$STATE[3] <<- kLinesandlabels[5]
      } else if (LIST_DATA$STATE[3] == '4') {
        LIST_DATA$STATE[3] <<- kLinesandlabels[6]
      } else if (LIST_DATA$STATE[3] == '3') {
        LIST_DATA$STATE[3] <<- kLinesandlabels[7]
      } else {
        LIST_DATA$STATE[3] <<- kLinesandlabels[8]
      }
    }
    # enables tabs after loading file
    shinyjs::enable("startoff")
    shinyjs::enable("startoff2")
    shinyjs::reset("filetable")
    removeCssClass(selector = "a[data-value='qcOptions']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='mainplot']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='filenorm']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='genelists']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='sorttool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='ratiotool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='clustertool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='cdftool']", class = "inactiveLink")
    
    gts <- LIST_DATA$table_file %>% group_by(set) %>%
      summarise(number_of_genes = n_distinct(gene)) %>%
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
      LIST_DATA$gene_file[["Complete"]]$use %>% summarise(total_number_distinct_genes = n_distinct(gene))
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
      gg <- LIST_DATA$gene_info %>% filter(gene_list != "Complete") %>%
        select(., gene_list, count) %>%
        distinct()
      ggg <- NULL
      for (i in names(LIST_DATA$gene_file)[-1]) {
        ggg <-
          c(
            ggg,
            sapply(LIST_DATA$gene_file[i], "[[", "full") %>% bind_cols(.) %>% suppressMessages() %>% n_distinct(1)
          )
      }
      ggg <- mutate(gg, "total_in_file" = ggg)
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
      LIST_DATA$gene_info %>% filter(gene_list != "Complete") %>%
      select(., gene_list, count) %>%
      distinct()
    ggg <- NULL
    for (i in names(LIST_DATA$gene_file)[-1]) {
      ggg <-
        c(
          ggg,
          sapply(LIST_DATA$gene_file[2], "[[", "full") %>% bind_cols(.) %>% suppressMessages() %>% n_distinct(1)
        )
    }
    ggg <- mutate(gg, "total_in_file" = ggg)
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
    updateSelectInput(session, "selectsave",
                      choices = names(LIST_DATA$gene_file))
    if( LIST_DATA$STATE[1] != 0){
      LIST_DATA$STATE[1] <<- .5
    }
  })
  
  # save functions, gene list ----
  output$downloadGeneList <- downloadHandler(
    filename = function() {
      paste0(input$selectsave, "_",
             Sys.Date(), "_",
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
          paste(LIST_DATA$gene_file[[input$selectsave]]$info)
        )))
      if (input$selectsave == "Complete") {
        # save $use
        new_comments2 <-
          LIST_DATA$gene_file[[input$selectsave]]$use$gene
      } else {
        # intersect $use with full and save values with comments on values
        new_comments2 <-
          LIST_DATA$gene_file[["1k_hg19_10k_isolate_bed"]]$full %>%
          semi_join(., LIST_DATA$gene_file[["1k_hg19_10k_isolate_bed"]]$use)
      }
      write_lines(new_comments, file)
      write_tsv(new_comments2,
                file,
                col_names = FALSE,
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
                   list_data_frame <- Active_list_data(LIST_DATA)
                   if (!is_empty(list_data_frame)) {
                     LIST_DATA$gene_info <<- rows_update(LIST_DATA$gene_info,
                                                         list_data_frame %>% 
                                                           distinct(set,gene_list,plot_set), 
                                                         by=c("set","gene_list"))
                     reactive_values$Apply_Math <- ApplyMath(
                       list_data_frame,
                       input$myMath,
                       input$selectplotnrom,
                       as.numeric(input$selectplotBinNorm)
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
                          n_distinct(list_data_frame$gene_list) == 1 | 
                          input$switchttest == "by files" & n_distinct(list_data_frame$set) == 1) {
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
                                                      input$selectttestpaired)
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
                input$sliderplotBinRange
    )
    Y_Axis_numbers <-
      c(input$numericYRangeLow,input$numericYRangeHigh)
    
    if(sum(reactive_values$mymath == mymath) != 7){
      reactive_values$mymath <- mymath
      reactive_values$myplot <- mymath
    } else if(sum(reactive_values$Y_Axis_numbers == Y_Axis_numbers) != 2){
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
    reactive_values$Plot_controler <-
      GGplotLineDot(
        reactive_values$Apply_Math,
        input$sliderplotBinRange,
        reactive_values$Plot_Options,
        reactive_values$Y_Axis_numbers,
        reactive_values$Lines_Labels_List,
        input$checkboxsmooth, Plot_Options_ttest,
        input$checkboxlog2,
        Y_Axis_Label,
        input$sliderplotOccupancy
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
  observeEvent(input$dropcolor, ignoreInit = T, {
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
                actionButton("actionmyrgb", "Update color",width = 100)
              )
            )
          ),
          box(
            title =  "Set Plot Options",
            width = 8,
            status = "navy",
            solidHeader = T,
            pickerInput("selectgenelistoptions", "", width = 300, choices = distinct(LIST_DATA$gene_info,gene_list)$gene_list,
                        selected = distinct(LIST_DATA$gene_info,gene_list)$gene_list[1]),
            pickerInput("selectdataoption", "", choices = distinct(LIST_DATA$gene_info,set)$set, 
                        selected = distinct(LIST_DATA$gene_info,set)$set[1]),
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
  
  # update display selected item info ----
  observeEvent(c(input$selectdataoption, input$selectgenelistoptions),
               ignoreInit = TRUE,
               {
                 if (LIST_DATA$STATE[1] == 0) {
                   return()
                 }
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
        
        my_sel <- LIST_DATA$gene_info %>% 
          dplyr::filter(gene_list == input$selectgenelistoptions) %>% 
          dplyr::select(mycol)
        reactive_values$Picker_controler <-
          paste("color", unlist(my_sel), sep = ":")
        
        if (!is.null(reactive_values$Apply_Math)) {
          reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA$gene_info)
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
        mutate(plot_set = str_replace(LIST_DATA$gene_info$plot_set,paste0("^",input$selectdataoption),input$textnickname))
      LIST_DATA$table_file <<- LIST_DATA$table_file %>%
        dplyr::mutate(set = if_else(set == input$selectdataoption,
                                    input$textnickname, set))
      reactive_values$myplot <- LIST_DATA$gene_info
      
      ff <- distinct(LIST_DATA$table_file, set)$set
      updatePickerInput(session,
                        "selectdataoption",
                        choices = ff,selected = input$textnickname)
      LIST_DATA$STATE[2] <<- -10
      reactive_values$Picker_controler <- 
        c(names(LIST_DATA$gene_file), distinct(LIST_DATA$table_file, set)$set)
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
          status = "navy",
          solidHeader = T,
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
                     value = 15,
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
                     value = 20,
                     min = 0,
                     max = 100
                   )
                 ),
                 div(
                   style = "padding:2px; display:inline-block; text-align:center;",
                   numericInput(
                     "numericbody2",
                     "4|3 bin",
                     value = 40,
                     min = 0,
                     max = 100
                   )
                 ),
                 div(
                   style = "padding:2px; display:inline-block; text-align:center;",
                   numericInput(
                     "numerictes",
                     "pA bin",
                     value = 45,
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
                     value = 100,
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
                     value = 5,
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
                  selected = "green"
                ),
                selectInput(
                  inputId = 'selecttssline',
                  label = 'TSS line type',
                  choices = c("dotted", "solid"),
                  selected = "dotted"
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
                  selected = "black"
                ),
                selectInput(
                  inputId = 'selectbody1line',
                  label = '5|4 line type',
                  choices = c("dotted", "solid"),
                  selected = "solid"
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
                  selected = "black"
                ),
                selectInput(
                  inputId = 'selectbody2line',
                  label = '4|3 line type',
                  choices = c("dotted", "solid"),
                  selected = "solid"
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
                  selected = "red"
                ),
                selectInput(
                  inputId = 'selecttesline',
                  label = 'TES line type',
                  choices = c("dotted", "solid"),
                  selected = "dotted"
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
                  value = 2,
                  min = .5,
                  max = 10,
                  step = .5
                ),
                numericInput(
                  inputId = 'selectfontsizex',
                  "Set X axis font size",
                  value = 13,
                  min = 1,
                  max = 30,
                  step = 1
                ),
                numericInput(
                  inputId = 'selectfontsizey',
                  "Set Y axis font size",
                  value = 13,
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
                  value = 2.5,
                  min = .5,
                  max = 10,
                  step = .5
                ),
                numericInput(
                  inputId = 'selectlegendsize',
                  "Set plot line size",
                  value = 10,
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
    if (input$leftSideTabs == "mainplot") {
      reactive_values$Picker_controler <- 
        c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
      if(LIST_DATA$STATE[1] == 0){
        LIST_DATA$STATE[1] <<- 1
        reactive_values$droplinesandlabels <- 1
      }
    }
    if (input$leftSideTabs == "sorttool"){
      if(input$sortSamples == "select sample(s)"){
        updateSliderInput(
          session,
          "slidersortbinrange",
          min = LIST_DATA$x_plot_range[1],
          max = LIST_DATA$x_plot_range[2],
          value = LIST_DATA$x_plot_range
        )
      }
      SGL <- input$sortGeneList
      updatePickerInput(session, "sortGeneList",
                        choices = c(distinct(LIST_DATA$gene_info, gene_list)$gene_list),
                        selected = SGL)
      STS <- input$sortSamples
      updatePickerInput(session, "sortSamples",
                        choices = c(distinct(LIST_DATA$gene_info, set)$set),
                        selected = STS)
      shinyjs::hide("hidesortplots1")
      shinyjs::hide("hidesortplots2")
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
      my_pos <- LIST_DATA$x_plot_range[2] * 2
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
        input$selectttestlinesize
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
    if(n_distinct(LIST_DATA$gene_info$gene_list) > 1){
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
        actionButton("actionttest","Apply")
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
      list_data_frame <- Active_list_data(LIST_DATA)
      if (!is_empty(list_data_frame)) {
        if (sum(ttest_values == reactive_values$ttest_values) != 7){
          reactive_values$ttest_values <- ttest_values
          if (input$switchttest != "none"){
            if(input$switchttest == "by lists" & n_distinct(list_data_frame$gene_list) == 1){
              showModal(modalDialog(
                title = "Information message",
                paste("Only 1 gene list active, replot with more than 1"),
                size = "s",
                easyClose = TRUE
              ))
              return()
            } else if (input$switchttest == "by files" & n_distinct(list_data_frame$set) == 1) {
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
                                                          input$selectttestpaired)
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
  output$plotcluster <- renderPlot({
    reactive_values$Plot_controler_cluster
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
                                              sep = ":"))
            )
          ))
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Filter"))){
        output$DynamicGenePicker_sort <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Filter")]
        })
        shinyjs::show("showpickersort")
      } else if(!any(str_detect(names(LIST_DATA$gene_file),"^Filter"))){
        shinyjs::hide("showpickersort")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Gene_List_"))){
        output$DynamicGenePicker_comparisons <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Gene_List_")]
        })
        shinyjs::show("showpickercomparisons")
      } else if(!any(str_detect(names(LIST_DATA$gene_file),"^Gene_List_"))){
        shinyjs::hide("showpickercomparisons")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Ratio_"))){
        output$DynamicGenePicker_ratio <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Ratio_")]
        })
        shinyjs::show("showpickerratio")
      } else if(!any(str_detect(names(LIST_DATA$gene_file),"^Ratio_"))){
        shinyjs::hide("showpickerratio")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Group_|^Cluster_"))){
        output$DynamicGenePicker_clusters <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Group_|^Cluster_")]
        })
        shinyjs::show("showpickercluster")
      } else if(!any(str_detect(names(LIST_DATA$gene_file),"^Group_|^Cluster_"))){
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
        pickerlist[!str_detect(names(LIST_DATA$gene_file),"^Filter|^Gene_List_|^Ratio_|^Group_|^Cluster_|^CDF")]
      })
      
    })
  
  # sort sum tool action ----
  observeEvent(input$actionsorttool, {
    print("sort tool")
    shinyjs::hide('actionsortdatatable')
    shinyjs::hide('sorttable')
    
    if (input$slidersortpercent < 50 &
        input$selectsorttop == "Middle%") {
      updateSliderInput(session, "slidersortpercent", value = 50)
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- FilterTop(
                     LIST_DATA,
                     input$sortGeneList,
                     input$sortSamples,
                     input$slidersortbinrange[1],
                     input$slidersortbinrange[2],
                     input$slidersortpercent,
                     input$selectsorttop,
                     input$checkboxfilterall
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      if(LIST_DATA$STATE[2] == 0){
        LIST_DATA$STATE[2] <<- -2
      }
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
      shinyjs::show('actionsortdatatable')
      if (any(grep("^Filter", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxsort <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$use$gene),
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
    if (!ol %in% names(LIST_DATA$gene_file)) {
      ol <- grep("^Filter", names(LIST_DATA$gene_file), value = TRUE)
      reactive_values$pickerfile_controler <- input$sortSamples
    } else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "sortGeneList",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
  })
  
  # sort % numeric controller ----
  observeEvent(c(input$numericsortmin,input$numericsortmax), ignoreInit = TRUE,{
    if (!is.numeric(input$numericsortmin)) {
      updateNumericInput(session, "numericsortmin", value = 1)
    }
    if (!is.numeric(input$numericsortmax)) {
      updateNumericInput(session, "numericsortmax", value = 99.5)
    }
    if (input$numericsortmin < 0 | input$numericsortmin > input$numericsortmax) {
      updateNumericInput(session, "numericsortmin", value = 1)
    }
    if (input$numericsortmax < input$numericsortmin | input$numericsortmax > 100) {
      updateNumericInput(session, "numericsortmax", value = 99.5)
    }
  })
  
  # sort min max between % tool action ----
  observeEvent(input$actionsortper, ignoreInit = TRUE, {
    print("sort % tool")
    
    sortmin <- FilterPer(LIST_DATA, 
                         input$sortGeneList,
                         input$sortSamples,
                         input$slidersortbinrange,
                         input$checkboxfilterall,
                         c(input$numericsortmin,input$numericsortmax),
                         input$selectsortper,
                         input$checkboxfilterall)
    
    if(!is_empty(sortmin$sortplot)){
      LIST_DATA <<- sortmin
      if(LIST_DATA$STATE[2] == 0){
        LIST_DATA$STATE[2] <<- -2
      } 
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
      if (any(grep("^Filter", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxsort <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$use$gene),
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
      reactive_values$Plot_controler_sort_min <-
        ggplot()
      reactive_values$Plot_controler_sort_max <-
        ggplot()
    }
    ol <- input$sortGeneList
    if (!ol %in% names(LIST_DATA$gene_file)) {
      ol <- grep("^Filter", names(LIST_DATA$gene_file), value = TRUE)
      reactive_values$pickerfile_controler <- input$sortSamples
    } else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "sortGeneList",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
  })
  
  # Filter gene list show data table ----
  observeEvent(input$actionsortdatatable, ignoreInit = TRUE, {
    print("show data table")
    if (any(grep("^Filter", names(LIST_DATA$gene_file)) > 0)) {
      newnames <-
        gsub("(.{20})", "\\1... ", names(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$full))
      dt <- datatable(
        rowid_to_column(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$full),
        rownames = FALSE,
        colnames = c("ID", strtrim(newnames, 30)),
        class = 'cell-border stripe compact',
        filter = 'top',
        caption = LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$info,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          scrollY = TRUE,
          autoWidth = FALSE,
          width = 5,
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
        
      ) %>% formatPercentage(names(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$full)[-1])
    } else {
      dt <- datatable(
        LIST_DATA$gene_file[[1]]$empty,
        rownames = FALSE,
        colnames = "none",
        options = list(searching = FALSE)
      )
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   output$sorttable <- DT::renderDataTable(dt)
                 })
    shinyjs::hide('actionsortdatatable')
    shinyjs::show('sorttable')
  })
  
  # sort tool gene list $use ----
  observeEvent(input$sorttable_rows_all,
               ignoreInit = TRUE,
               ignoreNULL = TRUE,
               {
                 oldname <-
                   last(grep("^Filter ", names(LIST_DATA$gene_file)))
                 oldname1 <- names(LIST_DATA$gene_file)[oldname] %>% str_split_fixed(., " n = ",n=2)
                 newname <-
                   paste(oldname1[1], "n =",
                         length(input$sorttable_rows_all))
                 if (newname != names(LIST_DATA$gene_file)[oldname]) {
                   print("sort filter $use")
                   names(LIST_DATA$gene_file)[oldname] <<- newname
                   LIST_DATA$gene_file[[newname]]$use <<-
                     tibble(gene = LIST_DATA$gene_file[[newname]]$full$gene[input$sorttable_rows_all])
                   ol <- input$sortGeneList
                   if (!ol %in% names(LIST_DATA$gene_file)) {
                     ol <- newname
                     reactive_values$Picker_controler <- 
                       c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
                     reactive_values$pickerfile_controler <-
                       input$sortSamples
                   } else {
                     reactive_values$pickerfile_controler <- ""
                     reactive_values$Picker_controler <- 
                       c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
                   }
                   updateSelectInput(
                     session,
                     "sortGeneList",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   output$valueboxsort <- renderValueBox({
                     valueBox(
                       n_distinct(LIST_DATA$gene_file[[newname]]$use),
                       "Gene List Filter",
                       icon = icon("list"),
                       color = "green"
                     )
                   })
                 }
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
                           }")
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
      menuItem("CDF Tools", tabName = "cdftool", icon = icon("ruler-combined"))
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
                      accept = c('.table'),
                      multiple = TRUE
                    ),
                    helpText("load windowed bedGraph file(s)"),
                    br(),
                    DT::dataTableOutput('loadedfilestable'),
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
                             selectInput("selectsave", "", choices = "Complete", ),
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
                    accept = c('.txt', 'bed')
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
            boxDropdownItem("Update sample color", id = "dropcolor", icon = icon("palette")),
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
                awesomeCheckbox("checkboxlog2", label = "log2")),
              column(
                3,
                numericInput("numericYRangeLow", label = "Plot Y min:", value = 0)
              ),
              column(
                3,
                numericInput("numericYRangeHigh", label = "Plot Y max:", value = 0)
              ),
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
          # hidden(
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
          # )
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
                status = "primary",
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
                status = "primary",
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
                status = "primary",
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
                title = "Clusters/Groups",
                width = 6,
                status = "primary",
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
          status = "purple",
          solidHeader = TRUE,
          title = "QC Options",
          fileInput(
            "filetable",
            label = "",
            accept = c('.table'),
            multiple = TRUE
          )
        )
      ),
      tabItem(
        # filenorm ----
        tabName = "filenorm",
        box(
          status = "purple",
          solidHeader = TRUE,
          title = "QC Options",
          fileInput(
            "filetable",
            label = "",
            accept = c('.table'),
            multiple = TRUE
          )
        )
      ),
      tabItem(
        # genelists ----
        tabName = "genelists",
        box(
          status = "purple",
          solidHeader = TRUE,
          title = "QC Options",
          fileInput(
            "filetable",
            label = "",
            accept = c('.table'),
            multiple = TRUE
          )
        )
      ),
      tabItem(
        # filter/sort tools ----
        tabName = "sorttool",
        div(
          id = "enablemainsort",
          box(
            title = "Filter by sum rank",
            solidHeader = T,
            width = 6,
            status = "navy",
            fluidRow(column(
              12,
              sliderInput(
                "slidersortpercent",
                label = "% select:",
                post = "%",
                min = 1,
                max = 100,
                value = 75
              )
            ),
            ),
            fluidRow(column(
              6,
              pickerInput(
                "selectsorttop",
                "Filter Options",
                choices = c("Top%", "Middle%", "Bottom%"),
                selected = "Middle%"
              )
            )),
            fluidRow(align="center",
                     actionButton("actionsorttool", "filter sum")
            )
          ),
          box(
            title = "filter by percentile distribution",
            solidHeader = T,
            width = 6,
            status = "navy",
            fluidRow(column(
              6,
              numericInputIcon("numericsortmin",
                               "min", 
                               value = "1",
                               max = "100", min="0",
                               step = "1",
                               icon = icon("percent")
              )
              ),
              column(
                6,
                numericInputIcon("numericsortmax",
                                 "max", 
                                 value = "99.5",
                                 max = "100", min="0",
                                 step = "1",
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
                  min = 0,
                  max = 80,
                  value = c(0, 80)
                )
              )
              ),
              column(width = 6,
              pickerInput("sortSamples", label = "select sample(s)",
                          choices = "select sample(s)",selected = "select sample(s)",
                          multiple = TRUE),
              div(
                style = "margin-bottom: -20px;",
                awesomeCheckbox("checkboxfilterall","Filter on all if any", value = TRUE)
              )
              )
              ),
          div(
            id = "hidesorttable",
            box(
              title = "Filter Table",
              solidHeader = T,
              width = 12,
              status = "navy",
              actionButton("actionsortdatatable", "Show gene list"),
              helpText("All filtering applied to gene list usage elsewhere"),
              DT::dataTableOutput('sorttable')
            )
          ),
          valueBoxOutput("valueboxsort")
        )
      ),
      tabItem(
        # raito tool ----
        tabName = "ratiotool",
        box(
          status = "purple",
          solidHeader = TRUE,
          title = "QC Options",
          fileInput(
            "filetable",
            label = "",
            accept = c('.table'),
            multiple = TRUE
          )
        )
      ),
      tabItem(
        # cluster tools ----
        tabName = "clustertool",
        box(
          status = "purple",
          solidHeader = TRUE,
          title = "QC Options",
          fileInput(
            "filetable",
            label = "",
            accept = c('.table'),
            multiple = TRUE
          )
        )
      ),
      tabItem(
        # cdf ----
        tabName = "cdftool",
        box(
          status = "purple",
          solidHeader = TRUE,
          title = "QC Options",
          fileInput(
            "filetable",
            label = "",
            accept = c('.table'),
            multiple = TRUE
          )
        )
      )
    )
  ),
  controlbar = dashboardControlbar(),
  title = "DashboardPage"
)

# exicute ----
shinyApp(ui = ui, server = server)