# Created by Benjamin Erickson BBErickson@gmail.com

# load/save/remove data tab
#
# data options and QC tab
#     stats on files, distribution of mean signals (5,4,3), number of 0,s na's, peaks, ???
#     ability to find and remove outliers and filter (add to server R RNAseq filtering outliers )
# plot tab , start with empty UI, apply math fun, plot fun, lines and labels funs,
#   plot box, file and gene list selections, plot math and axis options, color and font size, ordering, lines and labels, ttest,
# set observeEvents for all items, icon size/color?
# work on adding other boxes and sizes and functions

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
    "DT"
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
  
  output$user <- renderUser({
    dashboardUser(
      name = "BenTools V9.a",
      image = "ben head.jpg",
      title = "Benjamin Erickson",
      subtitle = "BBErickson@gmail.com"
    )
  })
  # disables tabs on start
  addCssClass(selector = "a[data-value='qcOptions']", class = "inactiveLink")
  # addCssClass(selector = "a[data-value='mainplot']", class = "inactiveLink")
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
    updatePickerInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_file),
      selected = names(LIST_DATA$gene_file)[1]
    )
    ff <- distinct(LIST_DATA$table_file, set)$set
    updatePickerInput(session,
                      "selectdataoption",
                      choices = ff)
    # first time starting
    if (LIST_DATA$STATE[1] == 0) {
      LIST_DATA$STATE[1] <<- 1
      shinyjs::show("startoff")
      shinyjs::show("startoff2")
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
    updatePickerInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_file),
      selected = last(names(LIST_DATA$gene_file))
    )
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
  
  # plots when action button is pressed ----
  observeEvent(input$actionmyplot, ignoreInit = TRUE, {
    print("plot button")
    
    shinyjs::hide("actionmyplotshow")
  })
  
  # updates plot ----
  observeEvent(input$actionButtonUpDatePlot, ignoreInit = T, {
    print("actionButtonUpDatePlot")
    updateBoxSidebar(id = "sidebarmath")
  })
  
  # opens color select diolog box
  observeEvent(input$dropcolor, ignoreInit = T, {
    showModal(modalDialog(
      title = "Information message",
      " Don't forget to save the gene list for future use",
      size = "l",
      easyClose = F,
      footer = tagList(
        box(
          width = 12,
          status = "primary",
          solidHeader = T,
          collapsible = FALSE,
          collapsed = FALSE,
          div(
            style = "padding-left: 25px; display:inline-block;",
            selectInput(
              selectize = T,
              "selectlineslabels",
              width = "200px",
              label = "quick set lines and labels",
              choices = c("Choose one" = "",
                          kLinesandlabels)
            )
          ),
          column(12,
                 div(
                   style = "padding:2px; display:inline-block;",
                   numericInput(
                     "numerictss",
                     "TSS bin",
                     value = 15,
                     min = 0,
                     max = 100
                   )
                 ),
                 div(
                   style = "padding:2px; display:inline-block;",
                   textInput("numerictssname", value = "TSS", label = "lable",width = "50px"),
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
                   style = "padding:2px; display:inline-block;",
                   numericInput(
                     "numericbody2",
                     "4|3 bin",
                     value = 40,
                     min = 0,
                     max = 100
                   )
                 ),
                 div(
                   style = "padding:2px; display:inline-block;",
                   numericInput(
                     "numerictes",
                     "pA bin",
                     value = 45,
                     min = 0,
                     max = 100
                   )
                 ),
                 div(
                   style = "padding:2px; display:inline-block;",
                   textInput("numerictesname", value = "pA", label = "lable",width = "50px"),
                 ),
                 div(
                   style = "padding:2px; display:inline-block;",
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
                   style = "padding:2px 8px 2px 2px; display:inline-block;",
                   numericInput(
                     "numericlabelspaceing",
                     "every bin",
                     value = 5,
                     min = 0,
                     max = 100
                   )
                 ),
                 actionButton("actionlineslabels", "UPDATE PLOT")
          ),
          helpText("For 543 style 0 > TSS < 5|4 < 4|3 < pA < max bin"),
          div(
            textInput("landlnames", "", label = "Yaxis labels"),
            textInput("landlposition", "", label = "Yaxis lable position (numbers only)")
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
                  label = 'TSS line and lable color',
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
                  label = '5|4 line and lable color',
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
                  label = '4|3 line and lable color',
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
                  label = 'TES line and lable color',
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
        modalButton("Cancel"),
        actionButton("ok", "OK")
      )
    ))
  })
  
  observeEvent(input$ok, {
    colorModal()
    removeModal()
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
      menuItem("Compare Lists", tabName = "genelists", icon = icon("cogs")),
      menuItem("Filter Tool", tabName = "sorttool", icon = icon("cogs")),
      menuItem("Ratio Tool", tabName = "ratiotool", icon = icon("cogs")),
      menuItem("Cluster Tools", tabName = "clustertool", icon = icon("cogs")),
      menuItem("CDF Tools", tabName = "cdftool", icon = icon("cogs"))
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
          sidebar = boxSidebar(
            id = "sidebarmath",
            width = 50,
            icon = icon("cogs"),
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
              inputId = "actionButtonUpDatePlot",
              label = "Update Plot",
              style = "unite",
              color = "default",
              size = "sm"
            )
            )
            )
          ),
          shinycssloaders::withSpinner(plotOutput("plot"), type = 4),
          hidden(
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
          )
        )),
        fluidRow(box(
          width = 12, 
          status = "navy",
          solidHeader = TRUE,
          title = "", 
          dropdownMenu = boxDropdown(
            boxDropdownItem("Update sample color", id = "dropcolor", icon = icon("palette")),
            dropdownDivider(),
            boxDropdownItem("Lines and Labels", id = "droplinesandlabels", icon = icon("chart-bar")),
            dropdownDivider(),
            boxDropdownItem("t-Test", id = "dropttest", icon = icon("chart-line")),
            dropdownDivider(),
            boxDropdownItem("Font and line size", id = "dropfontsize", icon = icon("font"))
          ),
            box(title = "Main",
                           width = 6,
                           status = "primary",
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
          sidebar = boxSidebar(
            id = "sidebarGenelistColor",
            icon = icon("palette"),
            width = 50,
            pickerInput(
              inputId = "kbrewer",
              label = "color brewer theme",
              choices = c("select", kBrewerList),
              selected = "select"
            ),
            actionButton("BttnNewColor", "Set color same as Complete")
          ),
          fluidRow(
            align = "center",
            
            pickerInput("selectdataoption", "",
                        width = 300, choices = "Load data file"),
            pickerInput(
              "selectgenelistoptions",
              "",
              width = 300,
              choices = "Complete"
            ),
            actionButton("actionremovegene", "Remove Gene list")
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
        # sort tools ----
        tabName = "sorttool",
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