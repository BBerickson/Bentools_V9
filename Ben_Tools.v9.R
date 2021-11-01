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
      menuItem(
        "Load Data",
        tabName = "loaddata",
        icon = icon("file-import")
      ),
      menuItem(
        "QC/Options",
        tabName = "qcOptions",
        icon = icon("clipboard-check")
      ),
      menuItem("Plot", tabName = "mainplot", icon = icon("chart-area")),
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
        # mainplot ----
        tabName = "mainplot",
        fluidRow(box(
          width = 12, 
          status = "navy",
          solidHeader = TRUE,
          title = "", dropdownMenu = boxDropdown(colourInput("colourhex1", "Select color HEX"),
                                                 boxDropdownItem(colourInput("colourhex2", "Select color HEX"))),
          sidebar = boxSidebar(
            id = "sidebarcolor",
            width = 25,
            colourInput("colourhex", "Select color HEX")
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