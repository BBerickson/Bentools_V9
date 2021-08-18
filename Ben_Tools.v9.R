# Created by Benjamin Erickson BBErickson@gmail.com

# load/save data tab
#   
#   load other gene list
#   remove gene list (excluded Main list)
#   reset app
#   save gene lists, options for list, bed, tool
# data options and QC tab
# plot tab
# set observeEvents for all items, icon size/color?
# work on adding other boxes and sizes and functions

source("R_scripts/helpers.R",local = TRUE)

# run load needed packages using my_packages(x) ----
suppressPackageStartupMessages(my_packages(
  c("tidyverse",
    "shiny",
    "shinydashboard",
    "shinydashboardPlus",
    "shinyWidgets",
    "shinyjs",
    "RColorBrewer")
))

source("R_scripts/functions.R",local = TRUE)

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 500MB. ----
options(shiny.maxRequestSize = 500 * 1024 ^ 2)

LIST_DATA <<- list(
  table_file = NULL,
  # gene bin score set
  gene_file = NULL,
  # holds $Compleat genes from files and $gene file(s)
  gene_info = NULL,
  # for holding meta data gene file(s) [c("gene_list", "set", "color", plot?, legend)]
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
  addCssClass(selector = "a[data-value='mainplot']", class = "inactiveLink")
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
    }
    # enables tabs after loading file
    shinyjs::enable("startoff")
    shinyjs::reset("filetable")
    removeCssClass(selector = "a[data-value='qcOptions']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='mainplot']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='filenorm']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='genelists']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='sorttool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='ratiotool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='clustertool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='cdftool']", class = "inactiveLink")
  })
  
  
  # loads gene list file ----
  observeEvent(input$filegene1, {
    print("load gene file")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   # load info, update select boxes, switching works and changing info and plotting
                   LD <- LoadGeneFile(
                     input$filegene1$datapath,
                     input$filegene1$name,
                     LIST_DATA
                   )
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
  })
  
}

# UI -----
ui <- dashboardPage(
  skin = "purple-light",
  options = list(sidebarExpandOnHover = F),
  header = dashboardHeader(userOutput("user")),
  sidebar = dashboardSidebar(id = "sidebar", minified = TRUE, collapsed = TRUE,
                             tags$head(tags$style(".inactiveLink {
                            pointer-events: none;
                           cursor: default;
                           }")), # disables tabs on start
                             sidebarMenu(
                               id = "leftSideTabs",
                               menuItem("Load Data", tabName = "loaddata", icon = icon("file-import")),
                               menuItem("QC/Options", tabName = "qcOptions", icon = icon("clipboard-check")),
                               menuItem("Plot", tabName = "mainplot", icon = icon("area-chart")),
                               menuItem("Norm data", tabName = "filenorm", icon = icon("files-o")),
                               menuItem("Compare Lists", tabName = "genelists", icon = icon("gears")),
                               menuItem("Filter Tool", tabName = "sorttool", icon = icon("gears")),
                               menuItem("Ratio Tool", tabName = "ratiotool", icon = icon("gears")),
                               menuItem("Cluster Tools", tabName = "clustertool", icon = icon("gears")),
                               menuItem("CDF Tools", tabName = "cdftool", icon = icon("gears"))
                             )
  ),
  body = dashboardBody(
    useShinyjs(),
    tabItems(
      # load data tab ----
      tabItem(tabName = "loaddata",
              box(status = "navy",
                  solidHeader = TRUE,
                  title = "Load .table/URL.txt file",
                  width = 6,
                  style = "height: 150px;",
                  align="center",
                  fileInput(
                    "filetable", width = "75%",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ),
                  helpText("load windowed bedGraph file(s)")
              ),
              hidden(div(
                id = "startoff",
              box(
                title = "Load Gene list",
                width = 6,
                style = "height: 150px;" ,
                solidHeader = TRUE,
                status = "navy",
                align="center",
                fileInput("filegene1", width = "75%",
                          label = "",
                          accept = c('.txt')),
                helpText("load gene list")
              )))
      ),
      tabItem(tabName = "qcOptions",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  sidebar = boxSidebar(
                    id = "sidebarGenelistColor",
                    icon = icon("palette"),
                    width = 50,
                    pickerInput(inputId = "kbrewer",
                                label = "color brewer theme",
                                choices = c("select", kBrewerList),
                                selected = "select"
                    ),
                    actionButton("BttnNewColor", "Set color same as Compleat")
                  ),
                  fluidRow(align="center",
                           
                           pickerInput("selectdataoption", "", 
                                       width = 300, choices = "Load data file"),
                           pickerInput("selectgenelistoptions", "",
                                       width = 300, choices = "Compleat"),
                           actionButton("actionremovegene", "Remove Gene list"))
              )
      ),
      tabItem(tabName = "mainplot",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "filenorm",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "genelists",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "sorttool",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "ratiotool",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "clustertool",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "cdftool",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      )
    )
  ),
  controlbar = dashboardControlbar(),
  title = "DashboardPage"
)

# exicute ----
shinyApp(ui = ui, server = server)