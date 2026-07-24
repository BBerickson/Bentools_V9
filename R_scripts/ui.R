# Shiny UI. Sourced by app.R; defines the object 'ui'.

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
                    accept = c('text/plain', 'application/gzip', 'application/x-gzip', '.txt', '.tsv', '.bed.gz', '.bed'),
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
                  checkboxInput("mygroup",
                                label = "Plot Group:",
                                value = FALSE
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
                helpText("!! Work in progress !!"),
                fluidRow(
                  column(
                    4,
                    div(style = "margin-top: 20px;",
                    awesomeCheckbox("checkboxlog2Violin", 
                                    label = "Log2 transform",
                                    value = TRUE)
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
                                 value = 9)
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

