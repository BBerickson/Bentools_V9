
# remove search bar
# Filter gene list show data table ----
observeEvent(input$actionsortdatatable, ignoreInit = TRUE, {
  print("show data table")
  if (any(grep("^Filter", names(LIST_DATA$gene_file)) > 0)) {
    newnames <-
      gsub("(.{20})", "\\1... ", names(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$full))
    dt <- datatable(
      LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$full,
      colnames = c("RANK",strtrim(newnames, 30)),
      class = 'cell-border stripe compact',
      filter = 'top',
      caption = LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$info,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        scrollY = TRUE,
        autoWidth = F,
        width = 5,
        columnDefs = list(
          list(className = 'dt-center ', targets = "_all"),
          list(targets = 0, width = 2),
          list(
            targets = 1,
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


# plot testing
library(tidyverse)

tt <- LIST_DATA$table_file %>% group_by(set,bin) %>% summarise(score=mean(score),.groups = "drop") %>% 
  mutate(bin=as.numeric(bin))
ttt <- LIST_DATA$gene_info %>% filter(gene_list == "Complete") %>% select(set,mycol) %>% inner_join(tt,.)



p2 <- ggline(ttt,x="bin",y="score",plot_type = "l",color = "set",palette = "mycol",size = 3)
p2 + scale_x_discrete(breaks = seq(1,by = 10,length.out = (80 / 10) + 1))
ggpar(p2,
      legend = "right", legend.title = "Complete",
      font.legend = c(12)) + 
  scale_x_discrete(breaks = seq(1,by = 10,length.out = (79 / 10)+1),
                   labels = seq(1,by = 10,length.out = (79 / 10)+1))





library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
shinydashboardPlusGallery()

shinyApp(
  ui = dashboardPage(
    header = dashboardHeader(),
    body = dashboardBody(
      box(
        title = "Update box sidebar", 
        closable = TRUE, 
        width = 12,
        height = "500px",
        solidHeader = FALSE, 
        collapsible = TRUE,
        actionButton("update", "Toggle card sidebar"),
        sidebar = boxSidebar(
          id = "mycardsidebar",
          width = 25,
          sliderInput(
            "obs", 
            "Number of observations:",
            min = 0, 
            max = 1000, 
            value = 500
          )
        ),
        plotOutput("distPlot")
      )
    ),
    sidebar = dashboardSidebar()
  ),
  server = function(input, output, session) {
    observe(print(input$mycardsidebar))
    
    output$distPlot <- renderPlot({
      hist(rnorm(input$obs))
    })
    
    observeEvent(input$update, {
      updateBoxSidebar("mycardsidebar")
    })
    
  }
)