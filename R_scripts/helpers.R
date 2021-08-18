
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

# color Brewer set that is active to use in plot remove all yellows ----
kListColorSet <- brewer.pal(11, kBrewerList[2]) %>% grep("#FF",.,value = T,invert = T)

# math options available ----
kMathOptions <- c("mean", "sum", "median")

# basic test for valid rgb color and type and
# switches hex <--> rgb, able to apply tint colors
RgbToHex <- function(x, 
                     convert = "hex", 
                     tint = FALSE) {
  if(str_count(x, ",") == 2){
    # hex <- rgb
    myhex <- try(rgb(str_split_fixed(x,",",n=3),
                     maxColorValue = 255), silent = TRUE)
    if("try-error" %in% class(myhex)){
      myhex <- "#000000"
      mygrb <- "0,0,0"
    } else {
      myrgb <- x
      if (is.numeric(tint) & between(tint,0,1)) {
        myrgb <- as.numeric(str_split_fixed(myrgb,",",n=3))
        myrgb <-
          paste(round(myrgb + (255 - myrgb) * tint), collapse = ",")
        myhex <- rgb(str_split_fixed(myrgb,",",n=3),
                     maxColorValue = 255)
      } 
    }
  } else {
    # rgb <- hex (tint)
    myrgb <- try(col2rgb(x), silent = TRUE)
    if("try-error" %in% class(myrgb)){
      myhex <- "#000000"
      myrgb <- "0,0,0"
    } else {
      if (is.numeric(tint) & between(tint,0,1)) {
        myrgb <-
          paste(round(as.numeric(myrgb) + (255 - as.numeric(myrgb)) * tint), collapse = ",")
        myhex <- rgb(str_split_fixed(myrgb,",",n=3),
                     maxColorValue = 255)
        
      } else {
        myhex <- x
        myrgb <- paste(myrgb, collapse = ",")
      }
    }
  }
  if(convert == "hex"){
    return(myhex)
  } else {
    return(myrgb) 
  }
}

# reads in table file(s), tests, fills out info and returns list_data
LoadTableFile <-
  function(file_path,
           file_name,
           list_data) {
    my_color <- NULL
    my_landl <- NA
    legend_nickname <- NULL
    my_remote_file <- NULL
    file_count <- length(list_data$table_file)
    # shiny progress bar
    setProgress(1, detail = "start gathering information on file(s)")
    # tests if loading a file with a list of address to remote files, requires .url.txt in file name
    for (i in seq_along(file_name)) {
      if (str_detect(file_name[i], ".url.txt")) {
        num_col <-
          try(count_fields(file_path[i],
                           n_max = 1,
                           skip = 1,
                           tokenizer = tokenizer_delim(" ")),silent = T)
        if ("try-error" %in% class(num_col)) {
          showModal(modalDialog(
            title = "Information message",
            paste(file_name[i], "cant find file to load"),
            size = "s",
            easyClose = TRUE
          ))
          next()
        }
        if(num_col > 1){
          col_names <- c("file","type","nick","color")
          col_types <- cols(file=col_character(),
                            type=col_character(),
                            nick=col_character(),
                            color=col_character())
        } else {
          col_names <- c("file")
          col_types <- cols(file=col_character())
        }
        meta_data <- suppressMessages(read_delim(file_path[i],delim = " ", 
                                                 col_names = col_names,
                                                 col_types = col_types))
        if(num_col == 1){
          meta_data <-  meta_data %>% mutate(type = NA,
                                             nick = NA,
                                             color = NA)
        }
        my_remote_file <- c(my_remote_file, meta_data$file)
        my_color <- c(my_color, meta_data$color) 
        legend_nickname <- c(legend_nickname, str_replace(meta_data$nick,"\\.","_")) 
        my_landl <- last(meta_data$type)
        if(i > 1){
          file_path[i] <- NULL
        } else {
          file_path <- NULL
        }
      } else {
        my_color <- c(my_color, NA)
        legend_nickname <- c(legend_nickname, last(str_split(file_name[i],"/",simplify = T)) %>% 
                               str_remove(., ".table")) %>% str_replace("\\.","_")
      }
    }
    if (!is.null(my_remote_file)) {
      file_path <- c(file_path, my_remote_file)
    }
    # shiny progress bar
    setProgress(2, detail = "getting meta data")
    # loop each item in file_path
    for (x in seq_along(file_path)) {
      # gets number of columns in file, used to guess how to deal with file
      #  and checks if file exits
      num_col <-
        try(count_fields(file_path[x],
                         n_max = 1,
                         skip = 1,
                         tokenizer = tokenizer_tsv()),silent = T)
      if ("try-error" %in% class(num_col)) {
        showModal(modalDialog(
          title = "Information message",
          paste(file_name[x], "cant find file to load"),
          size = "s",
          easyClose = TRUE
        ))
        next()
      }
      # checking/creating meta data
      if(is.na(legend_nickname)[x]){
        legend_nickname[x] <- last(str_split(file_path[x],"/",simplify = T)) %>% 
          str_remove(., ".table") %>% str_replace("\\.","_")
      }
      if(file_count > 0 ){
        # checks if file with same name is in master list of lists
        if (legend_nickname[x] %in% list_data$gene_info$set) {
          showModal(modalDialog(
            title = "Information message",
            paste(legend_nickname[x], "has already been loaded"),
            size = "s",
            easyClose = TRUE
          ))
          next()
        }
      }
      if(is.na(my_landl)){
        my_landl <- str_extract(legend_nickname[x], "^(\\d)+") %>% 
          replace_na("none") 
      }
      
      if (is.na(my_color[x])) {
        my_color[x] <- sample(suppressWarnings(brewer.pal(11, sample(kBrewerList, size=1))) %>% 
                                grep("#FF",.,value = T,invert = T),size = 1)
      } else {
        my_color[x] <- RgbToHex(my_color[x], convert = "hex")
      }
      
      if (num_col > 6) {
        # guesses is in wide format
        # shiny progress bar
        setProgress(3, detail = "loading wide file and converting")
        tablefile <- suppressMessages(
          read_tsv(
            file_path[x],
            comment = "#",
            col_names = c("gene", 1:(num_col - 1)),
            skip = 1
          ) %>%
            gather(., bin, score, 2:(num_col))
        ) %>%
          dplyr::select(gene, bin, score) %>%
          dplyr::mutate(set = legend_nickname[x],
                        bin = as.numeric(bin),
                        score = as.numeric(score)) %>%
          na_if(Inf) %>%
          replace_na(list(score = 0)) %>% 
          distinct(gene,bin,.keep_all = T)
        # guesses is in long bedtools from (bed or bedGraph)
      } else {
        if (num_col == 4) {
          # settings for new style with meta data info
          col_names <- c("gene", "bin", "score", "set")
          # settings for reading in bedGraph file style
        } else if (num_col == 3) {
          col_names <- c("gene", "bin", "score")
        } else {
          showModal(
            modalDialog(
              title = "I dont know how to load this file",
              "I use binned coverage files: gene bin score optinal(set) ",
              size = "s",
              easyClose = TRUE
            )
          )
          next()
        }
        # shiny progress bar
        setProgress(3, detail = "downloading/reading in file")
        # reads in file
        tablefile <-
          suppressMessages(read_tsv(file_path[x],
                                    comment = "#",
                                    col_names = col_names)) %>%
          dplyr::mutate(set = legend_nickname[x]) %>% na_if(Inf) %>%
          replace_na(list(score = 0)) %>% distinct(gene,bin,.keep_all = T)
      }
      num_bins <- n_distinct(tablefile$bin)
      # shiny progress bar
      setProgress(4, detail = "Checking form problems")
      # checks the number of bins and gene naming skeam is consistent
      if (file_count > 0) {
        if (num_bins != list_data$x_plot_range[2]) {
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
        # test data is compatible with already loaded data
        gene_names <-
          semi_join(tablefile, list_data$gene_file[[1]]$use, by = "gene") %>% distinct(gene)
        if (n_distinct(gene_names$gene) == 0) {
          showModal(
            modalDialog(
              title = "Information message",
              " No genes in common ",
              size = "s",
              easyClose = TRUE
            )
          )
          break()
        } else {
          # make complete gene list
          gene_names <-
            full_join(tablefile, list_data$gene_file[[1]]$use, by = "gene") %>% 
            distinct(gene)
        }
      } else {
        gene_names <- distinct(tablefile, gene)
        list_data$x_plot_range <- c(1, num_bins)
        list_data$STATE[3] <- my_landl
      }
      # sets master gene list name
      my_name <- paste("Compleat\nn =", n_distinct(gene_names$gene))
      if (file_count > 0) {
        list_data$gene_info <- list_data$gene_info %>% 
          dplyr::mutate(gene_list = if_else(str_detect(gene_list,"^Compleat\nn =") , my_name, gene_list))
        names(list_data$gene_file)[1] <- my_name
      }
      # shiny progress bar
      setProgress(5, detail = "building data and adding to tools")
      if (list_data$STATE[2] == 0 &
          n_distinct(list_data$table_file$set) < 3) {
        oo <- legend_nickname[x]
      } else {
        oo <- "0"
      }
      # saves data in list of lists
      list_data$table_file <- distinct(bind_rows(list_data$table_file, tablefile))
      list_data$gene_file[[my_name]]$use <- gene_names
      list_data$gene_file[[my_name]]$info <- paste("full gene list",
                                                   Sys.Date())
      list_data$gene_info <- distinct(bind_rows(list_data$gene_info,tibble(
        gene_list = my_name,
        set = legend_nickname[x],
        mycol = my_color[x],
        onoff = oo,
        sub = " ",
        plot_set = " "
      )))
      
      file_count <- 1
      setProgress(6, detail = "done with this file")
    }
    return(list_data)
  }

# reads in 1 column gene list file(s), tests, fills out info and returns list_data
LoadGeneFile <-
  function(file_path,
           file_name,
           list_data) {
    my_remote_file <- NULL
    legend_nickname <- NULL
    if(length(list_data$table_file) < 1){
      print("needed  or does shiny handle? add message?")
      return(list_data)
    }
    # shiny progress bar
    setProgress(1, detail = "start gathering information on file(s)")
    # tests if loading a file with a list of address to remote files, requirs .url.txt in file name
    for (i in seq_along(file_name)) {
      if (str_detect(file_name[i], ".url.txt")) {
        meta_data <- suppressMessages(read_delim(file_path[i],delim = " ",col_names = FALSE))
        my_remote_file <- c(my_remote_file, meta_data$file)
        if(i > 1){
          file_path[i] <- NULL
        } else {
          file_path <- NULL
        }
      }
    }
    if (!is.null(my_remote_file)) {
      file_path <- c(file_path, my_remote_file)
    }
    # shiny progress bar
    setProgress(2, detail = "downloading/reading in file")
    # loop thourgh each item in file_path
    for (x in seq_along(file_path)) {
      # createing nickname
      legend_nickname[x] <- last(str_split(file_name[x],"/",simplify = T)) %>% 
        str_remove(., ".txt") %>% str_replace("\\.","_")
      # checks if file with same nickname has already been loaded
      if (legend_nickname[x] %in% names(list_data$gene_file)) {
        showModal(modalDialog(
          title = "Information message",
          paste(file_name[x], "has already been loaded"),
          size = "s",
          easyClose = TRUE
        ))
        next()
      }
      # gets number of columns in file, used to guess how to deal with file
      #  and checks if file exits
      num_col <-
        try(count_fields(file_path[x],
                         n_max = 1,
                         skip = 1,
                         tokenizer = tokenizer_tsv()),silent = T)
      if ("try-error" %in% class(num_col)) {
        showModal(modalDialog(
          title = "Information message",
          paste(file_name[x], "cant find file to load"),
          size = "s",
          easyClose = TRUE
        ))
        next
      }
      if(num_col == 1){
        #normal gene list
        col_names <- c("gene")
      } else {
        showModal(
          modalDialog(
            title = "I dont know how to load this file",
            "I expected a 1 colunm file",
            size = "s",
            easyClose = TRUE
          )
        )
        next()
      }
      # reads in file
      tablefile <-
        suppressMessages(read_tsv(file_path[x],
                                  comment = "#",
                                  col_names = col_names)) %>%
        distinct(gene)
      # shiny progress bar
      setProgress(3, detail = "Checking form problems")
      # checks gene list is a subset of what has been loaded
      gene_names <-
        semi_join(tablefile, list_data$gene_file[[1]]$use, by = "gene") %>% distinct(gene)
      # test data is compatible with already loaded data
      if (n_distinct(gene_names$gene) == 0) {
        showModal(
          modalDialog(
            title = "Information message",
            " couldn't find a exact match, checking for partial match ",
            size = "s",
            easyClose = TRUE
          )
        )
        # tries to grep lists and find matches
        # shiny progress bar
        setProgress(4, detail = "looking for gene name matches")
        tablefile <-
          distinct(tibble(gene = str_subset(
            list_data$gene_file[[1]]$use$gene, tablefile$gene
          )))
        if (n_distinct(tablefile$gene) == 0) {
          showModal(
            modalDialog(
              title = "Information message",
              " No genes found after pattern matching search",
              size = "s",
              easyClose = TRUE
            )
          )
          return()
        }
        showModal(
          modalDialog(
            title = "Information message",
            " Don't forget to save the gene list for future use",
            size = "s",
            easyClose = TRUE
          )
        )
      }
      # adds full n count to nickname
      my_name <- paste0(legend_nickname[x], "\nn = ", n_distinct(tablefile$gene))
      # preps meta data
      gene_info <- list_data$gene_info %>% 
        dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
        dplyr::mutate(gene_list = my_name, sub = " ", onoff = "0",
                      plot_set = " ")
      # shiny progress bar
      setProgress(5, detail = "building data and adding list")
      # saves data in list of lists
      list_data$gene_file[[my_name]]$use <- distinct(tablefile, gene)
      list_data$gene_file[[my_name]]$full <- tablefile
      list_data$gene_file[[my_name]]$info <-
        paste("Loaded gene list from file",
              legend_nickname[x],
              Sys.Date())
      list_data$gene_info <- bind_rows(list_data$gene_info, gene_info)
      
      setProgress(6, detail = "done with this file")
    }
    return(list_data)
  }
