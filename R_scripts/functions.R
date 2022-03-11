
# color Brewer set that is active to use in plot remove all yellows ----
kListColorSet <- brewer.pal(11, kBrewerList[2]) %>% grep("#FF",.,value = T,invert = T)

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

# finds first partial match to gene list input 
MatchGenes <- function(common_list, gene_list){
  if(str_detect(tablefile$gene[1],"_")){
    tablefile <-
      map(paste0(";", gene_list$gene,"\\|"), str_subset, string = common_list$gene) 
  } else{
    tablefile <-
      map(paste0("\\|", gene_list$gene,"$"), str_subset, string = common_list$gene) 
  }
  tablefile <- distinct(tibble(gene = unlist(tablefile)))
  return(tablefile)
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
      if ("try-error" %in% class(num_col) & num_col > 0) {
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
      
      # matrix file
      if(length(grep(".matrix.gz", file_name[x])) == 1){
        num_bins <-
          count_fields(file_path[x],
                       n_max = 1,
                       skip = 1,
                       tokenizer = tokenizer_tsv())
        tablefile <- suppressMessages(
          read_tsv(
            file_path[x],
            comment = "#",
            col_names = c("chr", "start", "end","gene", "value", "sign", 1:(num_bins - 6)),
            skip = 1)) %>%
          dplyr::select(-chr, -start, -end, -sign, -value) %>% 
          gather(., bin, score, 2:(num_bins-5)) %>%
          dplyr::mutate(bin = as.numeric(bin),
                        score = as.numeric(score),
                        set = legend_nickname[x]) %>%
          na_if(Inf) %>%
          replace_na(list(score = 0)) %>% 
          distinct(gene,bin,.keep_all = T)
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
            "I use binned coverage files with the following columns:\n
              gene bin score optinal(set) ",
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
      # checks the number of bins and gene naming is consistent
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
      if (file_count > 0) {
        list_data$gene_info <- list_data$gene_info %>% 
          dplyr::mutate(count = if_else(gene_list == "Complete" , paste("n =", n_distinct(gene_names$gene)), count))
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
      list_data$gene_file[["Complete"]]$use <- gene_names
      list_data$gene_file[["Complete"]]$info <- tibble(loaded_info = paste("full gene list",
                                                                           Sys.Date()))
      list_data$gene_info <- distinct(bind_rows(list_data$gene_info,tibble(
        gene_list = "Complete",
        count = paste("n =", n_distinct(gene_names$gene)),
        set = legend_nickname[x],
        mycol = my_color[x],
        onoff = oo,
        sub = " ",
        plot_set = " "
      )))
      # adding table file after gene list had been loaded
      for(gg in seq_along(list_data$gene_file)[-1]){
        if(list_data$gene_file[[gg]]$info$matching){
          new_gene_match <- MatchGenes(gene_names, list_data$gene_file[[gg]]$full %>% select(gene))
          if (n_distinct(new_gene_match$gene) != 0) {
            # fix name, fix info
            listname <- names(list_data$gene_file)[gg]
            list_data$gene_info <- list_data$gene_info %>%
              dplyr::mutate(count=if_else(gene_list == listname,paste("n =",n_distinct(new_gene_match$gene)), 
                                          count))
          }
        }
        # add missing data
        list_data$gene_info <- list_data$gene_info %>%
          dplyr::filter(gene_list == names(list_data$gene_file)[gg]) %>%
          dplyr::mutate(set = legend_nickname[x],
                        count = count[x],
                        mycol = my_color[x],
                        onoff = "0",
                        sub = " ",
                        plot_set = " ") %>%
          bind_rows(list_data$gene_info, .)
      }
      file_count <- 1
      setProgress(6, detail = "done with this file")
    }
    return(list_data)
  }

# reads in gene list file(s), tests, fills out info and returns list_data
LoadGeneFile <-
  function(file_path,
           file_name,
           list_data) {
    my_remote_file <- NULL
    legend_nickname <- NULL
    matching <- FALSE
    # shiny progress bar
    setProgress(1, detail = "start gathering information on file(s)")
    # tests if loading a file with a list of address to remote files, requires .url.txt in file name
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
    # loop though each item in file_path
    for (x in seq_along(file_path)) {
      # creating nickname
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
      if ("try-error" %in% class(num_col) & num_col > 0) {
        showModal(modalDialog(
          title = "Information message",
          paste(file_name[x], "cant find file to load"),
          size = "s",
          easyClose = TRUE
        ))
        next
      }
      if(num_col == 1 | str_detect(file_path[x], ".table")){
        #normal gene list
        col_names <- c("gene")
      } else if(str_detect(file_path[x], ".bed")){
        col_names <- 1:num_col
        col_names[4] <- "gene"
      } else {
        col_names <- TRUE
      }
      # reads in file
      tablefile <-
        suppressMessages(read_tsv(file_path[x],
                                  comment = "#",
                                  col_names = col_names)) 
      if(!"gene" %in% names(tablefile)){
        names(tablefile)[1] <- "gene"
      }
      tablefile <- tablefile %>% 
        distinct(gene,.keep_all = T) %>% mutate(gene=as.character(gene))
      # shiny progress bar
      setProgress(3, detail = "Checking form problems")
      # checks gene list is a subset of what has been loaded
      gene_names <- tablefile %>% select(gene) %>% 
        semi_join(., list_data$gene_file[["Complete"]]$use, by = "gene") %>% distinct(gene)
      # test data is compatible with already loaded data
      if (n_distinct(gene_names$gene) == 0) {
        showModal(
          modalDialog(
            title = "Information message",
            " couldn't find a exact match, checking for partial match
            This might take a few minutes",
            size = "s",
            easyClose = T
          )
        )
        # tries to grep lists and find matches
        # shiny progress bar
        setProgress(4, detail = "looking for gene name matches")
        gene_names <- MatchGenes(list_data$gene_file[["Complete"]]$use, tablefile %>% select(gene))
        if (n_distinct(gene_names$gene) == 0) {
          showModal(
            modalDialog(
              title = "Information message",
              " No genes found after pattern matching search",
              size = "s",
              easyClose = TRUE
            )
          )
          return()
        } else {
          matching = TRUE
          showModal(
            modalDialog(
              title = "Information message",
              paste("Found", n_distinct(gene_names$gene), "genes after pattern matching search"),
              size = "s",
              easyClose = TRUE
            )
          )
        }
      }
      # adds full n count to nickname
      my_name <- legend_nickname[x]
      # preps meta data
      gene_info <- list_data$gene_info %>% 
        dplyr::filter(gene_list == "Complete") %>% 
        dplyr::mutate(gene_list = my_name, 
                      count = paste("n =", n_distinct(gene_names$gene)),
                      sub = " ", onoff = "0",
                      plot_set = " ")
      # shiny progress bar
      setProgress(5, detail = "building data and adding list")
      # saves data in list of lists
      list_data$gene_file[[my_name]]$use <- distinct(gene_names, gene)
      list_data$gene_file[[my_name]]$full <- tablefile
      list_data$gene_file[[my_name]]$info <- tibble(loaded_info =
                                                      paste("Loaded gene list from file",
                                                            legend_nickname[x],
                                                            Sys.Date()),
                                                    matching = matching)
      list_data$gene_info <- bind_rows(list_data$gene_info, gene_info)
      
      setProgress(6, detail = "done with this file")
    }
    return(list_data)
  }

# color module dialog box update
colorModal <- function(list_data, tt){
  list_data$test <- tt
  return(list_data)
}

Active_list_data <-
  # input data list, output filtered and bound data with plot legend column
  function(list_data) {
    table_file <- list_data$table_file
    gene_file <- list_data$gene_file
    gene_info <- list_data$gene_info
    list_data_out <- NULL
    print("active data function")
    for ( i in names(gene_file) ){
      # checks to see if at least one file in list is active
      if (gene_info %>% dplyr::filter(gene_list == i & onoff != 0) %>% nrow() == 0) {
        next()
      } else {
        my_sel <- gene_info %>% dplyr::filter(gene_list == i & onoff != 0)
        tf <- table_file %>% 
          dplyr::filter(set %in% my_sel$onoff)
        gene_common <- tf %>% group_by(set) %>% distinct(gene) %>% ungroup()
        gene_common <- gene_common %>% 
          group_by(gene) %>% 
          filter(n_distinct(set)==n_distinct(gene_common$set)) %>% 
          distinct(gene) %>% ungroup()
        list_data_out[[i]] <-
          tf %>% 
          semi_join(., gene_common, by = "gene") %>%
          semi_join(., gene_file[[i]]$use, by = "gene") %>% 
          dplyr::mutate(., gene_list = i)
        # adds line brake at 20 character for legend spacing
        my_sel2 <- my_sel %>% dplyr::mutate(.,plot_set = paste(
          gsub("(.{20})", "\\1\n", set),
          gsub("(.{20})", "\\1\n", 
               str_split_fixed(i, "\nn = ", n=2)[,1]),
          paste0("n = ", n_distinct(list_data_out[[i]]$gene)),
          sep = '\n'
        )) %>% select(set,plot_set)
        list_data_out[[i]] <- list_data_out[[i]] %>% inner_join(.,my_sel2,by="set")
      }
    }
    return(bind_rows(list_data_out))
  }

ApplyMath <-
  function(list_data,
           use_math,
           relative_frequency,
           normbin) {
    setProgress(1, detail = paste("Gathering info"))
    # apply's math to data file
    if (relative_frequency == "rel gene frequency") {
      list_data <- list_data %>% group_by(plot_set, gene) %>%
        dplyr::mutate(score = score / abs(sum(score, na.rm = TRUE))) %>%
        ungroup()
    }
    list_data <- list_data %>% group_by(plot_set, bin) %>%
      summarise(value = get(use_math)(score, na.rm = T), .groups="drop")
    
    if (normbin > 0) {
      list_data <- list_data %>% 
        group_by(plot_set) %>%
        arrange(bin) %>%
        dplyr::mutate(value = value / abs(nth(value, normbin))) %>%
        ungroup()
    } else if (relative_frequency == "relative frequency") {
      list_data <- list_data %>%
        group_by(plot_set) %>%
        dplyr::mutate(value = value / abs(sum(value))) %>%
        ungroup()
    }
    return(list_data %>% dplyr::mutate(set=plot_set))
  }

# get min and max from apply math data set
YAxisValues <-
  function(apply_math,
           xBinRange,
           yBinRange = c(0, 100),
           log_2 = F) {
    tt <- group_by(apply_math, set) %>%
      dplyr::filter(bin %in% xBinRange[1]:xBinRange[2]) %>%
      ungroup() %>%
      summarise(min(value, na.rm = T), max(value, na.rm = T),.groups="drop") %>%
      unlist(., use.names = FALSE)
    tt <-
      c(tt[1] + (tt[1] * (yBinRange[1] / 100)), tt[2] + (tt[2] * ((yBinRange[2] -
                                                                     100) / 100)))
    if (log_2) {
      tt <- log2((abs(tt))^(sign(tt)))
    }
    tt
  }

# Sets y label fix
YAxisLabel <-
  function(use_math = "mean",
           relative_frequency = "none",
           norm_bin = 0,
           smoothed = F,
           log_2 = F) {
    use_y_label <- paste(use_math, "of bin counts")
    if (relative_frequency == "rel gene frequency") {
      use_y_label <- paste("RF per gene :", use_y_label)
    } else if (relative_frequency == "relative frequency") {
      use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                           "bins : RF")
    }
    if (norm_bin > 0) {
      if (relative_frequency == "rel gene frequency") {
        use_y_label <- paste(use_y_label, " : Norm bin ", norm_bin)
      } else {
        use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                             "bins : Normalize to bin ",
                             norm_bin)
      }
    }
    if (log_2) {
      use_y_label <- paste0("log2(", use_y_label, ")")
    }
    if (smoothed) {
      use_y_label <- paste0("smoothed(", use_y_label, ")")
    }
    use_y_label
  }

# gather relevant plot options from gene_info, outputs for ggplot
MakePlotOptionFrame <- function(gene_info) {
  print("plot options fun")
  # checks to see if at least one file in list is active
  if (gene_info %>% dplyr::filter(onoff != 0) %>% nrow() == 0) {
    return(NULL)
  } else {
    gene_info <- gene_info %>% dplyr::filter(onoff != 0)
    gene_info <- gene_info %>%
      dplyr::mutate(
        myline = 1,
        mydot = 0,
        mysizedot = 0.01,
        set = plot_set
      )
  }
  # tint if same color is used more then once
  ldf <- duplicated(gene_info$mycol)
  for (i in seq_along(gene_info$mycol)) {
    if (ldf[i]) {
      gene_info$mycol[i] <- RgbToHex(gene_info$mycol[i], convert = "hex", tint = log(i,10))
    }
  }
  return(gene_info)
}

# main ggplot function
GGplotLineDot <-
  function(list_long_data_frame,
           xBinRange,
           plot_options,
           yBinRange,
           line_list,
           use_smooth,
           plot_ttest,
           use_log2,
           use_y_label,
           plot_occupancy) {
    list_long_data_frame$set <- factor(list_long_data_frame$set, levels = plot_options$set)
    legend_space <- lengths(strsplit(sort(plot_options$set), "\n")) / 1.1
    if (use_log2) {
      gp <-
        ggplot(
          list_long_data_frame,
          aes(
            x = as.numeric(bin),
            y = log2(value),
            color = set,
            shape = set,
            size = set,
            linetype = set
          )
        )
    } else {
      gp <-
        ggplot(
          list_long_data_frame,
          aes(
            x = as.numeric(bin),
            y = value,
            color = set,
            shape = set,
            linetype = set
          )
        )
    }
    if (use_smooth) {
      gp <- gp +
        geom_smooth(se = FALSE,
                    size = line_list$mysize[2],
                    span = .2) 
    } else{
      gp <- gp +
        geom_line(size = line_list$mysize[2],alpha=0.8)
    }
    gp <- gp +
      scale_size_manual(values = plot_options$mysizedot) +
      scale_color_manual(values = plot_options$mycol) +
      scale_shape_manual(values = plot_options$mydot) +
      scale_linetype_manual(values = plot_options$myline) +
      ylab(use_y_label) +
      geom_vline(
        data = line_list$myline,
        aes(xintercept = use_virtical_line),
        size = line_list$mysize[1],
        linetype = line_list$myline$use_virtical_line_type,
        color = line_list$myline$use_virtical_line_color
      ) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()) +
      theme(axis.title.y = element_text(size =  line_list$mysize[4] + 4, margin = margin(2, 10, 2, 2))) +
      theme(axis.text.y = element_text(size = line_list$mysize[4],
                                       face = 'bold')) 
    if(!is_empty(LIST_DATA$ttest)){
      use_col_tt <- plot_ttest$options_main_tt$mycol
      use_line_tt <- plot_ttest$options_main_tt$myline
      names(use_col_tt) <- plot_ttest$options_main_tt$set
      names(use_line_tt) <- plot_ttest$options_main_tt$set
      gp2 <- ggplot(LIST_DATA$ttest, aes(y=p.value,x=bin,
                                         color=set,
                                         shape = set,
                                         linetype = set)) + 
        geom_line(size = line_list$mysize[2],alpha=0.8) +
        scale_color_manual(values = use_col_tt) +
        scale_linetype_manual(values = use_line_tt)+
        geom_hline(yintercept = plot_ttest$hlineTT,color="blue") + 
        theme_bw() +
        geom_vline(
          data = line_list$myline,
          aes(xintercept = use_virtical_line),
          size = line_list$mysize[1],
          linetype = line_list$myline$use_virtical_line_type,
          color = line_list$myline$use_virtical_line_color
        ) +
        xlab(paste(Sys.Date(), paste(unique(
          plot_options$sub
        ), collapse = ", "), collapse = ", ")) +
        ylab(plot_ttest$ylabTT)+
        scale_x_continuous(breaks = line_list$mybrakes[between(line_list$mybrakes, xBinRange[1], xBinRange[2])],
                           labels = line_list$mylabels[between(line_list$mybrakes, xBinRange[1], xBinRange[2])]) +
        theme(axis.title.x = element_text(size =  line_list$mysize[3], vjust = .5)) +
        theme(axis.title.y = element_text(size =  line_list$mysize[4], margin = margin(2, 10, 2, 2))) +
        theme(axis.text.y = element_text(size = line_list$mysize[4],
                                         face = 'bold'))+
        theme(
          axis.text.x = element_text(
            color = line_list$mycolors[between(line_list$mybrakes, xBinRange[1], xBinRange[2])],
            size = line_list$mysize[3],
            angle = -45,
            hjust = .1,
            vjust = .9,
            face = 'bold'
          )
        ) +
        theme(
          legend.title = element_blank(),
          legend.key = element_rect(size = line_list$mysize[5] / 2, color = 'white'),
          legend.key.height = unit(legend_space/1.2, "line"),
          legend.text = element_text(size = line_list$mysize[5]/1.2, face = 'bold')
        ) +
        coord_cartesian(xlim = xBinRange, ylim = plot_ttest$ylimTT)
      my_occupancy <- c(6 - plot_occupancy, plot_occupancy + .5)
      gp <- gp + coord_cartesian(xlim = xBinRange, ylim = unlist(yBinRange))
      suppressMessages(print(gp + gp2 + plot_layout(ncol = 1, heights = my_occupancy)))
      return(suppressMessages(gp+ gp2 + plot_layout(ncol = 1, heights = my_occupancy)))
    } else{
      gp <- gp + 
        xlab(paste(Sys.Date(), paste(unique(
          plot_options$sub
        ), collapse = ", "), collapse = ", ")) +
        scale_x_continuous(breaks = line_list$mybrakes[between(line_list$mybrakes, xBinRange[1], xBinRange[2])],
                           labels = line_list$mylabels[between(line_list$mybrakes, xBinRange[1], xBinRange[2])]) +
        theme(axis.title.x = element_text(size =  line_list$mysize[3], vjust = .5)) +
        theme(
          axis.text.x = element_text(
            #fix for coord_cartesian [between(line_list$mybrakes, xBinRange[1], xBinRange[2])]
            color = line_list$mycolors[between(line_list$mybrakes, xBinRange[1], xBinRange[2])],
            size = line_list$mysize[3],
            angle = -45,
            hjust = .1,
            vjust = .9,
            face = 'bold'
          )
        ) +
        theme(
          legend.title = element_blank(),
          legend.key = element_rect(size = line_list$mysize[5] / 2, color = 'white'),
          legend.key.height = unit(legend_space, "line"),
          legend.text = element_text(size = line_list$mysize[5], face = 'bold')
        )  +
        coord_cartesian(xlim = xBinRange, ylim = unlist(yBinRange))
      
      suppressMessages(print(gp))
      return(suppressMessages(gp))
    }
  }

# Sets lines and labels
LinesLabelsListset <- function(body1bin = 20,
                               body2bin = 40,
                               tssbin = 15,
                               tesbin = 45,
                               binbp = 100,
                               totbins = 80,
                               everybin = 5,
                               tssname = "TSS",
                               tesname = "pA") {
  # I create this in steps bin 1 to next land mark (TSS TES) then go from there to next land mark until end of bins
  everybp <- everybin * binbp
  if (everybp > 0) {
    my_5prim <- NULL
    my_3prim <- NULL
    if (tssbin > 0) {
      # set up 1 to TSS'
      TSSname <- seq(-tssbin * binbp, 0, by = everybp)
      TSSloc <- seq(1,  by = everybin, length.out = length(TSSname))
      # make sure TSS is included
      if (any(TSSname == 0)) {
        TSSloc[TSSname == 0] <- tssbin + .5
        TSSname[TSSname == 0] <- tssname
      } else if(any(TSSloc == tssbin)){
        TSSname[TSSloc == tssbin] <- tssname
        TSSloc[TSSloc == tssbin] <- tssbin + .5
      } else {
        TSSloc <- sort(c(TSSloc, tssbin + .5))
        TSSname <- append(TSSname, tssname)
      }
      my_5prim <-
        tibble(lloc = TSSloc, lname = as.character(TSSname))
      # keep distance between TSS and numbers
      nextStart <- tssbin + everybin
      if (body1bin > 0 & tesbin > 0 & body2bin > 0) {
        nextEnd <- body1bin
      } else if (tesbin == 0) {
        nextEnd <- totbins
      } else {
        nextEnd <- 0
      }
      if (nextStart <= nextEnd) {
        TSSname1 <- seq(everybp, (nextEnd - tssbin) * binbp, by = everybp)
        TSSloc1 <-
          seq(nextStart,
              by = everybin,
              length.out = length(TSSname1))
        # make sure body brake is included
        if (!any(TSSloc1 == nextEnd)) {
          TSSname1 <- append(TSSname1, (nextEnd - tssbin) * binbp)
          TSSloc1 <- c(TSSloc1, nextEnd)
        }
        my_5prim2 <-
          tibble(lloc = TSSloc1, lname = as.character(TSSname1))
        my_5prim <-
          full_join(my_5prim, my_5prim2, by = c("lloc", "lname"))
      } else if (nextEnd > 0) {
        my_5prim2 <-
          tibble(lloc = nextEnd, lname = as.character((nextEnd - tssbin) * binbp))
        my_5prim <-
          full_join(my_5prim, my_5prim2, by = c("lloc", "lname"))
      }
    }
    if (tesbin > 0) {
      # next to TES'
      if (body1bin > 0 & tssbin > 0 & body2bin > 0) {
        nextStart <- body2bin
      } else if (tssbin == 0) {
        nextStart <- 1
      } else {
        nextStart <- totbins
      }
      if (nextStart <= tesbin) {
        TESname <-  abs(seq((nextStart - tesbin) * binbp, 0, by = everybp))
        if(nextStart == 1){
          TESname <- TESname + binbp
        }
        TESloc <-
          seq(nextStart,
              by = everybin,
              length.out = length(TESname))
        # make sure TES is included
        if (any(TESname == 0)) {
          TESloc[TESname == 0] <- tesbin + .5
          TESname[TESname == 0] <- tesname
        } else if(any(TESloc == tesbin)){
          TESname[TESloc == tesbin] <- tesname
          TESloc[TESloc == tesbin] <- tesbin + .5
        }else {
          TESloc <- sort(c(TESloc, tesbin + .5))
          TESname <- append(TESname, tesname)
        }
        my_3prim <-
          tibble(lloc = TESloc, lname = as.character(TESname))
      } else {
        my_3prim <- tibble(lloc = tesbin + .5, lname = tesname)
      }
      # TES to end
      nextStart <- tesbin + everybin
      if (nextStart < totbins) {
        TESname1 <- seq(everybp, (totbins - tesbin) * binbp, by = everybp)
        TESloc1 <-
          seq(nextStart,
              by = everybin,
              length.out = length(TESname1))
        if (!any(TESloc1 == totbins)) {
          TESname1 <- append(TESname1, (totbins - tesbin) * binbp)
          TESloc1 <- c(TESloc1, totbins)
        }
        my_3prim2 <-
          tibble(lloc = TESloc1, lname = as.character(TESname1))
        my_3prim <-
          full_join(my_3prim, my_3prim2, by = c("lloc", "lname"))
      }
    }
    if (!is.null(my_5prim) & !is.null(my_3prim)) {
      my_53prim <-
        arrange(full_join(my_5prim, my_3prim, by = c("lloc", "lname")), lloc)
      use_plot_breaks <- my_53prim$lloc
      use_plot_breaks_labels <- my_53prim$lname
      use_plot_breaks_labels <-
        use_plot_breaks_labels[seq_along(use_plot_breaks)]
    } else if (!is.null(my_5prim)) {
      use_plot_breaks <- my_5prim$lloc
      use_plot_breaks_labels <- my_5prim$lname
      use_plot_breaks_labels <-
        use_plot_breaks_labels[seq_along(use_plot_breaks)]
    } else if (!is.null(my_3prim)) {
      use_plot_breaks <- my_3prim$lloc
      use_plot_breaks_labels <- my_3prim$lname
      use_plot_breaks_labels <-
        use_plot_breaks_labels[seq_along(use_plot_breaks)]
    } else {
      # just print bin numbers
      use_plot_breaks <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin) + 1)
      use_plot_breaks_labels <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin) + 1)
    }
    # no bp bin labels
  } else {
    if (everybin > 0) {
      use_plot_breaks <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin) + 1)
      use_plot_breaks_labels <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin) + 1)
    } else {
      use_plot_breaks <- .5
      use_plot_breaks_labels <- "none"
    }
  }
  # virtical line set up
  use_plot_breaks <- na_if(use_plot_breaks, 0.5)
  use_plot_breaks_labels <-
    use_plot_breaks_labels[!is.na(use_plot_breaks)]
  use_plot_breaks <- use_plot_breaks[!is.na(use_plot_breaks)]
  list(mybrakes = use_plot_breaks,
       mylabels = use_plot_breaks_labels)
}

# Sets plot lines and labels colors
LinesLabelsListPlot <-
  function(body1bin,
           body1color,
           body1line,
           body2bin,
           body2color,
           body2line,
           tssbin,
           tsscolor,
           tssline,
           tesbin,
           tescolor,
           tesline,
           use_plot_breaks_labels,
           use_plot_breaks,
           vlinesize,
           linesize,
           fontsizex,
           fontsizey,
           legendsize,
           binsize,
           binspace) {
    print("lines and labels plot fun")
    if (length(use_plot_breaks_labels) > 0) {
      mycolors <- rep("black", length(use_plot_breaks))
      use_virtical_line <- c(NA, NA, NA, NA)
      if (tssbin > 0) {
        if(tssbin > 1){
          mod <- 0.5
        } else {
          mod <- 0
        }
        mycolors[which(use_plot_breaks == tssbin  + mod)] <- tsscolor
        use_virtical_line[1] <- tssbin  + mod
        if (tssbin < body1bin &
            body1bin < body2bin &
            body2bin < tesbin & tesbin <= last(use_plot_breaks)) {
          use_virtical_line[3:4] <- c(body1bin, body2bin)
        }
      }
      if (tesbin > 0) {
        if(tesbin > 1){
          mod <- 0.5
        } else {
          mod <- 0
        }
        mycolors[which(use_plot_breaks == tesbin  + mod)] <- tescolor
        use_virtical_line[2] <- tesbin + mod
      }
    } else {
      use_plot_breaks <- mod
      use_plot_breaks_labels <- "none"
      use_virtical_line <- c(NA, NA, NA, NA)
    }
    # vertical line set up
    use_virtical_line_color <-
      c(tsscolor, tescolor, body1color, body2color)
    use_virtical_line_type <-
      c(tssline, tesline, body1line, body2line)
    use_plot_breaks <- na_if(use_plot_breaks, 0.5)
    use_virtical_line <- na_if(use_virtical_line, 0.5)
    use_plot_breaks_labels <-
      use_plot_breaks_labels[!is.na(use_plot_breaks)]
    use_plot_breaks <- use_plot_breaks[!is.na(use_plot_breaks)]
    use_virtical_line_type <-
      use_virtical_line_type[!is.na(use_virtical_line)]
    use_virtical_line_color <-
      use_virtical_line_color[!is.na(use_virtical_line)]
    use_virtical_line <-
      use_virtical_line[!is.na(use_virtical_line)]
    list(
      myline = virtical_line_data_frame <- data.frame(
        use_virtical_line,
        use_virtical_line_type,
        use_virtical_line_color,
        stringsAsFactors = FALSE
      ),
      mycolors = mycolors,
      mybrakes = use_plot_breaks,
      mylabels = use_plot_breaks_labels,
      mysize = c(vlinesize, linesize, fontsizex, fontsizey, legendsize),
      myset = c(body1bin, body2bin, tssbin, tesbin, binsize, binspace)
    )
  }

# lines and labels preset helper
LinesLabelsPreSet <- function(mytype) {
  # 5|4, 4|3, tss, pA, bp/bin, max bins, every bin
  if (mytype == kLinesandlabels[1]) {
    tt <- c(20, 40, 15, 45, 100, LIST_DATA$x_plot_range[2], 10)
  } else if (mytype == kLinesandlabels[2]) {
    tt <- c(0, 0, 40, 0, 25, LIST_DATA$x_plot_range[2], 20)
  } else if (mytype == kLinesandlabels[3]) {
    tt <- c(0, 0, 0, 10, 100, LIST_DATA$x_plot_range[2], 10)
  } else if (mytype == kLinesandlabels[4]) {
    tt <- c(0, 0, 5, 0, 50, LIST_DATA$x_plot_range[2], 10)
  } else if (mytype == kLinesandlabels[5]) {
    tt <- c(0,
            0,
            floor(LIST_DATA$x_plot_range[2] * .25),
            0,
            100,
            LIST_DATA$x_plot_range[2],
            10)
  } else if (mytype == kLinesandlabels[6]) {
    tt <- c(
      0,
      0,
      floor(LIST_DATA$x_plot_range[2] * .25),
      ceiling(LIST_DATA$x_plot_range[2] * .75),
      100,
      LIST_DATA$x_plot_range[2],
      10
    )
  } else if (mytype == kLinesandlabels[7]) {
    tt <- c(0,
            0,
            0,
            ceiling(LIST_DATA$x_plot_range[2] * .75),
            100,
            LIST_DATA$x_plot_range[2],
            10)
  } else {
    tt <-
      c(
        ceiling(LIST_DATA$x_plot_range[2] * .33),
        floor(LIST_DATA$x_plot_range[2] * .66),
        floor(LIST_DATA$x_plot_range[2] * .25),
        ceiling(LIST_DATA$x_plot_range[2] * .75),
        100,
        LIST_DATA$x_plot_range[2],
        ceiling(LIST_DATA$x_plot_range[2] * .1)
      )
  }
  tt
}

# records check box on/off
CheckBoxOnOff <- function(check_box, list_data) {
  if(!all(is.na(names(check_box)))){
    list_data <- full_join(list_data,check_box,by=c("set","gene_list")) %>%
      dplyr::filter(!is.na(set)) %>% 
      dplyr::mutate(onoff=if_else(is.na(onoff.y),"0",set)) %>% 
      dplyr::select(-onoff.y,-onoff.x) %>%
      distinct()
  }
  list_data
}

# applys t.test to active data
ApplyTtest <-
  function(list_data,
           switchttest,
           use_tmath,
           switchttesttype,
           padjust,my_alt, my_exact, my_paired) {
    # t.test comparing files in same gene list
    ttest <- NULL
    if(switchttest == "by lists" & n_distinct(list_data$gene_list) > 1){
      list_data <- list_data %>% 
        rename(set=gene_list,gene_list=set)
    }
    n_test <- list_data %>% group_by(gene_list) %>% 
      summarise(n=n_distinct(set), .groups= "drop")
    if(switchttest != "none"){
      for(i in distinct(list_data, gene_list)$gene_list){
        if(n_test %>% dplyr::filter(gene_list == i) %>% dplyr::select(n) > 1){
          ttest[[i]] <- bind_rows(try_t_test(list_data %>% dplyr::filter(gene_list == i),i,
                                             use_tmath,switchttesttype,padjust,
                                             my_alt, noquote(my_exact), noquote(my_paired))) %>% 
            dplyr::mutate(mydot = 0,
                          myline = 1,
                          mycol = "#000000" )
        } 
      }
    } 
    bind_rows(ttest)
  }

# makes sure t test wont crash on an error
try_t_test <- function(db,my_set,my_math ="none",my_test="t.test",padjust="fdr",
                       alternative="two.sided", exact=FALSE, paired=FALSE){
  print("try t.test")
  exact <- if_else(exact=="TRUE",TRUE,FALSE)
  paired <- if_else(paired=="TRUE",TRUE,FALSE)
  combn(unique(db$set),2) -> my_comparisons
  my_comparisons2 <- list()
  db_out <- list()
  for(cc in 1:ncol(my_comparisons)){
    my_comparisons2[[cc]] <- (c(my_comparisons[1,cc],my_comparisons[2,cc]))
  }
  db <- spread(db,set,score) 
  for(i in my_comparisons2){
    db2 <- dplyr::select(db, gene,bin,all_of(i)) %>% 
      rename(score.x=all_of(names(.)[3]), score.y=all_of(names(.)[4])) 
    
    myTtest <- tibble(bin=NA,p.value=NA)
    
    for(t in unique(db2$bin)){
      x.score <- dplyr::filter(db2, bin ==t)
      y.score <- dplyr::filter(db2, bin ==t)
      kk <- try(get(my_test)(x.score$score.x,y.score$score.y,
                             alternative = alternative, 
                             exact=exact, 
                             paired=paired)$p.value)
      if("try-error" %in% class(kk) | !is.numeric(kk)){
        kk <-1
      }
      myTtest <- myTtest %>% add_row(bin = t, p.value=kk)
    }
    myTtest <- myTtest %>% dplyr::filter(!is.na(bin))
    if(padjust != "NO"){
      myTtest <- myTtest %>% dplyr::mutate(p.value=p.adjust(p.value,method = padjust))
    }
    if(my_math =="-log"){
      myTtest <- myTtest %>% dplyr::mutate(p.value=if_else(p.value==0,2.2e-16,p.value)) %>% dplyr::mutate(p.value=-log(p.value))
    } else if(my_math =="-log10"){
      myTtest <- myTtest %>% dplyr::mutate(p.value=if_else(p.value==0,2.2e-16,p.value)) %>% dplyr::mutate(p.value=-log10(p.value))
    }
    db_out[[str_c(i,collapse = "-")]] <- myTtest%>% 
      dplyr::mutate(., set = paste(
        paste0(gsub("(.{20})", "\\1\n", 
                    str_split_fixed(str_c(i,collapse = "-"), "\n",n=2)[,1]),
               gsub("(.{20})", "\\1\n", 
                    str_split_fixed(str_c(i,collapse = "-"), "\n",n=2)[,2])),
        gsub("(.{20})", "\\1\n", 
             str_split_fixed(my_set, "\n", n=2)[,1]),
        str_split_fixed(my_set, "\n", n=2)[,2],
        sep = '\n'
      ))
    
  }
  
  db_out
}

# gather relevant plot option data
MakePlotOptionttest <- function(list_data, Y_Axis_TT,my_ttest_log,hlineTT,pajust,ttype) {
  if(is_empty(list_data)){
    return(NULL)
  }
  print("plot options ttest fun")
  out_options <- list_data %>% select(-bin, -p.value) %>% distinct(set,mydot,myline,mycol)
  list_data_frame <- NULL
  
  ldf <- duplicated(out_options$mycol)
  for (i in seq_along(out_options$mycol)) {
    if (ldf[i]) {
      out_options$mycol[i] <-
        RgbToHex(out_options$mycol[i], convert = "hex", tint = log(i,10))
    }
  }
  
  list_data_frame$options_main_tt <- out_options
  list_data_frame$ylimTT <- Y_Axis_TT
  if(pajust != "none"){
    pp <- paste0("p.adjust:", pajust)
  } else{
    pp <- "p.value"
  }
  if(my_ttest_log == "-log"){
    list_data_frame$hlineTT <- -log(hlineTT)
    list_data_frame$ylabTT <- paste0(my_ttest_log,"(",pp,")"," ",ttype)
  } else if(my_ttest_log == "-log10"){
    list_data_frame$hlineTT <- -log10(hlineTT)
    list_data_frame$ylabTT <- paste0(my_ttest_log,"(",pp,")"," ",ttype)
  } else{
    list_data_frame$hlineTT <- hlineTT
    list_data_frame$ylabTT <- paste0(pp," ",ttype)
  }
  
  return(list_data_frame)
}

# YaxisValuetTest gathers values
YaxisValuetTest <- function(ttest, hlinettest, selectttestlog ){
  mm <- round(extendrange(range(ttest$p.value, na.rm = T,finite=T),f = .1),digits = 2)
  p_cutoff <- hlinettest
    if(selectttestlog == "-log"){
      p_cutoff <- -log(hlinettest)
    } else if(selectttestlog == "-log10"){
      p_cutoff <- -log10(hlinettest)
    }
    if(mm[1] > 0){
      mm[1] <- 0
    }
    if(mm[2] < p_cutoff){
      mm[2] <- p_cutoff
    }
  mm <- c(round(min(mm), 4),round(max(mm), 4))
  return(mm)
}

# sorts active gene list contain top % signal based on selected bins and file
FilterTop <-
  function(list_data,
           list_name,
           file_names,
           start_bin,
           end_bin,
           mynum,
           topbottom,
           my_filter_all) {
    if (is.null(file_names)) {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return(NULL)
    }
    lc <- 0
    outlist <- NULL
    lapply(file_names, function(j) {
      setProgress(lc + 1, detail = paste("sorting", j))
      apply_bins <-
        semi_join(dplyr::filter(list_data$table_file, set == j), 
                  list_data$gene_file[[list_name]]$use, by = 'gene')  
      apply_bins <- group_by(apply_bins, gene) %>%
        dplyr::filter(bin %in% start_bin:end_bin) %>%
        summarise(mysums = sum(score, na.rm = TRUE),.groups="drop") %>%
        mutate(myper = as.numeric(strtrim(cume_dist(mysums), 5))) %>%
        arrange(desc(mysums))
      gene_count <- nrow(apply_bins)
      if (topbottom == "Top%") {
        num2 <- c(1, ceiling(gene_count * (mynum / 100)))
      } else if (topbottom == "Middle%") {
        if (mynum == 100) {
          num2 <- c(1, gene_count)
        } else {
          med <- median(apply_bins$myper)
          num2 <-
            c(count(apply_bins, myper >= max(med, mynum / 100))[[2]][2],
              count(apply_bins, myper <= min(med, (100 - mynum) / 100))[[2]][1])
        }
      } else {
        num2 <-
          c(ceiling((gene_count + 1) - (gene_count * (mynum / 100))), gene_count)
      }
      if (any(is.na(num2))) {
        num2 <-
          c(ceiling((gene_count) - (gene_count * max(.5, mynum / 100))), ceiling(gene_count * max(.5, mynum / 100)))
      }
      outlist2 <- dplyr::mutate(apply_bins,!!j := myper) %>%
        dplyr::select(gene,!!j) %>%
        slice(num2[1]:num2[2])
      if (lc > 0) {
        if(my_filter_all){
          outlist <<- inner_join(outlist, outlist2, by = 'gene')
        } else{
          outlist <<- full_join(outlist, outlist2, by = 'gene')
        }
      } else {
        outlist <<- outlist2
      }
      lc <<- lc + 1
    })
    if (length(outlist$gene) == 0) {
      return(NULL)
    }
    old_name <- grep("^Filter", names(list_data$gene_file), value = T)
    if (length(old_name) > 3) {
      # remove old sort gene list keepint 4
      list_data$gene_file[[first(old_name)]] <- NULL
      list_data$gene_info <- dplyr::filter(list_data$gene_info,
                                           gene_list != first(old_name))
    }
    setProgress(lc + 2, detail = "building list")
    topbottom2 <- paste(str_remove(topbottom,"%"), paste0(mynum, "%"))
    nick_name <-
      strtrim(gsub("(.{30})",
                   "\\1... ",
                   paste0("Filter ",topbottom2, "\nn = ", n_distinct(outlist$gene))), 33)
    list_data$gene_file[[nick_name]]$full <- outlist
    list_data$gene_file[[nick_name]]$use <- dplyr::select(outlist, gene)
    list_data$gene_file[[nick_name]]$info <- tibble(loaded_info =
      paste(
        "Filter",
        topbottom2,
        "bins",
        start_bin,
        "to",
        end_bin,
        "from",
        list_name,
        paste(file_names, collapse = " "),
        Sys.Date(),
        list_data$gene_file[[list_name]]$info
      ),
      matching = FALSE)
    list_data$gene_info <- 
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>% 
                           dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                           dplyr::mutate(gene_list = nick_name,
                                         sub =  paste("Filter",
                                                      topbottom2,
                                                      "bins",
                                                      start_bin,
                                                      "to",
                                                      end_bin), 
                                         onoff = "0",
                                         count = paste0("n = ", n_distinct(outlist$gene)),
                                         plot_set = " ")))
    list_data
  }

# sort my percent
FilterPer <-
  function(list_data,
           list_name,
           file_names,
           start_end_bin,
           my_per,
           my_type,
           anyall) {
    if (is.null(file_names)) {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return(NULL)
    }
    p_funs <- map(my_per/100, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
      set_names(paste0("my_p_",seq_along(my_per)))
    gene_list <- list_data$gene_file[[list_name]]$use
    out_list <- list_data$table_file %>% 
      dplyr::filter(set %in% file_names) %>% 
      semi_join(.,gene_list,by="gene") %>% 
      dplyr::filter(bin %in% start_end_bin[1]:start_end_bin[2]) 
    
    out_per <- out_list %>%
      group_by(set,bin) %>% summarize_at(vars(score), p_funs) %>% ungroup() 
    
    if(my_type == "min%"){
      if(anyall){
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene) %>% dplyr::filter(all(score >= my_p_1)) %>% 
          ungroup() %>% distinct(gene)
      } else {
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene,set) %>% dplyr::filter(all(score >= my_p_1)) %>% 
          ungroup() %>% distinct(gene)
      }
      topbottom2 <- paste(str_remove(my_type,"%"), paste0(my_per[1], "%"))
    } else if(my_type == "max%"){
      if(anyall){
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene) %>% dplyr::filter(all(score <= my_p_2)) %>% 
          ungroup() %>% distinct(gene)
      } else {
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene,set) %>% dplyr::filter(all(score <= my_p_2)) %>% 
          ungroup() %>% distinct(gene)
      }
      topbottom2 <- paste(str_remove(my_type,"%"), paste0(my_per[2], "%"))
    } else {
      if(anyall){
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene) %>% dplyr::filter(all(score >= my_p_1 & score <= my_p_2)) %>% 
          ungroup() %>% distinct(gene)
      } else {
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene,set) %>% dplyr::filter(all(score >= my_p_1 & score <= my_p_2)) %>% 
          ungroup() %>% distinct(gene)
      }
      topbottom2 <- paste(paste(str_remove(my_type,"%"), paste0(my_per[1], "%")),paste0(my_per[2], "%"),collapse = " and ")
    }
    if (length(out_list$gene) == 0) {
      return(NULL)
    }
    old_names <- grep("^Filter", names(list_data$gene_file), value = T)
    if (length(old_names) > 3) {
      # remove old sort gene list keeping 4
      list_data$gene_file[[first(old_names)]] <- NULL
      list_data$gene_info <- dplyr::filter(list_data$gene_info,
                                           gene_list != first(old_names))
    }
    nick_name <-
      strtrim(gsub("(.{30})",
                   "\\1... ",
                   paste0("Filter Prob ",topbottom2, "\nn = ", n_distinct(out_list$gene))), 33)
    if(length(file_names) == 1){
      my_per <- seq(my_per[1],my_per[2]/2.1,length.out = 5)
      p_funs <- map(my_per/100, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
        set_names(my_per)
      my_list <- list_data$table_file %>% 
        dplyr::filter(set == file_names) %>% 
        semi_join(.,gene_list,by="gene") %>% 
        dplyr::filter(bin %in% start_end_bin[1]:start_end_bin[2]) 
      out_per1 <- my_list %>%
        group_by(set,bin) %>% summarize_at(vars(score), p_funs) %>% ungroup() %>% 
        gather(.,key = set,value = "my_p_1",-bin,-set)
      
      my_per <- seq(my_per[2]/2,my_per[2],length.out = 5)
      p_funs <- map(my_per/100, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
        set_names(my_per)
      
      out_per2 <- my_list %>%
        group_by(set,bin) %>% summarize_at(vars(score), p_funs) %>% ungroup() %>% 
        gather(.,key = set,value = "my_p_2",-bin,-set)
      out_list1 <- full_join(out_per1,out_per2,by=c("bin","set")) %>% 
        dplyr::mutate(set=paste0(round(as.numeric(set),2),"%"))
    } else {
      out_list1 <- list_data$table_file %>% 
        dplyr::filter(set %in% file_names) %>% 
        semi_join(.,gene_list,by="gene") %>% 
        full_join(.,out_per,by=c("bin","set")) %>% 
        replace_na(list(my_p_1 = 0, my_p_2 = 0)) %>% 
        dplyr::select(-gene,-score) 
    }
    
    list_data$gene_file[[nick_name]]$full <- out_list %>% dplyr::mutate(min=my_per[1],max=my_per[2])
    list_data$gene_file[[nick_name]]$use <- out_list
    list_data$sortplot <- out_list1
    list_data$gene_file[[nick_name]]$info <- tibble(loaded_info =
      paste(
        "Filter Prob:",
        topbottom2,
        "bins",
        start_end_bin[1],
        "to",
        start_end_bin[2],
        "from",
        list_name,
        paste(file_names, collapse = " "),
        Sys.Date(),
        list_data$gene_file[[list_name]]$info
      ),
      matching = FALSE)
    list_data$gene_info <-
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>%
                           dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>%
                           dplyr::mutate(gene_list = nick_name,
                                         sub =  paste("Filter Prob:",
                                                      topbottom2,
                                                      "bins",
                                                      start_end_bin[1],
                                                      "to",
                                                      start_end_bin[2]),
                                         onoff = "0",
                                         count = paste0("n = ", n_distinct(out_list$gene)),
                                         plot_set = " ")))
    
    list_data
  }

# filter peaks
FilterPeak <-
  function(list_data,
           list_name,
           file_names,
           start_end_bin_peak,
           start_end_bin_filter,
           my_type,
           anyall) {
    if (is.null(file_names)) {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return(NULL)
    }
    gene_list <- list_data$gene_file[[list_name]]$use
    out_list <- list_data$table_file %>% 
      dplyr::filter(set %in% file_names) %>% 
      semi_join(.,gene_list,by="gene") 
    
    my_filter <- out_list %>% 
      dplyr::filter(bin %in% start_end_bin_peak[1]:start_end_bin_peak[2]) %>% 
      group_by(set,gene) %>% summarise(score2=max(score,rm.na=T),.groups = "drop")
    out_list <- out_list %>% 
      dplyr::filter((bin %in% start_end_bin_filter[1]:start_end_bin_filter[2])) %>% 
      full_join(.,my_filter,by=c("gene","set")) 
    
    if(my_type == "peak"){
      if(anyall){
        out_gene <- out_list %>% 
          group_by(gene) %>% 
          dplyr::filter(all(score<=score2)) %>% 
          ungroup() %>% distinct(gene)
      } else {
        out_gene <- out_list %>%
          group_by(gene,set) %>% 
          dplyr::filter(all(score<=score2)) %>% 
          ungroup() %>% distinct(gene)
      }
    } else {
      if(anyall){
        out_gene <- out_list %>% 
          group_by(gene) %>% 
          dplyr::filter(!all(score<=score2)) %>% 
          ungroup() %>% distinct(gene)
      } else {
        out_gene <- out_list %>%
          group_by(gene,set) %>% 
          dplyr::filter(!all(score<=score2)) %>% 
          ungroup() %>% distinct(gene)
      }
    }
    if (length(out_gene$gene) == 0) {
      return(NULL)
    }
    old_names <- grep("^Filter", names(LIST_DATA$gene_file), value = T)
    if (length(old_names) > 3) {
      # remove old sort gene list keeping 4
      list_data$gene_file[[first(old_names)]] <- NULL
      list_data$gene_info <- dplyr::filter(list_data$gene_info,
                                           gene_list != first(old_names))
    }
    nick_name <-
      strtrim(gsub("(.{30})",
                   "\\1... ",
                   paste0("Filter ",my_type, "\nn = ", n_distinct(out_gene$gene))), 33)
    list_data$gene_file[[nick_name]]$full <- out_gene
    list_data$gene_file[[nick_name]]$use <- out_gene
    list_data$gene_file[[nick_name]]$info <- tibble(loaded_info =
                                                      paste(
                                                        "Filter peak:",
                                                        my_type,
                                                        "bins",
                                                        start_end_bin_peak[1],
                                                        "to",
                                                        start_end_bin_peak[2],
                                                        "filter bins",
                                                        start_end_bin_filter[1],
                                                        "to",
                                                        start_end_bin_filter[2],
                                                        "from",
                                                        list_name,
                                                        paste(file_names, collapse = " "),
                                                        Sys.Date(),
                                                        list_data$gene_file[[list_name]]$info
                                                      ),
                                                    matching = FALSE)
    list_data$gene_info <-
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>%
                           dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>%
                           dplyr::mutate(gene_list = nick_name,
                                         sub =  paste("Filter peak:",
                                                      my_type,
                                                      "peak bins"),
                                         onoff = "0",
                                         count = paste0("n = ", n_distinct(out_gene$gene)),
                                         plot_set = " ")))
    
    list_data
  }

# make a new normalized file by dividing one file by the other
MakeNormFile <-
  function(list_data,
           nom,
           dnom,
           gbyg,
           divzerofix,
           addfiles,
           nickname) {
    # check 2 files have been selected
    if (nom == "" | dnom == "") {
      return(NULL)
    }
    # set up tool info and progress bar
    myname <- "bin_by_bin"
    # get data files
    nd <- list_data$table_file %>% dplyr::filter(set == nom | set== dnom) %>% 
      replace_na(., list(score = 0))
    if(nchar(nickname)<1){
      nickname <- paste(nom, addfiles, dnom,sep = " ")
    }
    # if min/2 find Na's and 0's, and replace
    if (divzerofix) {
      nd <- nd %>% 
        dplyr::mutate(score=if_else(set == dnom & score == 0, NA_real_, score))
      myname <- paste0(myname, "_0->min/2")
      new_min_for_dom <-
        min(nd$score, na.rm = TRUE) / 2
      nd <-
        replace_na(nd, list(score = new_min_for_dom))
    }
    # files numbers are replaced with mean of bins if applied
    if (gbyg != "bin by bin") {
      myname <- "mean_of_bins"
      if (divzerofix) {
        myname <- paste0(myname, "_0->min/2")
      }
      nd <- nd %>% 
        group_by(bin, set) %>%
        dplyr::mutate(score = mean(score, na.rm = TRUE)) %>% ungroup()
    }
    # applies custom norm factor(s)
    legend_nickname <- nickname
    if (addfiles == "+") {
      new_gene_list <- full_join(dplyr::filter(nd,set == nom), 
                                 dplyr::filter(nd,set == dnom), by = c("gene", "bin")) %>% 
        replace_na(., list(score = 0))
      new_gene_list <- transmute(
        new_gene_list,
        gene = gene,
        bin = bin,
        set = legend_nickname,
        score = score.x + score.y
      )
    } else if(addfiles == "-"){
      new_gene_list <- full_join(dplyr::filter(nd,set == nom), 
                                 dplyr::filter(nd,set == dnom), by = c("gene", "bin")) %>% 
        replace_na(., list(score = 0))
      new_gene_list <- transmute(
        new_gene_list,
        gene = gene,
        bin = bin,
        set = legend_nickname,
        score = score.x - score.y
      )
    } else {
      
      new_gene_list <- full_join(dplyr::filter(nd,set == nom), 
                                 dplyr::filter(nd,set == dnom), by = c("gene", "bin")) %>% 
        replace_na(., list(score = 0))
      # make gene list and do math
      legend_nickname <- paste0(nickname, ": ", myname)
      new_gene_list <- transmute(
        new_gene_list,
        gene = gene,
        bin = bin,
        set = legend_nickname,
        score = score.x / score.y
      ) %>% na_if(Inf)
    }
    # output test
    if (n_distinct(new_gene_list$gene) < 1) {
      showModal(
        modalDialog(
          title = "Information message",
          " No genes left, try replacing Inf and/or bin by bin",
          size = "s",
          easyClose = TRUE
        )
      )
      return(NULL)
    }
    # adds meta data 
    list_data$table_file <- dplyr::filter(list_data$table_file, set != legend_nickname)
    list_data$gene_info <- dplyr::filter(list_data$gene_info, set != legend_nickname)
    list_data$table_file <- bind_rows(list_data$table_file, new_gene_list)
    list_data$gene_info <- bind_rows(list_data$gene_info,
                                     tibble(
                                       gene_list = names(list_data$gene_file)[1],
                                       set = legend_nickname,
                                       mycol = sample(suppressWarnings(brewer.pal(11, sample(kBrewerList,size=1)))[-c(4:7)],size = 1),
                                       onoff = "0",
                                       sub = " ",
                                       plot_set = " "
                                     ))
    return(list_data)
  }

# Total, exclusive and intersected gene lists
IntersectGeneLists <-
  function(list_data, list_name) {
    if (is.null(list_name)) {
      return(NULL)
    }
    setProgress(1, detail = paste("building list"))
    outlist <- NULL
    # grab selected gene list(s)
    lapply(list_name, function(j) {
      outlist[[j]] <<- list_data$gene_file[[j]]$use %>% dplyr::select(gene)
    })
    # collapses into one list
    outlist <- bind_rows(outlist)
    # remove any pre used data
    list_data$gene_file <- list_data$gene_file[!str_detect(names(list_data$gene_file),"^Gene_List_")]
    list_data$gene_info <- list_data$gene_info %>% dplyr::filter(!str_detect(gene_list,"^Gene_List_"))
    
    # record for info
    setProgress(2, detail = paste("building Total list"))
    if (n_distinct(outlist$gene) > 0) {
      nick_name1 <-
        paste("Gene_List_Total\nn =", n_distinct(outlist$gene))
      # record for info
      list_data$gene_file[[nick_name1]]$full <- distinct(outlist)
      list_data$gene_file[[nick_name1]]$use <- distinct(dplyr::select(outlist, gene))
      list_data$gene_file[[nick_name1]]$info <- tibble(loaded_info =
        paste("Gene_List_Total",
              "from",
              paste(list_name, collapse = " and "),
              Sys.Date()),
        matching = FALSE)
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste("Gene_List_Total"), 
                                           onoff = "0",
                                           plot_set = " ")))
    }
    setProgress(3, detail = paste("building intersect list"))
    intersected <- dplyr::filter(outlist, duplicated(gene))
    if (n_distinct(intersected$gene) > 0) {
      nick_name1 <-
        paste("Gene_List_intersect\nn =", n_distinct(intersected$gene))
      # record for info
      list_data$gene_file[[nick_name1]]$full <- intersected
      list_data$gene_file[[nick_name1]]$use <- dplyr::select(intersected, gene)
      list_data$gene_file[[nick_name1]]$info <- tibble(loaded_info =
        paste("Gene_List_intersect",
              "from",
              paste(list_name, collapse = " and "),
              Sys.Date()),
      matching = FALSE)
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste("Gene_List_intersect"), 
                                           onoff = "0",
                                           plot_set = " ")))
      setProgress(4, detail = paste("building exclusive list"))
      exclusive <- anti_join(outlist, intersected, by = "gene")
      if (n_distinct(exclusive$gene) == 0) {
        exclusive <- distinct(outlist)
      }
      nick_name1 <-
        paste("Gene_List_exclusive\nn =", n_distinct(exclusive$gene))
      # record for info
      list_data$gene_file[[nick_name1]]$full <- exclusive
      list_data$gene_file[[nick_name1]]$use <- dplyr::select(exclusive, gene)
      list_data$gene_file[[nick_name1]]$info <- tibble(loaded_info =
        paste("Gene_List_exclusive",
              "from",
              paste(list_name, collapse = " and "),
              Sys.Date()),
        matching = FALSE)
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste("Gene_List_exclusive"), 
                                           onoff = "0",
                                           plot_set = " ")))
      
    }
    list_data
  }

# creates clusters using hclust.vector, from one active file
FindClusters <- function(list_data,
                         list_name,
                         clusterfile,
                         start_bin,
                         end_bin) {
  if (clusterfile == "") {
    showModal(modalDialog(
      title = "Information message",
      paste("No file selected to work on"),
      size = "s",
      easyClose = TRUE
    ))
    return(NULL)
  }
  setProgress(1, detail = paste("gathering data"))
  df <-
      semi_join(dplyr::filter(list_data$table_file, set == clusterfile), 
                list_data$gene_file[[list_name]]$use, by = 'gene') 
    setProgress(2, detail = "hierarchical clustering using ward method")
    list_data$clust <- list()
    list_data$clust$cm <-
      hclust.vector(as.data.frame(spread(df, bin, score))[, c((start_bin:end_bin) + 2)], method = "ward")
  list_data$clust$use <- distinct(df, gene)
  list_data
}

# Pull out the number of groups of clusters 2-10
ClusterNumList <- function(list_data,
                           list_name,
                           clusterfile,
                           start_bin,
                           end_bin,
                           num) {
  if (is_empty(list_data$clust) | clusterfile == "") {
    showModal(modalDialog(
      title = "Information message",
      paste("No file selected to work on"),
      size = "s",
      easyClose = TRUE
    ))
    return(NULL)
  }
  setProgress(3, detail = "spliting into clusters")
  list_data$gene_file <- list_data$gene_file[!str_detect(names(list_data$gene_file),"^Cluster_")]
  list_data$gene_info <- list_data$gene_info %>% dplyr::filter(!str_detect(gene_list,"^Cluster_"))
  gene_list <-
      dplyr::mutate(list_data$clust$use, cm = cutree(list_data$clust$cm, num))
  for (nn in 1:num) {
    outlist <- dplyr::filter(gene_list, cm == nn)
    nick_name <-
      paste(paste0("Cluster_", nn, "\nn ="), n_distinct(outlist$gene))
    list_data$gene_file[[nick_name]]$full <- outlist
    list_data$gene_file[[nick_name]]$use <- dplyr::select(outlist, gene)
    list_data$gene_file[[nick_name]]$info <- tibble(loaded_info =
      paste(
        nick_name,
        "bins",
        start_bin,
        "to",
        end_bin,
        "from",
        list_name,
        clusterfile,
        num,
        "Cluster_",
        "total",
        Sys.Date()
      ),
    matching = FALSE)
    list_data$gene_info <- 
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>% 
                           dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                           dplyr::mutate(gene_list = nick_name,
                                         sub =  paste(
                                           "Cluster_",
                                           "by",
                                           num,
                                           "from",
                                           list_name
                                         ), 
                                         onoff = "0",
                                         count = paste0("n = ", n_distinct(outlist$gene)),
                                         plot_set = " ")))
    setProgress(4, detail = paste("finishing cluster", nn))
  }
  list_data
}

# a[1]/b[2] or (a[1]/a[2])/(b[1]/b[2]) make gene list
CompareRatios <-
  function(list_data,
           list_name,
           ratio1file,
           ratio2file,
           start1_bin,
           end1_bin,
           start2_bin,
           end2_bin,
           num,
           divzerofix,
           normbin = 0) {
    if (ratio1file == "") {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return()
    }
    outlist <- NULL
    if (ratio2file == "None" | ratio2file == "") {
      if(sum(start2_bin,end2_bin) == 0){
        showModal(modalDialog(
          title = "Information message",
          paste("no file or bins to compare to"),
          size = "s",
          easyClose = TRUE
        ))
        return()
      }
      ratiofile <- ratio1file
    } else {
      ratiofile <- c(ratio1file, ratio2file)
    }
    lc <- 0
    lapply(ratiofile, function(j) {
      df <-
        semi_join(dplyr::filter(list_data$table_file, set == j), 
                  list_data$gene_file[[list_name]]$use, by = 'gene') 
      if (normbin > 0) {
        df <- group_by(df, gene) %>%
          arrange(bin) %>% 
          dplyr::mutate(score = score / nth(score, normbin))
      }
      df <- group_by(df, gene) %>%
        summarise(sum1 = sum(score[start1_bin:end1_bin],	na.rm = T),
                  sum2 = sum(score[start2_bin:end2_bin],	na.rm = T),.groups="drop") %>%
        ungroup()
      # if min/2 find Na's and 0's, and replace
      if(sum(start2_bin,end2_bin) == 0){
        df$sum2 <- 1
      }
      if (divzerofix) {
        df$sum2 <- na_if(df$sum2, 0)
        new_min <-
          min(df$sum2, na.rm = TRUE) / 2
        df <-
          replace_na(df, list(sum2 = new_min))
      }
      lc <<- lc + 1
      outlist[[lc]] <<-
        transmute(df, gene = gene, Ratio = sum1 / sum2) %>%
        na_if(Inf) %>% dplyr::select(gene, Ratio)
      
      if (lc > 1) {
        if (divzerofix) {
          outlist[[2]]$Ratio <- na_if(outlist[[2]]$Ratio, 0)
          outlist[[2]] <-
            replace_na(outlist[[2]], list(Ratio = new_min))
          outlist[[1]] <<-
            inner_join(outlist[[1]], outlist[[2]], by = 'gene') %>%
            transmute(gene = gene, Ratio = Ratio.x / Ratio.y) %>%
            na_if(Inf)  %>% dplyr::select(gene, Ratio)
        } else {
          outlist[[1]] <<-
            inner_join(outlist[[1]], outlist[[2]], by = 'gene') %>%
            transmute(gene = gene, Ratio = Ratio.x / Ratio.y) %>%
            na_if(Inf)  %>% dplyr::select(gene, Ratio)
        }
      }
    })
    #remove old info
    list_data$gene_file <- list_data$gene_file[!str_detect(names(list_data$gene_file),"^Ratio_")]
    list_data$gene_info <- list_data$gene_info %>% dplyr::filter(!str_detect(gene_list,"^Ratio_"))
    if(num < 0){
      num <- 1/num
    }
    nick_name <- NULL
    if(num != 0){
      upratio <- dplyr::filter(outlist[[1]], Ratio > num)
    } else {
      upratio <- NULL
    }
    
    if (n_distinct(upratio$gene) > 0) {
      nick_name1 <-
        paste("Ratio_Up_file1\nn =", n_distinct(upratio$gene))
      nick_name <- c(nick_name, nick_name1)
      list_data$gene_file[[nick_name1]]$full <- upratio %>% dplyr::mutate(set=nick_name1)
      list_data$gene_file[[nick_name1]]$use <- dplyr::select(upratio, gene)
      list_data$gene_file[[nick_name1]]$info <- tibble(loaded_info =
        paste(
          "Ratio_Up_file1",
          ratio2file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "]/",
          ratio1file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "] ",
          "fold change cut off",
          num,
          divzerofix,
          "from",
          list_name,
          "gene list",
          Sys.Date()
        ),
        matching = FALSE)
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste(
                                             "\nRatio_Up_file1",
                                             "fold change cut off",
                                             num,
                                             divzerofix,
                                             "from",
                                             list_name
                                           ), 
                                           onoff = "0",
                                           count = paste0("n = ", n_distinct(upratio$gene)),
                                           plot_set = " ")))
    }
    if(num != 0){
      upratio <- dplyr::filter(outlist[[1]], Ratio < 1 / num & Ratio != 0)
    } else {
      upratio <- NULL
    } 
    if (n_distinct(upratio$gene) > 0) {
      nick_name2 <-
        paste("Ratio_Down_file1\nn =", n_distinct(upratio$gene))
      nick_name <- c(nick_name, nick_name2)
      list_data$gene_file[[nick_name2]]$full <- upratio %>% dplyr::mutate(set=nick_name2)
      list_data$gene_file[[nick_name2]]$use <- dplyr::select(upratio, gene)
      list_data$gene_file[[nick_name2]]$info <- tibble(loaded_info =
        paste(
          "Ratio_Up_file1",
          ratio2file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "]/",
          ratio1file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "] ",
          "fold change cut off",
          num,
          divzerofix,
          "from",
          list_name,
          "gene list",
          Sys.Date()
        ),
        matching = FALSE)
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name2,
                                           sub =  paste(
                                             "Ratio_Down_file1",
                                             "fold change cut off",
                                             num,
                                             divzerofix,
                                             "from",
                                             list_name
                                           ), 
                                           onoff = "0",
                                           count = paste0("n = ", n_distinct(upratio$gene)),
                                           plot_set = " ")))
    }
    if(num != 0){
      upratio <-
        dplyr::filter(outlist[[1]], Ratio <= num &
                        Ratio >= 1 / num | Ratio == 0)
    } else {
      upratio <- outlist[[1]]
    }
    
    if (n_distinct(upratio$gene) > 0) {
      nick_name3 <-
        paste("Ratio_No_Diff\nn =", n_distinct(upratio$gene))
      nick_name <- c(nick_name, nick_name3)
      list_data$gene_file[[nick_name3]]$full <- upratio %>% dplyr::mutate(set=nick_name3)
      list_data$gene_file[[nick_name3]]$use <- dplyr::select(upratio, gene)
      list_data$gene_file[[nick_name3]]$info <- tibble(loaded_info =
        paste(
          "Ratio_Up_file1",
          ratio2file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "]/",
          ratio1file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "] ",
          "fold change cut off",
          num,
          divzerofix,
          "from",
          list_name,
          "gene list",
          Sys.Date()
        ),
      matching = FALSE)
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name3,
                                           sub =  paste(
                                             "Ratio_No_Diff",
                                             "fold change cut off",
                                             num,
                                             divzerofix,
                                             "from",
                                             list_name
                                           ), 
                                           onoff = "0",
                                           count = paste0("n = ", n_distinct(upratio$gene)),
                                           plot_set = " ")))
    }
    list_data$boxRatio <- NULL
    for (nn in nick_name) {
      list_data$boxRatio <- bind_rows(list_data$boxRatio,list_data$gene_file[[nn]]$full)
    }
    list_data
  }

# Cumulative Distribution data prep "PI / EI"
CumulativeDistribution <-
  function(list_data,
           onoffs,
           start1_bin,
           end1_bin,
           start2_bin,
           end2_bin) {
    print("cdf function")
    if (is.null(onoffs)) {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return()
    }
    # remove old data sets
    list_data$gene_file <- list_data$gene_file[!str_detect(names(list_data$gene_file),"^CDF")]
    list_data$gene_info <- list_data$gene_info %>% dplyr::filter(!str_detect(gene_list,"^CDF"))
    outlist <- NULL
    for (list_name in names(onoffs)) {
      setProgress(1, detail = paste("dividing one by the other"))
      # Complete within gene list and sum regions
      tf <- dplyr::filter(list_data$table_file, set %in% onoffs[[list_name]])
      gene_common <- tf %>% group_by(set) %>% distinct(gene) %>% ungroup()
      gene_common <- gene_common %>% 
        group_by(gene) %>% 
        filter(n_distinct(set)==n_distinct(gene_common$set)) %>% 
        distinct(gene) %>% ungroup()
      tf <- tf %>% 
        semi_join(., gene_common, by = "gene")
      outlist[[list_name]] <-
        semi_join(tf, 
                  list_data$gene_file[[list_name]]$use, by = 'gene') %>% 
        group_by(gene,set) %>%
        summarise(sum1 = mean(score[start1_bin:end1_bin],	na.rm = T),
                  sum2 = mean(score[start2_bin:end2_bin],	na.rm = T),.groups="drop") %>%
        dplyr::mutate(., value = sum1 / sum2) %>%
        dplyr::mutate(value=log2(value)) %>% 
        na_if(Inf) %>% na_if(-Inf) %>% group_by(gene) %>% 
        dplyr::filter(!any(is.na(value))) %>% ungroup() %>% 
        group_by(., set) %>%
        arrange(value) %>%
        dplyr::transmute(
          gene=gene,
          bin = row_number(),
          set = set,
          plot_set = paste(list_name, "-", set),
          value = value
        ) %>%
        ungroup()
    }
    
    # unlist and binds all together
    outlist <- bind_rows(outlist) %>% distinct() %>% arrange(bin)
    
    setProgress(2, detail = paste("building list"))
    # removes top and bottom %
    if (sum(start1_bin, end1_bin) > sum(start2_bin, end2_bin)) {
      use_header <- "Log2 EI Cumulative plot"
    } else {
      use_header <- "Log2 PI Cumulative plot"
    }
    if (n_distinct(outlist$gene) > 0) {
      nick_name1 <- paste0("CDF ", use_header)
      list_data$gene_file[[nick_name1]]$full <- outlist
      list_data$gene_file[[nick_name1]]$use <- outlist %>% select(gene)
      list_data$gene_file[[nick_name1]]$info <- tibble(loaded_info =
        paste(
          use_header,
          "CDF",
          "bins",
          start1_bin,
          "to",
          end1_bin,
          "/",
          start2_bin,
          "to",
          end2_bin,
          "from",
          names(onoffs),
          "gene list(s)",
          paste(distinct(outlist, plot_set), collapse = " "),
          Sys.Date()
        ),
        matching = FALSE)
    } else {
      nick_name1 <- paste("CDF n = 0")
    }
    for (list_name in names(onoffs)) {
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == list_name &
                                             set %in% onoffs[[list_name]]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste(
                                             "CDF from",
                                             list_name
                                           ), 
                                           onoff = "0",
                                           count = paste("n =", outlist %>% 
                                                           summarise(n=n_distinct(bin))),
                                           plot_set = paste(list_name, "-", set),
                                           myheader = use_header)))
      setProgress(5, detail = "finishing up")
    }
    list_data
  }
# CDG ggplot function
GGplotC <-
  function(df2,
           plot_options,
           use_header,
           plot_range) {
    print("ggplot CDF")
    use_col <- plot_options$mycol
    names(use_col) <- plot_options$set
    legend_space <- max(1, (lengths(strsplit(
      plot_options$set, "\n"
    ))))
    gp <- ggplot(df2, aes(value, color = set)) +
      stat_ecdf(show.legend = TRUE, size = 1.8) +
      scale_color_manual(name = "Sample", values = use_col) +
      ylab("Fraction of genes") +
      ggtitle(use_header) +
      theme_bw() +
      theme(legend.title = element_blank()) +
      theme(axis.title.y = element_text(size =  15)) +
      theme(axis.title.x = element_text(size =  13, vjust = .5)) +
      theme(axis.text.x = element_text(
        size = 12,
        angle = -45,
        hjust = .1,
        vjust = .9,
        face = 'bold'
      )) +
      theme(
        legend.title = element_blank(),
        legend.key = element_rect(size = 5, color = 'white'),
        legend.key.height = unit(legend_space, "line"),
        legend.text = element_text(size = 10)
      ) +
      coord_cartesian(xlim = plot_range)
    suppressMessages(print(gp))
    return(suppressMessages(gp))
  }
