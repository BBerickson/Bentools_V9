
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
  # print("gene match fun")
  if(str_detect(gene_list$gene[1],"|")){
    tablefile <-
      map(paste0(";", gene_list$gene,"\\|"), str_subset, string = common_list$gene) %>% 
      setNames(gene_list$gene)
  } else if(str_detect(gene_list$gene[1],"_")){
    tablefile <-
      map(paste0(";", gene_list$gene,"\\|"), str_subset, string = common_list$gene) %>% 
      setNames(gene_list$gene)
  }else {
    tablefile <-
      map(paste0("\\|", gene_list$gene,"$"), str_subset, string = common_list$gene) %>% 
      setNames(gene_list$gene)
  }
  tablefile <- map_df(tablefile, ~as.data.frame(.x), .id="gene") %>% 
    distinct(gene,.keep_all = T) %>% 
    full_join(.,gene_list,by="gene") %>% 
    transmute(org_gene = gene, gene=.x)
  return(tablefile)
}

# lines and labels preset helper for older table files
LinesLabelsPreSetGuess <- function(mytype) {
  # type,binsize,upstream,downstream,body,unscaled5prime,unscaled3prime
  if (mytype == "543") {
    tt <- c(543, 100, 1500, 3500, 2000, 500, 500)
  } else if (mytype == "5" | mytype == "TSS") {
    tt <- c(5, 25, 1000, 1000, 0,0,0)
  } else if (mytype == "5L"| mytype == "TSS") {
    tt <- c(5, 100, 100, 30000, 0,0,0)
  } else if (mytype == "3" | mytype == "TES") {
    tt <- c(3, 100, 1000, 9000, 0,0,0)
  } else if (mytype == "PI") {
    tt <- c(543, 400, 400, 400, 400, 0, 0)
  } else {
    tt <- c(0, 100, 1500, 3500, 2000, 500, 500)
  }
  tt
}

# takes info from file type and number of bins to pre set tools sliders
SlidersSetsInfo <- function(slider_breaks, type){
  # 5Min, 5Max, 3Min, 3Max
  num_bins <- max(slider_breaks$mybrakes)
  if (num_bins == 80 & type == '543') { 
    setsliders <- slider_breaks$mylabels[c(15,19,20,47)] 
  } else if (num_bins == 80 & (type == '5'| type == "TSS")) {
    setsliders <- c(slider_breaks$mylabels[c(21,61)],NA,NA)
  } else if (num_bins == 2 | type == 'PI') {
    setsliders <- slider_breaks$mylabels[c(1,1,2,2)]
  } else if (num_bins == 205 & (type == '5L'| type == "TSS")) {
    setsliders <- slider_breaks$mylabels[c(15,19,20,47)]
  } else if (type == '3'| type == "TES") {
    setsliders <- c(num_bins/2,
                    num_bins,num_bins/2+1,num_bins)
  } else {
    setsliders <- c(
      floor(num_bins / 5.5),
      floor(num_bins / 4.4),
      floor(num_bins / 4.4) + 1,
      floor(num_bins / 1.77))
  }
  setsliders
}

# sort out info on file/url
PrepMetaFile <-
  function(file_path,
           file_name) {
    # print("PrepMetaFile")
    # tests if loading a file with a list of address to remote files, requires .url.txt in file name
    if (str_detect(file_name, ".url.txt")) {
      num_col <-
        try(count_fields(
          file_path,
          n_max = 1,
          skip = 1,
          tokenizer = tokenizer_delim(" ")
        ),
        silent = T)
      if ("try-error" %in% class(num_col)) {
        showModal(modalDialog(
          title = "Information message",
          paste(file_name, "cant find file to load"),
          size = "s",
          easyClose = TRUE
        ))
        return(NULL)
      }
      if (num_col > 1) {
        col_names <- c("filepath", "type", "nick", "color")
        col_types <- cols(
          filepath = col_character(),
          type = col_character(),
          nick = col_character(),
          color = col_character()
        )
      } else {
        col_names <- c("filepath")
        col_types <- cols(filepath = col_character())
      }
      meta_data <-
        suppressMessages(read_delim(
          file_path,
          delim = " ",
          col_names = col_names,
          col_types = col_types
        )) 
      if (num_col == 1) {
        meta_data <-
          meta_data %>% mutate(
            nick = str_split(filepath, "/", simplify = T) %>%
              as_tibble(., name_repair =
                          "unique") %>%
              select(last_col()) %>% unlist() %>%
              str_remove(., ".table|.table.gz|.matrix.gz") %>% str_replace("\\.", "_"),
            type = str_extract(nick, "^(\\d)+") %>%
              replace_na("none"),
            color = sample(
              suppressWarnings(brewer.pal(11, sample(
                kBrewerList, size = 1
              ))) %>%
                grep("#FF", ., value = T, invert = T),
              size = nrow(meta_data)
            )
          )
      } else {
        meta_data <-
          meta_data %>%
          mutate(nick = str_replace(nick, "\\.", "_"))
      }
    } else {
      meta_data <-  tibble(
        filepath = file_path,
        nick = last(str_split(file_name, "/", simplify = T)) %>%
          str_remove(., ".table|.matrix.gz") %>% str_replace("\\.", "_"),
        type = str_extract(nick, "^(\\d)+") %>%
          replace_na("none"),
        color = sample(
          suppressWarnings(brewer.pal(11, sample(
            kBrewerList, size = 1
          ))) %>%
            grep("#FF", ., value = T, invert = T),
          size = 1
        )
      )
    }
    meta_data
  }

# test file is there and what type
tableTestbin <- function(meta_data){
  # get info on file to help know what type it is
  # print("tableTestbin")
  num_bins <-
    try(count_fields(meta_data$filepath,
                     n_max = 1,
                     tokenizer = tokenizer_tsv()),silent = T)
  binning <- LinesLabelsPreSetGuess(meta_data$type)
  # test if file can be loaded in
  if ("try-error" %in% class(num_bins)) {
    showModal(modalDialog(
      title = "Information message",
      paste(meta_data$nick, "cant find file to load"),
      size = "s",
      easyClose = TRUE
    ))
    return()
  }
  # check if table file has meta data
  if (str_detect(read_lines(meta_data$filepath,n_max = 1),"# meta=")){
    binning <- str_remove(trimws(read_lines(meta_data$filepath,n_max = 1)),"# meta=") %>% 
      str_split_fixed(.,",",n=7) %>% as.numeric()
  }
  # check if file is in wide format or deeptools matrix file
  if (num_bins == 1 | str_detect(meta_data$filepath, "matrix.gz$")) {
    num_bins <-
      count_fields(meta_data$filepath,
                   n_max = 1,
                   skip = 1,
                   tokenizer = tokenizer_tsv()) - 5
    col_names <- c("chr", "start", "end","gene", "value", "sign", 1:(num_bins - 1))
    mylist <- c("bin size","upstream","downstream","body","unscaled 5 prime","unscaled 3 prime")
    
    meta <- read_tsv(meta_data$filepath,n_max = 1,col_names = F)
    mm <- meta %>% str_remove_all("[@{}]|\\]|\\[") %>% str_split(",",simplify = T) %>% 
      str_replace_all(.,fixed('\"'),"") 
    type <- mm[str_which(mm,"ref point")] %>% str_replace_all(., "ref point:", "")
    if(type == "TSS"){
      type <- 5
    } else if(type == "TES") {
      type <- 3
    } else{
      type <- 543
    }
    header <- type
    for(i in mylist){
      header <- c(header,mm[str_which(mm,i)] %>% 
                        str_remove(.,i) %>% 
                        str_replace_all(., "[a-zA-Z: ]", ""))
    }
    binning <- header %>% as.numeric()
  } else if (num_bins > 6 ){
    col_names <- c("gene", 1:(num_bins - 1))
  } else if (num_bins == 4) {
    # settings for new style with meta data info
    col_names <- c("gene", "bin", "score", "set")
    # settings for reading in bedGraph file style
  } else if (num_bins == 3) {
    col_names <- c("gene", "bin", "score")
  } else {
    col_names <- NULL
  }
  list(num_bins = num_bins, col_names = col_names, binning = binning)
}

LoadTableFile <-
  function(meta_data,
           bin_colname) {
    # print("LoadTableFile")
    # wide file
    if (bin_colname$num_bin > 6) {
      tablefile <- suppressMessages(
        read_tsv(
          meta_data$filepath,
          comment = "#",
          col_names = bin_colname$col_names,
          skip = 1
        )
      ) %>%
        select(gene, num_range(
          range = 1:(bin_colname$num_bins - 1),
          prefix = ""
        )) %>%
        gather(., bin, score, 2:(bin_colname$num_bins))
      # long file
    } else {
      tablefile <-
        suppressMessages(read_tsv(
          meta_data$filepath,
          comment = "#",
          col_names = bin_colname$col_names
        ))
    }
    tablefile %>%
      dplyr::mutate(bin = as.numeric(bin),
                    score = as.numeric(score),
                    set = meta_data$nick) %>%
      na_if(Inf) %>%
      replace_na(list(score = 0)) %>%
      distinct(gene, bin, .keep_all = T)
  }

# reads in gene list file(s), tests, fills out info and returns list_data
LoadGeneFile <-
  function(file_path,
           file_name,
           list_data) {
    legend_nickname <- last(str_split(file_name,"/",simplify = T)) %>% 
        str_remove(., ".txt") %>% str_replace("\\.","_")
      # checks if file with same nickname has already been loaded
      if (legend_nickname %in% names(list_data$gene_file)) {
        showModal(modalDialog(
          title = "Information message",
          paste(file_name, "has already been loaded"),
          size = "s",
          easyClose = TRUE
        ))
        return()
      }
      # gets number of columns in file, used to guess how to deal with file
      #  and checks if file exits
      num_col <-
        try(count_fields(file_path,
                         n_max = 1,
                         skip = 1,
                         tokenizer = tokenizer_tsv()),silent = T)
      if ("try-error" %in% class(num_col) & num_col > 0) {
        showModal(modalDialog(
          title = "Information message",
          paste(file_name, "cant find file to load"),
          size = "s",
          easyClose = TRUE
        ))
        return()
      }
      if(num_col == 1){
        #normal gene list
        col_names <- c("gene")
      } else if(str_detect(file_path, ".bed")){
        col_names <- 1:num_col
        col_names[4] <- "gene"
      } else {
        col_names <- FALSE
      }
      # reads in file
      tablefile <-
        suppressMessages(read_tsv(file_path,
                                  comment = "#",
                                  col_names = col_names)) 
      if(!"gene" %in% names(tablefile)){
        names(tablefile)[1] <- "gene"
      }
      tablefile <- tablefile %>% 
        distinct(gene,.keep_all = T) %>% mutate(gene=as.character(gene))
      # checks gene list is a subset of what has been loaded
      gene_names <- tablefile %>% select(gene) %>% 
        semi_join(., list_data$gene_file$Complete$full, by = "gene") %>% distinct(gene)
      # test data is compatible with already loaded data
      if (n_distinct(gene_names$gene, na.rm = T) == 0) {
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
        gene_names <- MatchGenes(list_data$gene_file$Complete$full, tablefile %>% select(gene))
        if (n_distinct(gene_names$gene, na.rm = T) == 0) {
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
          showModal(
            modalDialog(
              title = "Information message",
              paste("Found", n_distinct(gene_names$gene, na.rm = T), "genes after pattern matching search"),
              size = "s",
              easyClose = TRUE
            )
          )
        }
      }
      # adds full n count to nickname
      my_name <- legend_nickname
      # preps meta data
      gene_info <- list_data$gene_info %>% 
        dplyr::filter(gene_list == "Complete") %>% 
        dplyr::mutate(gene_list = my_name, 
                      count = paste("n =", n_distinct(gene_names$gene, na.rm = T)),
                      sub = " ", onoff = "0",
                      plot_set = " ")
      # saves data in list of lists
      list_data$gene_file[[my_name]]$full <- distinct(gene_names)
      list_data$gene_file[[my_name]]$info <- tibble(loaded_info =
                                                      paste("Loaded gene list from file",
                                                            legend_nickname,
                                                            Sys.Date()),
                                                    save_name = gsub(" ", "_", str_squish(paste(legend_nickname,
                                                                      Sys.Date(), sep ="_"))),
                                                    col_info = "loaded file"
                                                    )
      list_data$gene_info <- bind_rows(list_data$gene_info, gene_info)
      
    return(list_data)
  }

# input data list, output filtered and bound data with plot legend column
Active_list_data <-
  function(list_data, group="none") {
    table_file <- list_data$table_file
    gene_file <- list_data$gene_file
    gene_info <- list_data$gene_info
    list_data_out <- NULL
    # print("active data function")
    for ( i in names(gene_file) ){
      # checks to see if at least one file in list is active
      if (gene_info %>% dplyr::filter(gene_list == i & onoff != 0) %>% nrow() == 0) {
        next()
      } else {
        if(group != "none"){
          gene_info <- gene_info %>% group_by(group,gene_list) %>% 
            mutate(onoff=if_else(gene_list == i & any(onoff != 0),set,onoff)) %>% 
            ungroup()
        }
        my_sel <- gene_info %>% dplyr::filter(gene_list == i & onoff != 0)
        tf <- table_file %>% 
          dplyr::filter(set %in% my_sel$onoff)
        gene_common <- tf %>% group_by(set) %>% distinct(gene) %>% ungroup()
        gene_common <- gene_common %>% 
          group_by(gene) %>% filter(n_distinct(set)==length(my_sel$set)) %>% 
          ungroup() %>% distinct(gene)
        list_data_out[[i]] <- tf %>% 
          semi_join(., gene_common, by = "gene") %>%
          semi_join(., gene_file[[i]]$full, by = "gene") %>% 
          dplyr::mutate(., gene_list = i)
        # test for empty results
        if(is_empty(list_data_out[[i]]$gene)){
          list_data_out[[i]] <- tf %>% 
            semi_join(., gene_file[[i]]$full, by = "gene") %>% 
            complete(.,my_sel %>% select(set),bin,gene) %>% 
            replace_na(list(score = 0)) %>% 
            dplyr::mutate(., gene_list = i)
        }
        # adds line brake at 20 character for legend spacing
        my_sel2 <- my_sel %>% dplyr::mutate(.,plot_set = paste(
          gsub("(.{20})", "\\1\n", set),
          gsub("(.{20})", "\\1\n", 
               str_split_fixed(i, "\nn = ", n=2)[,1]),
          paste0("n = ", n_distinct(list_data_out[[i]]$gene, na.rm = T)),
          sep = '\n'
        ),group=paste(
          gsub("(.{20})", "\\1\n", group)
        )) %>% dplyr::select(set,plot_set,group)
        list_data_out[[i]] <- list_data_out[[i]] %>% inner_join(.,my_sel2,by="set")
      }
    }
    return(bind_rows(list_data_out))
  }

# apply math to data file and get ready to plot
ApplyMath <-
  function(list_data,
           use_math,
           relative_frequency,
           normbin,
           group="none") {
    # print("applymath fun")
    # normalize per gene relative frequency
    if(is.empty(normbin)){
      normbin <- 0
    }
    if (relative_frequency == "rel gene frequency") {
      list_data <- list_data %>% group_by(plot_set, gene) %>%
        dplyr::mutate(score = score / abs(sum(score, na.rm = TRUE))) %>%
        ungroup()
    }
    # apply mean/median/sum/var
    list_data <- list_data %>% group_by(set, plot_set, bin, group, gene_list) %>%
      summarise(value = get(use_math)(score, na.rm = T), .groups="drop")
    # norm to bin or overall relative frequency
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
    # finish making file ready for ggplot
    if(group == "groups only"){
      list_data <- select(list_data, -plot_set) %>% 
        right_join(.,distinct(list_data,group,gene_list,.keep_all = T) %>% 
                                select(group, plot_set, gene_list),by=c("group","gene_list")) %>% 
        separate(plot_set,c("cc","set2"),"\n",extra = "merge",remove = F) %>% dplyr::select(-cc) %>%
        dplyr::mutate(group=if_else(set != group, paste(group,set2,sep = "\n"),plot_set)) %>%
        transmute(set = plot_set,plot_set=group,bin=bin,value=value) %>% 
         group_by(set, bin,plot_set) %>% 
        summarise(min=min(value),max=max(value),
                  value = mean(value), .groups = "drop") 
      } else if (group == "groups and single"){
        list_data2 <- select(list_data, -plot_set) %>% 
          right_join(.,distinct(list_data,group,gene_list,.keep_all = T) %>% 
                       select(group, plot_set, gene_list),by=c("group","gene_list")) %>% 
          separate(plot_set,c("cc","set2"),"\n",extra = "merge",remove = F) %>% dplyr::select(-cc) %>%
          dplyr::mutate(group=if_else(set != group, paste(group,set2,sep = "\n"),plot_set)) %>%
          transmute(set = plot_set,plot_set=group,bin=bin,value=value) %>% 
          group_by(set, bin,plot_set) %>% 
          summarise(min=min(value),max=max(value),
                    value = mean(value), .groups = "drop") 
        list_data <- list_data %>% 
          mutate(set=plot_set,min=value,max=value) %>% 
            bind_rows(.,list_data2)
      } else {
        list_data <- list_data %>% 
          mutate(set=plot_set,min=value,max=value)
    }
    return(list_data)
  }

# get min and max from apply math data set
YAxisValues <-
  function(apply_math,
           xBinRange,
           log_2 = F,
           yBinRange = c(0, 100)) {
    tt <- group_by(apply_math, set) %>%
      dplyr::filter(bin %in% xBinRange[1]:xBinRange[2]) %>%
      ungroup() %>%
      summarise(min(min, na.rm = T), max(max, na.rm = T),.groups="drop") %>%
      unlist(., use.names = FALSE)
    tt <-
      c(tt[1] + (tt[1] * (yBinRange[1] / 100)), tt[2] + (tt[2] * ((yBinRange[2] -
                                                                     100) / 100)))
    if (log_2) {
      tt <- log2((abs(tt))^(sign(tt)))
    }
    tt <- c(round(min(tt), 4),round(max(tt), 4))
    tt
  }

# Sets y label fix
YAxisLabel <-
  function(use_math = "mean",
           relative_frequency = "none",
           norm_bin = "NA",
           smoothed = F,
           log_2 = F) {
    use_y_label <- paste(use_math, "of bin counts")
    if (relative_frequency == "rel gene frequency") {
      use_y_label <- paste("RF per gene :", use_y_label)
    } else if (relative_frequency == "relative frequency") {
      use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                           "bins : RF")
    }
    if (norm_bin != "NA") {
      if (relative_frequency == "rel gene frequency") {
        use_y_label <- paste(use_y_label, " : Norm ", norm_bin)
      } else {
        use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                             "bins : Normalize to ",
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
  # print("plot options fun")
  # checks to see if at least one file in list is active
  if (gene_info %>% dplyr::filter(onoff != 0) %>% nrow() == 0) {
    return(NULL)
  } else {
    gene_info <- gene_info %>%
      dplyr::mutate(
        myline = 1,
        set=plot_set
      )
  }
  # tint if same color is used more then once
  ldf <- duplicated(gene_info[c("mycol","onoff")]) & gene_info$onoff != 0
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
           plot_occupancy,
           auc = FALSE) {
    plot_options <- list_long_data_frame %>% 
      distinct(set,plot_set) %>% right_join(plot_options,.,by="set") %>% 
      dplyr::rename(plot_set=plot_set.y) %>% dplyr::select(-plot_set.x) %>% 
      dplyr::mutate(set = plot_set) %>% distinct(.,set,.keep_all = T)
    list_long_data_frame <- list_long_data_frame %>% 
      dplyr::mutate(set = plot_set)
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
            linetype = set
          )
        )
    }
    if (use_smooth) {
      gp <- gp +
        geom_smooth(se = FALSE,
                    size = line_list$mysize[2],
                    span = .2, alpha=line_list$mysize[6]) 
    } else{
      gp <- gp +
        geom_line(size = line_list$mysize[2],alpha=line_list$mysize[6])
    }
    gp <- gp + 
      geom_ribbon(aes(ymin=min,ymax=max,fill=set),linetype=0,alpha = 0.5) +
      scale_color_manual(values = plot_options$mycol) +
      scale_fill_manual(values = plot_options$mycol) +
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
    if(auc){
      plot_options <- list_long_data_frame %>% 
        group_by(set) %>% 
        summarise(AUC=paste("AUC =",round(AUC(value),2))) %>% 
        inner_join(.,plot_options,by="set") 
      gp <- gp + 
        annotate("text",
                 x=floor(seq(xBinRange[1],xBinRange[2],length=nrow(plot_options))),
                 y= Inf,vjust = 1.5, hjust = "inward",
                 color = plot_options$mycol,
                 label = plot_options$AUC)
      
    }
    if(!is_empty(LIST_DATA$ttest)){
      use_col_tt <- plot_ttest$options_main_tt$mycol
      use_line_tt <- plot_ttest$options_main_tt$myline
      names(use_col_tt) <- plot_ttest$options_main_tt$set
      names(use_line_tt) <- plot_ttest$options_main_tt$set
      gp2 <- ggplot(LIST_DATA$ttest, aes(y=p.value,x=bin,
                                         color=set,
                                         linetype = set)) + 
        geom_line(size = line_list$mysize[2],alpha=line_list$mysize[6]) +
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

# gets y axis label landmarks
LinesLableLandmarks <- function(myinfo){
  # type, bp/bin, before, after, body, un5, un3, spacing
  # myinfo <- c(543,100,1500,3500,2000,500,500,500)
  tssbin <- myinfo[3]/myinfo[2]
  if (sum(myinfo[5:7]) > 0) {
    body1bin <- tssbin + myinfo[6]/myinfo[2]
    body2bin <- body1bin + myinfo[5]/myinfo[2]
  } else{
    body1bin <- 0
    body2bin <- 0
  }
  tesbin <- sum(myinfo[c(3,5:7)])/myinfo[2]
  floor(c(tssbin, tesbin, body1bin, body2bin, myinfo[8]/myinfo[2]))
}

# Sets lines and labels
LinesLabelsSet <- function(myinfo,
    totbins = 80,
    tssname = "TSS",
    tesname = "pA",
    slider = F) {
  # LinesLabelsSet(c(543,100,1500,3500,2000,500,500,500),slider = F)
  if (myinfo[8] > 0) {
    if(totbins > 2){
    if(slider){
      myinfo[8] <- myinfo[2]
    }
    landmarks <- LinesLableLandmarks(myinfo)
    # y asis locations and labels for 1:before
    if(landmarks[1] > 0){
      mod <- 0.5
    } else {
      mod <- 1
    }
    before <- seq(-myinfo[3], 0, by = myinfo[8])
    beforebins <- seq(1,  by = landmarks[5], length.out = length(before))
    # landmark1 is 5' end
    if(str_detect(myinfo[1], "^3|TES")){
      tssname <- tesname
      before <- abs(before)
    }
    # make sure landmark is included properly 
    myloc <- which(before == 0 | beforebins == landmarks[1])
    if(!is_empty(myloc)){
      myloc <- max(myloc)
    }
    if (any(before == 0) | any(beforebins == landmarks[1])) {
      beforebins[myloc] <- landmarks[1] + mod
      before[myloc] <- tssname
    } else {
      beforebins <- sort(c(beforebins, landmarks[1] + mod))
      before <- append(before, tssname)
    }
    
    # any other landmarks
    if(all(c(landmarks[3],landmarks[4]) > 0)){
      # test for unscaled5prime and unscaled3prime
      if(landmarks[3] > landmarks[1] & landmarks[4] > landmarks[3] & all(landmarks%%landmarks[5] == 0)){
        # landmark1 to unscaled5prime
        unscaled5prime <- seq(myinfo[8], (landmarks[3] - landmarks[1]) * myinfo[2], by = myinfo[8])
        unscaled5primebin <- seq(landmarks[1]+landmarks[5], by = landmarks[5], length.out = length(unscaled5prime))
        if(landmarks[1] == 0 & myinfo[8] == 1){
          unscaled5prime <- unscaled5prime + 1
          unscaled5primebin <- unscaled5primebin + 1
        }
        # make sure body brake is included
        if (!any(unscaled5primebin == landmarks[3])) {
          unscaled5prime <- append(unscaled5prime, (landmarks[3] - landmarks[1]) * myinfo[2])
          unscaled5primebin <- c(unscaled5primebin, landmarks[3])
        }
        # slider body
        if(slider){
          unscaled5prime <- c(paste0("5'unscale_",unscaled5prime), 
                              paste0("scaled_",seq_along((landmarks[3]+1):(landmarks[4]))))
          unscaled5primebin <- c(unscaled5primebin, (landmarks[3]+1):(landmarks[4]))
        }
        # unscaled3prime to last landmark
        unscaled3prime <-  abs(seq((landmarks[4] - landmarks[2]) * myinfo[2], 0, by = myinfo[8]))
        if(landmarks[4] == 1){
          unscaled3prime <- unscaled3prime + myinfo[2]
        }
        unscaled3primebin <- seq(landmarks[4]+1, by = landmarks[5], length.out = length(unscaled3prime))
        # make sure TES is included
        myloc <- which(unscaled3prime == 0 | unscaled3primebin == landmarks[2]+1)
        if(!is_empty(myloc)){
          myloc <- max(myloc)
        }
        if (any(unscaled3prime == 0) | any(unscaled3primebin == landmarks[2]+1)) {
          unscaled3primebin[myloc] <- landmarks[2] + .5
          unscaled3prime[myloc] <- tesname
        } else {
          unscaled3primebin <- sort(c(unscaled3primebin, landmarks[2] + .5))
          unscaled3prime[which(unscaled3primebin == landmarks[2] + .5)] <- tesname
        }
        if(slider){
          unscaled3prime[which(unscaled3prime != tesname)] <- 
            paste0("3'unscale_",unscaled3prime[which(unscaled3prime != tesname)])
        }
        before <- c(before,unscaled5prime,unscaled3prime)
        beforebins <- c(beforebins,unscaled5primebin,unscaled3primebin)
      } 
      tt <- c(myinfo[8], abs(totbins - landmarks[2]) * myinfo[2])
      TESname <- seq(min(tt),max(tt), by = myinfo[8])
      TESloc <-
        seq(landmarks[2] + landmarks[5],
            by = landmarks[5],
            length.out = length(TESname))
      # make sure last location is included
      if (!any(TESloc == totbins)) {
        TESname <- append(TESname, abs(totbins - landmarks[2]) * myinfo[2])
        TESloc <- c(TESloc, totbins)
      }
      # make sure TES is included
      if (!any(before == tesname)){
        before <- c(before,tesname,TESname)
        beforebins <- c(beforebins,landmarks[2] + .5,TESloc)
      } else {
        before <- c(before,TESname)
        beforebins <- c(beforebins,TESloc)
      }
      # just 5' or 3'
    } else if(myinfo[8] <= myinfo[4]){
      landmark  <- trunc(last(beforebins)) + landmarks[5] 
      TESname <- seq(myinfo[8], myinfo[4], by = myinfo[8])
      TESloc <-
        seq(landmark,
            by = landmarks[5],
            length.out = length(TESname))
      # make sure last location is included
      if (!any(TESloc == totbins)) {
        TESname <- append(TESname, myinfo[4])
        TESloc <- c(TESloc, totbins)
      }
      before <- c(before,TESname)
      beforebins <- c(beforebins,TESloc)
    } else {
      before <- c(before,NA)
      beforebins <- c(beforebins,totbins)
    }
    # put it all together
    use_plot_breaks <- beforebins
      use_plot_breaks_labels <- before
    } else{
      use_plot_breaks <- c(1,2)
      use_plot_breaks_labels <- c(tssname,tesname)
    }
  } else {
    # just print bin numbers
    use_plot_breaks <-
      seq(1,totbins,by=ceiling(totbins/10))
    use_plot_breaks_labels <- use_plot_breaks
  }
  list(mybrakes = use_plot_breaks,
       mylabels = use_plot_breaks_labels)
}

# Sets plot lines and labels colors
# makes plot ascetics  
LinesLabelsPlot <-
  function(myinfo,
           body1color,
           body1line,
           body2color,
           body2line,
           tsscolor,
           tssline,
           tescolor,
           tesline,
           use_plot_breaks_labels,
           use_plot_breaks,
           vlinesize,
           linesize,
           fontsizex,
           fontsizey,
           legendsize,
           myalpha) {
    # print("lines and labels plot fun")
    # myinfo <- c(543,100,1500,3500,2000,500,500,500)
    landmarks <- LinesLableLandmarks(myinfo)
    tssbin <- landmarks[1]
    tesbin <- landmarks[2]
    body1bin <- landmarks[3]
    body2bin <- landmarks[4]
    binspace <- landmarks[5]
    binsize <- myinfo[8]
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
      mysize = c(vlinesize, linesize, fontsizex, fontsizey, legendsize, myalpha),
      myset = c(body1bin, body2bin, tssbin, tesbin, binsize, binspace)
    )
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
           padjust,
           my_alt, 
           my_exact, 
           my_paired,
           group = "none") {
    # t.test comparing files in same gene list
    ttest <- NULL
    if(switchttest != "none"){
      if(group == "groups only"){
        list_data <- list_data %>% 
          dplyr::mutate(set= gsub("\n", "", group)) 
      }
      if(switchttest == "by lists" & n_distinct(list_data$gene_list, na.rm = T) > 1){
        list_data <- list_data %>% 
          rename(set=gene_list,gene_list=set)
      }
      n_test <- list_data %>% group_by(gene_list) %>% 
      summarise(n=n_distinct(set), .groups= "drop")
    
      for(i in n_test$gene_list){
        if(n_test %>% dplyr::filter(gene_list == i) %>% dplyr::select(n) > 1){
          ttest[[i]] <- bind_rows(try_t_test(list_data %>% dplyr::filter(gene_list == i), i,
                                             use_tmath,switchttesttype,padjust,
                                             my_alt, noquote(my_exact), noquote(my_paired))) %>% 
            dplyr::mutate(myline = 1,
                          mycol = "#000000" )
        } 
      }
    } 
    bind_rows(ttest)
  }

# makes sure t test wont crash on an error
try_t_test <- function(db,my_set,my_math ="none",my_test="t.test",padjust="fdr",
                       alternative="two.sided", exact=FALSE, paired=FALSE){
  # print("try t.test")
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
    db_out[[str_c(i,collapse = "-")]] <- myTtest %>% 
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
  # print("plot options ttest fun")
  out_options <- list_data %>% select(-bin, -p.value) %>% distinct(set,myline,mycol)
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
           start_end_bin,
           start_end_label,
           mynum,
           topbottom) {
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
      apply_bins <-
        semi_join(dplyr::filter(list_data$table_file, set == j), 
                  list_data$gene_file[[list_name]]$full, by = 'gene')  
      apply_bins <- group_by(apply_bins, gene) %>%
        dplyr::filter(bin %in% min(start_end_bin):max(start_end_bin)) %>%
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
        outlist <<- inner_join(outlist, outlist2, by = 'gene')
      } else {
        outlist <<- outlist2
      }
      lc <<- lc + 1
    })
    if (length(outlist$gene) == 0) {
      return(NULL)
    }
    old_names <- grep("^Filter", names(list_data$gene_file), value = T)
    if (length(old_names) > 3) {
      old_names <- old_names[!old_names %in% list_name]
      # remove old sort gene list keeping max 4
      list_data$gene_file[[first(old_names)]] <- NULL
      list_data$gene_info <- dplyr::filter(list_data$gene_info,
                                           gene_list != first(old_names))
    }
    topbottom2 <- paste(str_remove(topbottom,"%"), paste0(mynum, "%"))
    nick_name <-
      strtrim(gsub("(.{30})",
                   "\\1... ",
                   paste0("Filter ",topbottom2, "\nn = ", n_distinct(outlist$gene, na.rm = T))), 33)
    list_data$gene_file[[nick_name]]$full <- outlist
    list_data$gene_file[[nick_name]]$info <- tibble(loaded_info =
      paste(
        "Filter",
        topbottom2,
        start_end_label[1],
        "to",
        start_end_label[2],
        "from",
        list_name,
        paste(file_names, collapse = " "),
        Sys.Date(),
        list_data$gene_file[[list_name]]$info
      ),
      save_name = gsub(" ", "_", paste("Filter",str_remove(topbottom,"%"), Sys.Date(), sep = "_")),
      col_info = "gene [ % rank(s) ]"
      )
    list_data$gene_info <- 
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>% 
                           dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                           dplyr::mutate(gene_list = nick_name,
                                         sub =  paste("Filter",
                                                      topbottom2,
                                                      start_end_label[1],
                                                      "to",
                                                      start_end_label[2]), 
                                         onoff = "0",
                                         count = paste0("n = ", n_distinct(outlist$gene, na.rm = T)),
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
           start_end_label) {
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
    gene_list <- list_data$gene_file[[list_name]]$full
    out_list <- list_data$table_file %>% 
      dplyr::filter(set %in% file_names) %>% 
      semi_join(.,gene_list,by="gene") %>% 
      dplyr::filter(bin %in% start_end_bin[1]:start_end_bin[2]) 
    
    out_per <- out_list %>%
      group_by(set,bin) %>% summarize_at(vars(score), p_funs) %>% ungroup()
    # if all == 0 finds the first point signal starts
    my_step <- .1
    while(sum(out_per$my_p_1,na.rm = T) == 0 & my_per[1] < 100){
      my_per[1] <- my_per[1] + my_step
      if(my_per[1] > 100){
        my_per[1] <- 99.9
      }
      p_funs <- map(my_per[1]/100, ~partial(quantile, probs = .x, na.rm = TRUE)) %>%
        set_names(paste0("my_p_",seq_along(my_per[1])))
      out_per <- out_list %>%
        group_by(set,bin) %>% summarize_at(vars(score), p_funs) %>% ungroup() %>% 
        rows_update(out_per,.,by=c("set","bin"))
      if(my_step < 1){
        my_step <- my_step + .1
      } else if(between(my_step,1,10)){
        my_step <- my_step + 2
      }else if(between(my_step,10,50)){
        my_step <- my_step + 10
      }else if(between(my_step,50,.100)){
        my_step <- my_step + 5
      } else {
        my_step <- my_step + 1
      }
      if(my_per[1] >= 99.9){
        my_per[1] <- 100
      }
    }
    if(my_type == "min%"){
      out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene) %>% dplyr::filter(all(score >= my_p_1)) %>%
        filter(n_distinct(set)==length(file_names)) %>%
          ungroup() %>% distinct(gene)
      topbottom2 <- paste(str_remove(my_type,"%"), paste0(my_per[1], "%"))
    } else if(my_type == "max%"){
      out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene) %>% dplyr::filter(all(score <= my_p_2)) %>%
        filter(n_distinct(set)==length(file_names)) %>%
          ungroup() %>% distinct(gene)
      topbottom2 <- paste(str_remove(my_type,"%"), paste0(my_per[2], "%"))
    } else {
      out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene) %>% dplyr::filter(all(score >= my_p_1 & score <= my_p_2)) %>% 
        filter(n_distinct(set)==length(file_names)) %>%
          ungroup() %>% distinct(gene)
      topbottom2 <- paste(paste(str_remove(my_type,"%"), paste0(my_per[1], "%")),paste0(my_per[2], "%"),collapse = " and ")
    }
    if (length(out_list$gene) == 0) {
      return(NULL)
    }
    old_names <- grep("^Filter", names(list_data$gene_file), value = T)
    if (length(old_names) > 3) {
      old_names <- old_names[!old_names %in% list_name]
      # remove old sort gene list keeping 4
      list_data$gene_file[[first(old_names)]] <- NULL
      list_data$gene_info <- dplyr::filter(list_data$gene_info,
                                           gene_list != first(old_names))
    }
    nick_name <-
      strtrim(gsub("(.{30})",
                   "\\1... ",
                   paste0("Filter Prob ",topbottom2, "\nn = ", n_distinct(out_list$gene, na.rm = T))), 33)
    if(length(file_names) == 1){
      my_per2 <- seq(my_per[1],my_per[2]/2.1,length.out = 5)
      p_funs <- map(my_per2/100, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
        set_names(my_per2)
      my_list <- list_data$table_file %>% 
        dplyr::filter(set == file_names) %>% 
        semi_join(.,gene_list,by="gene") %>% 
        dplyr::filter(bin %in% start_end_bin[1]:start_end_bin[2]) 
      out_per1 <- my_list %>%
        group_by(set,bin) %>% summarize_at(vars(score), p_funs) %>% ungroup() %>% 
        dplyr::mutate(across(where(is.double),as.character)) %>% 
        gather(.,key = set,value = "my_p_1",-bin,-set) %>% 
        dplyr::mutate(my_p_1=as.double(my_p_1))
      
      my_per2 <- seq(my_per[2]/2,my_per[2],length.out = 5)
      p_funs <- map(my_per2/100, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
        set_names(my_per2)
      
      out_per2 <- my_list %>%
        group_by(set,bin) %>% summarize_at(vars(score), p_funs) %>% ungroup() %>% 
        dplyr::mutate(across(where(is.double),as.character)) %>% 
        gather(.,key = set,value = "my_p_2",-bin,-set) %>% 
        dplyr::mutate(my_p_2=as.double(my_p_2))
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
    list_data$sortplot <- out_list1 %>% dplyr::mutate(set = gsub("(.{15})", "\\1\n", set))
    list_data$gene_file[[nick_name]]$info <- tibble(loaded_info =
      paste(
        "Filter Prob:",
        topbottom2,
        start_end_label[1],
        "to",
        start_end_label[2],
        "from",
        list_name,
        paste(file_names, collapse = " "),
        Sys.Date(),
        list_data$gene_file[[list_name]]$info
      ),
      save_name = gsub(" ", "_", paste("Filter",str_remove(topbottom2,"%"), Sys.Date(), sep = "_")),
      col_info = "gene"
      )
    list_data$gene_info <-
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>%
                           dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>%
                           dplyr::mutate(gene_list = nick_name,
                                         sub =  paste("Filter Prob:",
                                                      topbottom2,
                                                      start_end_label[1],
                                                      "to",
                                                      start_end_label[2]),
                                         onoff = "0",
                                         count = paste0("n = ", n_distinct(out_list$gene, na.rm = T)),
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
           start_end_label_peak,
           start_end_label_filter) {
    if (is.null(file_names)) {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return(NULL)
    }
    gene_list <- list_data$gene_file[[list_name]]$full
    out_list <- list_data$table_file %>% 
      dplyr::filter(set %in% file_names) %>% 
      semi_join(.,gene_list,by="gene") %>% 
      mutate(score=abs(score)) 
      
    my_filter <- out_list %>% 
      dplyr::filter(bin %in% start_end_bin_peak[1]:start_end_bin_peak[2]) %>% 
      group_by(set,gene) %>% summarise(score2=max(score,rm.na=T),.groups = "drop")
    out_list <- out_list %>% 
      dplyr::filter((bin %in% start_end_bin_filter[1]:start_end_bin_filter[2])) %>% 
      full_join(.,my_filter,by=c("gene","set")) 
    
    if(my_type == "peak"){
      out_gene <- out_list %>% 
          group_by(gene) %>% 
          filter(n_distinct(set)==length(file_names)) %>% 
          dplyr::filter(all(score<=score2)) %>% 
          ungroup() %>% distinct(gene)
    } else {
      out_gene <- out_list %>% 
          group_by(gene) %>% 
          dplyr::filter(!all(score<=score2)) %>%
          filter(n_distinct(set)==length(file_names)) %>%
          ungroup() %>% distinct(gene)
    }
    if (length(out_gene$gene) == 0) {
      return(NULL)
    }
    old_names <- grep("^Filter", names(list_data$gene_file), value = T)
    if (length(old_names) > 3) {
      old_names <- old_names[!old_names %in% list_name]
      # remove old sort gene list keeping 4
      list_data$gene_file[[first(old_names)]] <- NULL
      list_data$gene_info <- dplyr::filter(list_data$gene_info,
                                           gene_list != first(old_names))
    }
    nick_name <-
      strtrim(gsub("(.{30})",
                   "\\1... ",
                   paste0("Filter ",my_type, "\nn = ", n_distinct(out_gene$gene, na.rm = T))), 33)
    list_data$gene_file[[nick_name]]$full <- out_gene
    list_data$gene_file[[nick_name]]$info <- tibble(loaded_info =
                                                      paste(
                                                        "Filter:",
                                                        my_type,
                                                        "peak",
                                                        start_end_label_peak[1],
                                                        "to",
                                                        start_end_label_peak[2],
                                                        "filter peaks",
                                                        start_end_label_filter[1],
                                                        "to",
                                                        start_end_label_filter[2],
                                                        "from",
                                                        list_name,
                                                        paste(file_names, collapse = " "),
                                                        Sys.Date(),
                                                        list_data$gene_file[[list_name]]$info
                                                      ),
                                                    save_name = gsub(" ", "_", paste("Filter", my_type, Sys.Date(), sep = "_")),
                                                    col_info = "gene"
                                                    )
    list_data$gene_info <-
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>%
                           dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>%
                           dplyr::mutate(gene_list = nick_name,
                                         sub =  paste("Filter:",
                                                      my_type,
                                                      start_end_label_filter[1],
                                                      "to",
                                                      start_end_label_filter[2]),
                                         onoff = "0",
                                         count = paste0("n = ", n_distinct(out_gene$gene, na.rm = T)),
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
    if(dnom != "-1"){
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
    } else {
      if(nchar(nickname)<1){
        nickname <- paste0(nom,"*",dnom)
      }
      legend_nickname <- nickname
      new_gene_list <- list_data$table_file %>% dplyr::filter(set == nom) %>%
      replace_na(., list(score = 0)) %>% 
      dplyr::mutate(score=score*-1,set=legend_nickname)
    } 
    # adds meta data 
    list_data$table_file <- dplyr::filter(list_data$table_file, set != legend_nickname)
    list_data$gene_info <- dplyr::filter(list_data$gene_info, set != legend_nickname)
    list_data$table_file <- bind_rows(list_data$table_file, new_gene_list)
    list_data$gene_info <- bind_rows(list_data$gene_info,
                                     tibble(
                                       gene_list = names(list_data$gene_file)[1],
                                       count = distinct(list_data$gene_info %>% 
                                                          filter(gene_list == names(list_data$gene_file)[1]),count)$count,
                                       set = legend_nickname,
                                       group = legend_nickname,
                                       mycol = sample(suppressWarnings(brewer.pal(11, sample(kBrewerList,size=1)))[-c(4:7)],size = 1),
                                       onoff = "0",
                                       sub = " ",
                                       plot_set = " "
                                     ))
    return(list_data)
  }

# Total, antijoin and innerjoined gene lists
IntersectGeneLists <-
  function(list_data, list_name) {
    if (is.null(list_name)) {
      return(NULL)
    }
    outlist <- NULL
    # grab selected gene list(s)
    lapply(list_name, function(j) {
      outlist[[j]] <<- list_data$gene_file[[j]]$full %>% 
        dplyr::select(gene) %>% 
        dplyr::mutate(set=j)
    })
    # collapses into one list
    outlist <- bind_rows(outlist)
    # remove any pre used data
    list_data$gene_file <- list_data$gene_file[!str_detect(names(list_data$gene_file),"^Gene_List_")]
    list_data$gene_info <- list_data$gene_info %>% dplyr::filter(!str_detect(gene_list,"^Gene_List_"))
    
    # record for info
    if (n_distinct(outlist$gene) > 0) {
      nick_name1 <-
        paste("Gene_List_Total\nn =", n_distinct(outlist$gene, na.rm = T))
      # record for info
      list_data$gene_file[[nick_name1]]$full <- distinct(outlist,gene)
      list_data$gene_file[[nick_name1]]$info <- tibble(loaded_info =
        paste("Gene_List_Total",
              "from",
              paste(list_name, collapse = " and "),
              Sys.Date()),
        save_name = gsub(" ", "_", paste("Gene_List_Total_n=", n_distinct(outlist$gene, na.rm = T), Sys.Date(), sep = "_")),
        col_info = "gene"
        )
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste("Gene_List_Total"), 
                                           onoff = "0",
                                           plot_set = " ")))
    }
    innerjoined <- outlist %>% group_by(gene) %>% 
      filter(n_distinct(set)==length(list_name)) %>% 
      distinct(gene)
    if (n_distinct(innerjoined$gene) > 0) {
      nick_name1 <-
        paste("Gene_List_innerjoin\nn =", n_distinct(innerjoined$gene, na.rm = T))
      # record for info
      list_data$gene_file[[nick_name1]]$full <- innerjoined
      list_data$gene_file[[nick_name1]]$info <- tibble(loaded_info =
        paste("Gene_List_innerjoin",
              "from",
              paste(list_name, collapse = " and "),
              Sys.Date()),
        save_name = gsub(" ", "_", paste("Gene_List_innerjoin_n=", n_distinct(outlist$gene, na.rm = T), Sys.Date(), sep = "_")),
        col_info = "gene"
        )
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste("Gene_List_innerjoin"), 
                                           onoff = "0",
                                           plot_set = " ")))
      antijoin <- anti_join(distinct(outlist,gene), innerjoined, by = "gene")
      if (n_distinct(antijoin$gene) == 0) {
        antijoin <- distinct(outlist)
      }
      nick_name1 <-
        paste("Gene_List_antijoin\nn =", n_distinct(antijoin$gene, na.rm = T))
      # record for info
      list_data$gene_file[[nick_name1]]$full <- antijoin
      list_data$gene_file[[nick_name1]]$info <- tibble(loaded_info =
        paste("Gene_List_antijoin",
              "from",
              paste(list_name, collapse = " and "),
              Sys.Date()),
        save_name = gsub(" ", "_", paste("Gene_List_antijoin_n=", n_distinct(outlist$gene, na.rm = T), Sys.Date(), sep = "_")),
        col_info = "gene"
        )
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste("Gene_List_antijoin"), 
                                           onoff = "0",
                                           plot_set = " ")))
      
    }
    list_data
  }

# creates clusters using hclust.vector, from one active file
FindClusters <- function(list_data,
                         list_name,
                         clusterfile,
                         start_end_bin) {
  # print("find clusters")
  if (clusterfile == "") {
    showModal(modalDialog(
      title = "Information message",
      paste("No file selected to work on"),
      size = "s",
      easyClose = TRUE
    ))
    return(NULL)
  }
  df <-
      semi_join(dplyr::filter(list_data$table_file, set == clusterfile), 
                list_data$gene_file[[list_name]]$full, by = 'gene') 
    list_data$clust <- list()
    list_data$clust$cm <-
      hclust.vector(as.data.frame(spread(df, bin, score))[, c((start_end_bin[1]:start_end_bin[2]) + 2)], method = "ward")
  list_data$clust$full <- distinct(df, gene)
  list_data
}

# Pull out the number of groups of clusters 2-10
ClusterNumList <- function(list_data,
                           list_name,
                           clusterfile,
                           start_end_label,
                           my_num) {
  # print("cutree")
  if (is_empty(list_data$clust) | clusterfile == "") {
    showModal(modalDialog(
      title = "Information message",
      paste("No file selected to work on"),
      size = "s",
      easyClose = TRUE
    ))
    return(NULL)
  }
  if (n_distinct(LIST_DATA$clust$full) < as.numeric(my_num)) {
    showModal(modalDialog(
      title = "Information message",
      paste("Can't make more clusters than number of genes"),
      size = "s",
      easyClose = TRUE
    ))
    return(NULL)
  }
  list_data$gene_file <- list_data$gene_file[!str_detect(names(list_data$gene_file),"^Cluster_")]
  list_data$gene_info <- list_data$gene_info %>% dplyr::filter(!str_detect(gene_list,"^Cluster_"))
  gene_list <-
      dplyr::mutate(list_data$clust$full, cm = cutree(list_data$clust$cm, my_num))
  for (nn in 1:my_num) {
    outlist <- dplyr::filter(gene_list, cm == nn)
    nick_name <-
      paste(paste0("Cluster_", nn, "\nn ="), n_distinct(outlist$gene, na.rm = T))
    list_data$gene_file[[nick_name]]$full <- dplyr::select(outlist, gene)
    list_data$gene_file[[nick_name]]$info <- tibble(loaded_info =
      paste(
        nick_name,
        start_end_label[1],
        "to",
        start_end_label[2],
        "from",
        list_name,
        clusterfile,
        my_num,
        "Cluster_",
        "total",
        Sys.Date()
      ),
      save_name = gsub(" ", "_", paste("Cluster", nn, "of", my_num, Sys.Date(), sep = "_")),
      col_info = "gene"
      )
    list_data$gene_info <- 
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>% 
                           dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                           dplyr::mutate(gene_list = nick_name,
                                         sub =  paste(
                                           "Cluster_",
                                           "by",
                                           my_num,
                                           "from",
                                           list_name
                                         ), 
                                         onoff = "0",
                                         count = paste0("n = ", n_distinct(outlist$gene, na.rm = T)),
                                         plot_set = " ")))
  }
  list_data
}

# a[1]/b[2] or (a[1]/a[2])/(b[1]/b[2]) make gene list
CompareRatios <-
  function(list_data,
           list_name,
           ratio1file,
           ratio2file,
           startend1_bin,
           startend1_label,
           startend2_bin,
           startend2_label,
           my_num,
           divzerofix,
           normbin = "NA",
           normlabel = "NA") {
    if (ratio1file == "") {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return()
    }
    if(length(startend2_bin) < 2){
      startend2_bin <- c(0,0)
    }
    start1_bin <- startend1_bin[1]
    start2_bin <- startend2_bin[1]
    end1_bin <- startend1_bin[2]
    end2_bin <- startend2_bin[2]
    outlist <- NULL
    if (ratio2file == "None" | ratio2file == "") {
      if(start2_bin == 0){
        showModal(modalDialog(
          title = "Information message",
          paste("no file or bins to compare to"),
          size = "s",
          easyClose = TRUE
        ))
        return()
      }
      ratiofile <- ratio1file
      ratio2file <- paste0(ratio1file,
      "[", startend1_label[1], ":", startend1_label[2],
      "]/[", startend2_label[1], ":", startend2_label[2], "]")
    } else {
      ratiofile <- c(ratio1file, ratio2file)
      ratio2file <- paste0(ratio2file,
             "[", startend1_label[1], ":", startend1_label[2],
             "]/[", startend2_label[1], ":", startend2_label[2], "]/",
             ratio1file,
             "[", startend1_label[1], ":", startend1_label[2],
             "]/[", startend2_label[1], ":", startend2_label[2], "]")
    }
    lc <- 0
    lapply(ratiofile, function(j) {
      df <-
        semi_join(dplyr::filter(list_data$table_file, set == j), 
                  list_data$gene_file[[list_name]]$full, by = 'gene') 
      if (normlabel != "NA") {
        df <- group_by(df, gene) %>%
          arrange(bin) %>% 
          dplyr::mutate(score = score / nth(score, normbin))
      }
      df <- group_by(df, gene) %>%
        summarise(sum1 = sum(score[start1_bin:end1_bin],	na.rm = T),
                  sum2 = sum(score[start2_bin:end2_bin],	na.rm = T),.groups="drop") %>%
        ungroup()
      # if min/2 find Na's and 0's, and replace
      if(start2_bin == 0){
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
    if(my_num < 0){
      my_num <- 1/my_num
    }
    nick_name <- NULL
    if(my_num != 0){
      upratio <- dplyr::filter(outlist[[1]], Ratio > my_num)
    } else {
      upratio <- NULL
    }
    
    if (n_distinct(upratio$gene) > 0) {
      nick_name1 <-
        paste("Ratio_Up_file1\nn =", n_distinct(upratio$gene, na.rm = T))
      nick_name <- c(nick_name, nick_name1)
      list_data$gene_file[[nick_name1]]$full <- upratio 
      list_data$gene_file[[nick_name1]]$info <- tibble(loaded_info =
        paste(
          "Ratio_Up_file1",
          ratio2file,
          "fold change cut off",
          my_num,
          divzerofix,
          "from",
          list_name,
          "gene list",
          Sys.Date()
        ),
        save_name = gsub(" ", "_", paste("Ratio_Up_file1_fold_cut_off", my_num, Sys.Date(), sep = "_")),
        col_info = "gene file2/file1"
        )
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste(
                                             "\nRatio_Up_file1",
                                             "fold change cut off",
                                             my_num,
                                             divzerofix,
                                             "from",
                                             list_name
                                           ), 
                                           onoff = "0",
                                           count = paste0("n = ", n_distinct(upratio$gene, na.rm = T)),
                                           plot_set = " ")))
    }
    if(my_num != 0){
      upratio <- dplyr::filter(outlist[[1]], Ratio < 1 / my_num & Ratio != 0)
    } else {
      upratio <- NULL
    } 
    if (n_distinct(upratio$gene) > 0) {
      nick_name2 <-
        paste("Ratio_Down_file1\nn =", n_distinct(upratio$gene, na.rm = T))
      nick_name <- c(nick_name, nick_name2)
      list_data$gene_file[[nick_name2]]$full <- upratio
      list_data$gene_file[[nick_name2]]$info <- tibble(loaded_info =
        paste(
          "Ratio_Down_file1",
          ratio2file,
          "fold change cut off",
          my_num,
          divzerofix,
          "from",
          list_name,
          "gene list",
          Sys.Date()
        ),
        save_name = gsub(" ", "_", paste("Ratio_Down_file1_fold_cut_off", my_num, Sys.Date(), sep = "_")),
        col_info = "gene file2/file1"
        )
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name2,
                                           sub =  paste(
                                             "Ratio_Down_file1",
                                             "fold change cut off",
                                             my_num,
                                             divzerofix,
                                             "from",
                                             list_name
                                           ), 
                                           onoff = "0",
                                           count = paste0("n = ", n_distinct(upratio$gene, na.rm = T)),
                                           plot_set = " ")))
    }
    if(my_num != 0){
      upratio <-
        dplyr::filter(outlist[[1]], Ratio <= my_num &
                        Ratio >= 1 / my_num | Ratio == 0)
    } else {
      upratio <- outlist[[1]]
    }
    
    if (n_distinct(upratio$gene) > 0) {
      nick_name3 <-
        paste("Ratio_No_Diff\nn =", n_distinct(upratio$gene, na.rm = T))
      nick_name <- c(nick_name, nick_name3)
      list_data$gene_file[[nick_name3]]$full <- upratio
      list_data$gene_file[[nick_name3]]$info <- tibble(loaded_info =
        paste(
          "Ratio_No_Diff",
          ratio2file,
          "fold change cut off",
          my_num,
          divzerofix,
          "from",
          list_name,
          "gene list",
          Sys.Date()
        ),
        save_name = gsub(" ", "_", paste("Ratio_No_Diff_fold_cut_off", my_num, Sys.Date(), sep = "_")),
        col_info = "gene file2/file1"
        )
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name3,
                                           sub =  paste(
                                             "Ratio_No_Diff",
                                             "fold change cut off",
                                             my_num,
                                             divzerofix,
                                             "from",
                                             list_name
                                           ), 
                                           onoff = "0",
                                           count = paste0("n = ", n_distinct(upratio$gene, na.rm = T)),
                                           plot_set = " ")))
    }
    list_data$boxRatio <- NULL
    for (nn in nick_name) {
      
      list_data$boxRatio <- bind_rows(list_data$boxRatio,list_data$gene_file[[nn]]$full %>% dplyr::mutate(set=nn))
    }
    list_data
  }

# Cumulative Distribution data prep "PI / EI"
CumulativeDistribution <-
  function(list_data,
           onoffs,
           startend1_bin,
           startend1_label,
           startend2_bin,
           startend2_label) {
    # print("cdf function")
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
      # Complete within gene list and sum regions
      tf <- dplyr::filter(list_data$table_file, set %in% onoffs[[list_name]]) %>% 
        semi_join(., 
                  list_data$gene_file[[list_name]]$full, by = 'gene') 
      tf <- tf %>% dplyr::filter(bin == 1) %>% 
        group_by(gene) %>% 
        filter(n_distinct(set)==n_distinct(tf$set)) %>% 
        distinct(gene) %>% ungroup() %>% 
        semi_join(tf, ., by = 'gene')
      outlist[[list_name]] <- tf %>% 
        group_by(gene,set) %>%
        summarise(sum1 = mean(score[startend1_bin[1]:startend1_bin[2]],	na.rm = T),
                  sum2 = mean(score[startend2_bin[1]:startend2_bin[2]],	na.rm = T),.groups="drop") %>%
        dplyr::mutate(., value = sum1 / sum2) %>%
        dplyr::mutate(value=log2(value)) %>% 
        na_if(Inf) %>% na_if(-Inf) %>% group_by(gene) %>% 
        dplyr::filter(!any(is.na(value))) %>% ungroup() %>% 
        group_by(., set) %>%
        arrange(value) %>%
        dplyr::transmute(
          gene = gene,
          bin = row_number(),
          set = set,
          plot_set = paste(list_name, "-", gsub("(.{20})", "\\1\n", set)),
          value = value
        ) %>%
        ungroup()
    }
    
    # unlist and binds all together
    outlist <- bind_rows(outlist) %>% distinct() %>% arrange(bin)
    # removes top and bottom %
    if (sum(startend1_bin) > sum(startend2_bin)) {
      use_header <- "Log2 EI Cumulative plot"
    } else {
      use_header <- "Log2 PI Cumulative plot"
    }
    if (n_distinct(outlist$gene) > 0) {
      nick_name1 <- paste0("CDF ", use_header)
      list_data$gene_file[[nick_name1]]$full <- outlist
      list_data$gene_file[[nick_name1]]$info <- tibble(loaded_info =
        paste(
          use_header,
          "CDF",
          startend1_label[1],
          "to",
          startend1_label[2],
          "/",
          startend2_label[1],
          "to",
          startend2_label[2],
          "from",
          names(onoffs),
          "gene list(s)",
          paste(distinct(outlist, plot_set), collapse = " "),
          Sys.Date()
        ),
        save_name = gsub(" ", "_", paste("CDF", Sys.Date(), sep = "_")),
        col_info = "gene Rank_order sample, plot_legend, PI/EI"
        )
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
                                           count = paste("n =", outlist %>% dplyr::filter(grepl(list_name,plot_set)) %>% 
                                                           summarise(n=n_distinct(bin))),
                                           plot_set = paste(list_name, "-", gsub("(.{20})", "\\1\n", set)),
                                           myheader = use_header)))
    }
    list_data
  }
# CDG ggplot function
GGplotC <-
  function(df2,
           plot_options,
           use_header,
           plot_range) {
    # print("ggplot CDF")
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

# AUC using the trapezoidal rule
AUC <- function(y){
  n <- length(y)
  0.5*(y[1]+y[n]+2*sum(y[-c(1,n)]))
}

# ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
