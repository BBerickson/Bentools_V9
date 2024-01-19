
# Brewer color sets to be available ----
kBrewerList <-
  c("Set1","Paired","Dark2","Spectral")

# lines and labels types ----
kLinesandlabels <- c(
  "543",
  "5",
  "3",
  "5L",
  "PI"
)

# math options available ----
kMathOptions <- c("mean", "sum", "median")

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

# settings for tss tes names
LinesLabelsSetNames <- function(mytype){
  if(mytype == "543"){
    myname <- c("TSS","pA")
  } else if(mytype == "5"){
    myname <- c("TSS","")
  } else if(mytype == "5L"){
    myname <- c("TSS","")
  } else if(mytype == "3"){
    myname <- c("","pA")
  } else if(mytype == "PI"){
    myname <- c("TSS","body")
  } else {
    myname <- c("start","end")
  }
}  
