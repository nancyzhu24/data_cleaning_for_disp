# R function to retrieve ATC codes

library(rvest)
library(stringr)

whoatc <- function(code, level = all){
  #if(!as.numeric(level) %in% c(1:5)){return("Please include a level from 1 to 5")}
  atclink <- "http://www.whocc.no/atc_ddd_index/"
  s <- html_session(atclink)
  atchtml <- html_form(s)[[2]] %>%
    set_values(code = code) %>%
    {submit_form(s, .)} %>%
    read_html() %>%
    html_nodes("td a, b a")
  codes <- atchtml %>%
    html_attr('href') %>%
    str_extract(regex("(?<=code=)\\w*"))
  desc <- atchtml %>%
    html_text(trim = TRUE)
  return(data.frame(code = codes, desc = desc, stringsAsFactors = FALSE))
}