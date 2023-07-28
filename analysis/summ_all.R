library(openxlsx)
library(tidyverse)
library(magrittr)
library(rio)
library(rlang)

xlsx <- import("K:\\projects\\alindsay\\Projects\\wwCOV\\results\\ncov-ww_upto_230728.csv")
summCols <- xlsx %>% select(A_x:XCG) %>% colnames() 
lineCols <- xlsx %>% select(A_y:'XBB.2.3* [Omicron (XBB.2.3.X)]') %>% colnames()

summ <- xlsx# %>% filter(coverage_x > 80)
summ %<>% mutate(Max.Lineage=pmax(!!!rlang::syms(lineCols),na.rm=TRUE))
#summ %<>% filter(Max.Lineage < 95)
all_na <- function(x) any(!is.na(x))
summ %<>% select_if(all_na)

df <- summ

max_list <- function(x) {
  max <- x[!is.na(x)]
  #max <- as.list(x %<>% select_if(~ !any(is.na(.))))
  max <- head(max[order(max,decreasing=TRUE)],3)
  #max <- max[order(names(max),decreasing=TRUE)][1:3]
  cols <- names(max)
  vals <- unlist(unname(max))

  out <- mapply(paste, cols, vals,  MoreArgs = list(sep = ";"))
  out <- paste(out,collapse=";")
  return (out)
}

summCols <- intersect(summCols, colnames(df))
lineCols <- intersect(lineCols, colnames(df))


df['SummStats'] = apply(df[summCols],1,max_list)
df['LineStats'] = apply(df[lineCols],1,max_list)

out <- df

out %<>% separate_wider_delim(SummStats, ";", names = c("Summ. 1", "Summ. 1 (%)","Summ. 2", "Summ. 2 (%)","Summ. 3", "Summ. 3 (%)"), too_few = "align_start")
out %<>% separate_wider_delim(LineStats, ";", names = c("Line. 1", "Line. 1 (%)","Line. 2", "Line. 2 (%)","Line. 3", "Line. 3 (%)"), too_few = "align_start")

write.xlsx(out,"K:\\projects\\alindsay\\Projects\\wwCOV\\results\\ncov-ww_upto_230728.xlsx")

