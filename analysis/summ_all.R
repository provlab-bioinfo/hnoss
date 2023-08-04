library(openxlsx)
library(tidyverse)
library(magrittr)
library(rio)
library(rlang)

xlsx <- import("K:\\projects\\alindsay\\Projects\\wwCOV\\results\\mixed_samples_big_BN2.csv")
summCols <- xlsx %>% select(A_x:Zeta) %>% colnames() 
lineCols <- xlsx %>% select(A_y:Z.1) %>% colnames()

summ <- xlsx %>% filter(coverage_x > 80)
summ %<>% mutate(Max.Lineage=pmax(!!!rlang::syms(lineCols),na.rm=TRUE))
summ %<>% filter(Max.Lineage < 95)
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

SummCols <- c("Summ. 1", "Summ. 1 (%)","Summ. 2", "Summ. 2 (%)","Summ. 3", "Summ. 3 (%)")
LineCols <- c("Line. 1", "Line. 1 (%)","Line. 2", "Line. 2 (%)","Line. 3", "Line. 3 (%)")

out %<>% separate_wider_delim(SummStats, ";", names = SummCols, too_few = "align_start")
out %<>% separate_wider_delim(LineStats, ";", names = LineCols, too_few = "align_start")

writeClipboard(paste(colnames(out),collapse="\n"))

keep_cols <- c("file","Key","fastaPath","collection_date","ct","resid_x","coverage_x","lineage","scorpio_call","ncov_quality",
               "Summ. 1","Summ. 1 (%)","Summ. 2","Summ. 2 (%)","Summ. 3","Summ. 3 (%)","Line. 1",
               "Line. 1 (%)","Line. 2","Line. 2 (%)","Line. 3","Line. 3 (%)")

# out %<>% select(any_of(keep_cols)) %>% rename(resid = resid_x, coverage = coverage_x)



write.xlsx(out,"K:\\projects\\alindsay\\Projects\\wwCOV\\results\\mixed_samples_big_BN4.xlsx")

