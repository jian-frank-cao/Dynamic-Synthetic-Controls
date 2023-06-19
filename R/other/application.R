## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)


## Data ------------------------------------------------------------------------
folder = "/Users/jiancao/Downloads/dataverse_files/"
data = haven::read_dta(paste0(folder, "all_files_in_csv_format/CNR_QJE_apple.csv"))
iphone4 = read.csv(paste0(folder, "iPhone4.csv"))


## Run -------------------------------------------------------------------------
data$date = data$date

for (i in 1:nrow(iphone4)) {
  target = iphone4[i,]
  models = strsplit(target$Model, split = ", ")
  res = data %>% filter(id %in% models)
}