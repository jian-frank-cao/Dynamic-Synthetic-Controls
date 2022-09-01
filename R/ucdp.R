library(checkpoint)
checkpoint("2022-04-01")
library(tidyverse)
temp = read.csv("/Users/jiancao/Downloads/GEDEvent_v22_1.csv", nrows = 1)
col_names = colnames(temp)
col_classes = rep("NULL", length(col_names))
col_classes[col_names %in% c("priogrid_gid", "best")] = "integer"
col_classes[col_names %in% c("date_start")] = "character"

data = data.table::fread("/Users/jiancao/Downloads/GEDEvent_v22_1.csv",
                select = c("best", "date_start", "country")) %>% 
  mutate(date = as.Date(date_start)) %>% 
  select(-date_start)
AllDates = expand.grid(country = unique(data$country), 
                       date = as.Date(as.Date("1989-01-01"):as.Date("2021-12-31"),
                                      origin = "1970-01-01"))
data = left_join(AllDates, data, by = c("country", "date"))
data[is.na(data)] = 0
data = data %>% 
  mutate(month = substr(date, 1, 7)) %>% 
  group_by(country, month) %>% 
  summarize(best = sum(best)) %>% 
  arrange(country, month) %>% 
  mutate(time = 1:((2021 - 1989 + 1)*12))

data %>% 
  filter(country != "Rwanda") %>%
  filter(time %in% 1:100) %>%
  ggplot(aes(x = time, y = best, color = country)) +
  geom_line() +
  theme(legend.position="none")
