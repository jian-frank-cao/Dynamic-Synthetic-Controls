## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)


## Data ------------------------------------------------------------------------
folder = "/Users/jiancao/Downloads/dataverse_files/all_files_in_csv_format/"
exchange = data.table::fread(paste0(folder, "exchange_long.csv"))
exchange$date = as.Date(exchange$date, format = "%d %b %y")


## Functions -------------------------------------------------------------------
# prepare one company
one_company = function(path, exchange){
  data = data.table::fread(path, select = 1:5)
  
  # prepare date
  data$date = as.Date(data$date, format = "%d %b %y")
  gc()
  
  # compute diff of log price
  data = left_join(data, exchange, by = c("country" = "country_string", 
                                        "date" = "date")) %>% 
    group_by(country, id) %>% 
    mutate(diff_log = log(price) - lag(log(price))) %>% 
    # mutate(diff_log = log(price*usdx) - lag(log(price*usdx))) %>% 
    ungroup()
  gc()
  
  selection = data %>% 
    group_by(id, country) %>% 
    summarise(ratio = max(price, na.rm = T)/min(price, na.rm = T)) %>% 
    filter(country != "ph" & ratio > 1 & ratio < 1000) %>% 
    select(-ratio)
  
  data = inner_join(data, selection, by = c("country", "id"))
  
  # filter outliers and compute average diff_log
  data = data %>% 
    # filter(diff_log >= quantile(diff_log, 0.25, na.rm = TRUE) -
    #          3*IQR(diff_log, na.rm = TRUE) &
    #          diff_log <= quantile(diff_log, 0.75, na.rm = TRUE) +
    #          3*IQR(diff_log, na.rm = TRUE)) %>% 
    group_by(country, date) %>% 
    summarise(avg_log = mean(diff_log, na.rm = TRUE))
  
  return(data)
}

# normalization
t.normalize = function(data, reference = NULL){
  if (is.null(reference)) {
    mu = mean(data)
    sigma = sd(data)
  }else{
    mu = mean(reference)
    sigma = sd(reference)
  }
  res = (data - mu)/sigma
  return(res)
}

# transform warping path W to weight
warp2weight = function(W){
  w = as.matrix(W)
  count = rep(1/colSums(w), nrow(w)) %>% 
    matrix(.,
           nrow = ncol(w),
           ncol = nrow(w)) %>% 
    t(.)
  weight = rowSums(w * count)
  
  return(weight)
}


## Run -------------------------------------------------------------------------
# prepare data
# df_apple = one_company(paste0(folder, "CNR_QJE_apple.csv"),
#                        exchange = exchange)
# df_apple$avg_log[is.na(df_apple$avg_log)] = 0
# gc()
# df_ikea = one_company(paste0(folder, "CNR_QJE_ikea.csv"),
#                        exchange = exchange)
# df_ikea$avg_log[is.na(df_ikea$avg_log)] = 0
# gc()
# df_handm = one_company(paste0(folder, "CNR_QJE_handm.csv"),
#                        exchange = exchange)
# df_handm$avg_log[is.na(df_handm$avg_log)] = 0
# gc()
# df_zara = one_company(paste0(folder, "CNR_QJE_zara.csv"),
#                        exchange = exchange)
# df_zara$avg_log[is.na(df_zara$avg_log)] = 0
# gc()
# 
# saveRDS(df_apple, "./data/df_apple.Rds")
# saveRDS(df_ikea, "./data/df_ikea.Rds")
# saveRDS(df_handm, "./data/df_handm.Rds")
# saveRDS(df_zara, "./data/df_zara.Rds")

df_apple = readRDS("./data/df_apple.Rds")
df_ikea = readRDS("./data/df_ikea.Rds")
df_handm = readRDS("./data/df_handm.Rds")
df_zara = readRDS("./data/df_zara.Rds")

# compute index
df = rbind(df_apple, df_ikea, df_handm, df_zara) %>% # df_ikea, 
  group_by(country, date) %>% 
  summarise(avg_log = mean(avg_log, na.rm = TRUE)) %>% 
  ungroup() %>% 
  group_by(country) %>% 
  mutate(cum_log = cumsum(avg_log),
         index = exp(cum_log))

# # plot index
# oecd = c("at", "ca", "dk", "fr", "de", "gr", "is", "ie", "be", 
#          "lu", "nl", "no", "pt", "es", "se", "tr", "it", "ch", 
#          "uk", "us")
# 
# df_oecd = df %>% filter(country %in% c("jp", oecd))
# 
# range_oecd = df_oecd %>% 
#   group_by(country) %>% 
#   summarise(min_date = min(date),
#             max_date = max(date))
# 
# select_country = range_oecd %>%
#   filter(min_date == as.Date("2009-04-01")) %>% 
#   .[["country"]]

# df %>%
#   filter(country %in% oecd) %>%
#   ggplot(aes(x = date, y = index, color = country)) +
#   geom_line()
# 
# 
# df_a = df %>% filter(country == "jp")
# df_b = df %>% filter(country == "uk")
# 
# 
# align = dtw::dtw(diff(t.normalize(df_a$index)),
#                  diff(t.normalize(df_b$index)),
#                  keep = TRUE,
#                  step.pattern = dtw::symmetricP2)
# dtw::dtwPlotThreeWay(align) 
# P = Matrix::sparseMatrix(align$index1,
#                          align$index2)
# W = warp2weight(P)
# a = fitted(forecast::ets(W))
# plot(ts(a))

select_country = exchange %>% 
  group_by(country_string) %>%
  summarise(sd = sd(usdx),
            avg = mean(usdx),
            ratio = sd/avg) %>% 
  filter(ratio > 0.1)

exchange %>%
  filter(country_string %in% select_country$country_string) %>%
  ggplot(aes(x = date, y = usdx, color = country_string)) +
  geom_line()




df_a = exchange %>% filter(country_string == "jp")
df_b = exchange %>% filter(country_string == "uk")


align = dtw::dtw(diff(t.normalize(df_a$index)),
                 diff(t.normalize(df_b$index)),
                 keep = TRUE,
                 step.pattern = dtw::symmetricP2)
dtw::dtwPlotThreeWay(align) 
P = Matrix::sparseMatrix(align$index1,
                         align$index2)
W = warp2weight(P)
a = fitted(forecast::ets(W))
plot(ts(a))

