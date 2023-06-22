## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)


## Data ------------------------------------------------------------------------
folder = "/Users/jiancao/Downloads/dataverse_files/all_files_in_csv_format/"
exchange = read.csv(paste0(folder, "exchange_long.csv"))


## Functions -------------------------------------------------------------------
# prepare one company
one_company = function(path, exchange){
  data = read.csv(path)
  
  # prepare data
  data$date = as.Date(data$date, format = "%d %b %y")
  exchange$date = as.Date(exchange$date, format = "%d %b %y")
  
  # compute diff of log price
  df = left_join(data, exchange, by = c("country" = "country_string", 
                                        "date" = "date")) %>% 
    mutate(price_usd = price*usdx,
           log_price_usd = log(price_usd)) %>% 
    group_by(country) %>% 
    mutate(diff_log = log_price_usd - lag(log_price_usd)) %>% 
    ungroup()
  
  # filter outliers and compute average diff_log
  df = df %>% 
    filter(diff_log > quantile(diff_log, 0.25, na.rm = TRUE) -
             3*IQR(diff_log, na.rm = TRUE) &
             diff_log < quantile(diff_log, 0.75, na.rm = TRUE) +
             3*IQR(diff_log, na.rm = TRUE)) %>% 
    group_by(country, date) %>% 
    summarise(avg_log = mean(diff_log, na.rm = TRUE))
  
  return(df)
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
df_apple = one_company(paste0(folder, "CNR_QJE_apple.csv"),
                       exchange = exchange)
gc()
df_ikea = one_company(paste0(folder, "CNR_QJE_ikea.csv"),
                       exchange = exchange)
gc()
df_handm = one_company(paste0(folder, "CNR_QJE_handm.csv"),
                       exchange = exchange)
gc()
df_zara = one_company(paste0(folder, "CNR_QJE_zara.csv"),
                       exchange = exchange)
gc()

# compute index
df = df %>% 
  group_by(country) %>% 
  mutate(cum_log = cumsum(avg_log),
         index = exp(cum_log))

# plot index
oecd = c("at", "be", "ca", "dk", "fr", "de", "gr", "is", "ie",
         "it", "lu", "nl", "no", "pt", "es", "se", "ch", "tr",
         "uk", "us")
df %>%
  filter(country %in% oecd) %>%
  ggplot(aes(x = date, y = index, color = country)) +
  geom_line()

countries = unique(df2$country)
df_a = df2 %>% filter(country == "jp")
df_b = df2 %>% filter(country == "us")


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

