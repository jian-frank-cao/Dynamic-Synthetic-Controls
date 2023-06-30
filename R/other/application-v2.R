## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)


## Data ------------------------------------------------------------------------
symbols = c("TIP", "INXG", "XRB", "IFLD", "ILB", "IBGM")
data_list <- list()

# Loop over the symbols and download the data
for(symbol in symbols) {
  data_list[[symbol]] <- quantmod::getSymbols(symbol, auto.assign = FALSE)
}

quantmod::getSymbols("AAPL")


# Install and load the package
install.packages("tidyquant")
library(tidyquant)

# Get data
data <- tq_get(c("AAPL", "GOOG", "FB"), 
               from = "2020-01-01", 
               to = "2020-12-31")


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

# plot index
oecd = c("at", "ca", "dk", "fr", "de", "gr", "is", "ie", "be",
         "lu", "nl", "no", "pt", "es", "se", "tr", "it", "ch",
         "uk", "us")

df_oecd = df %>% filter(country %in% c("jp", oecd))

range_oecd = df_oecd %>%
  group_by(country) %>%
  summarise(min_date = min(date),
            max_date = max(date))

select_country = range_oecd %>%
  filter(min_date == as.Date("2009-04-01")) %>%
  .[["country"]]

df %>%
  filter(country %in% oecd) %>%
  ggplot(aes(x = date, y = index, color = country)) +
  geom_line()


df_a = df %>% filter(country == "jp")
df_b = df %>% filter(country == "uk")


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


## Exchange --------------------------------------------------------------------
# countries = readRDS("./data/country_selected.Rds")
# 
# euro_diff = exchange %>% 
#   filter(country_string %in% countries) %>% 
#   group_by(country_string) %>% 
#   mutate(diff_euro = usdx - euro) %>% 
#   summarise(avg = mean(diff_euro, na.rm = T),
#             sd = sd(diff_euro, na.rm = T))
# 
# count_sd = data.frame(table(euro_diff$sd)) %>% 
#   arrange(desc(Freq))
# 
# 
# chf_countries = euro_diff %>% filter(sd == count_sd$Var[9]) %>% 
#   .[["country_string"]]
# 
# mad_countries = euro_diff %>% filter(sd == count_sd$Var[7]) %>% 
#   .[["country_string"]]
# 
# nok_countries = euro_diff %>% filter(sd == count_sd$Var[6]) %>% 
#   .[["country_string"]]
# 
# dkk_countries = euro_diff %>% filter(sd == count_sd$Var[5]) %>% 
#   .[["country_string"]]
# 
# gbp_countries = euro_diff %>% filter(sd %in% count_sd$Var[c(4,8)]) %>% 
#   .[["country_string"]]
# 
# nzd_countries = euro_diff %>% filter(sd == count_sd$Var[3]) %>% 
#   .[["country_string"]]
# 
# aud_countries = euro_diff %>% filter(sd == count_sd$Var[2]) %>% 
#   .[["country_string"]]
# 
# euro_countries = euro_diff %>% filter(sd == count_sd$Var[1]) %>% 
#   .[["country_string"]]
# 
# currency = data.frame(country = countries) %>% 
#   mutate(currency = case_when(country %in% euro_countries ~ "EURO",
#                               country %in% aud_countries ~ "AUD",
#                               country %in% nzd_countries ~ "NZD",
#                               country %in% gbp_countries ~ "GBP",
#                               country %in% dkk_countries ~ "DKK",
#                               country %in% nok_countries ~ "NOK",
#                               country %in% mad_countries ~ "MAD",
#                               country %in% chf_countries ~ "CHF",
#                               TRUE ~ country))
# 
# saveRDS(currency, "./data/currency.Rds")

currency = readRDS("./data/currency.Rds")

targets = currency %>% 
  group_by(currency) %>% 
  summarise(country = first(country))

df = inner_join(exchange, targets, by = c("country_string" = "country"))

df %>%
  filter(currency %in% c("EURO", "GBP", "AUD", "NZD",
                         "DKK", "NOK", "MAD", "CHF")) %>%
  ggplot(aes(x = date, y = usdx, color = currency)) +
  geom_line()


df_a = df %>% filter(currency == "EURO")
df_b = df %>% filter(currency == "AUD")


align = dtw::dtw(diff(t.normalize(df_a$usdx)),
                 diff(t.normalize(df_b$usdx)),
                 keep = TRUE,
                 step.pattern = dtw::symmetric2)
dtw::dtwPlotThreeWay(align) 
P = Matrix::sparseMatrix(align$index1,
                         align$index2)
W = warp2weight(P)
a = fitted(forecast::ets(W))
plot(ts(a))

