## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)


## Data ------------------------------------------------------------------------
folder = "/Users/jiancao/Downloads/dataverse_files/all_files_in_csv_format/"
data = read.csv(paste0(folder, "CNR_QJE_apple.csv"))
iphone4 = read.csv(paste0(folder, "iPhone4.csv"))
exchange = read.csv(paste0(folder, "exchange_long.csv"))


## Functions -------------------------------------------------------------------
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

minmax.normalize = function(data, reference = NULL){
  if (is.null(reference)) {
    minimum = min(data)
    maximum = max(data)
  }else{
    minimum = min(reference)
    maximum = max(reference)
  }
  res = (data - minimum)/(maximum - minimum)
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

# smooth speed profile
smooth_profile = function(W){
  n = length(W)
  for (i in 2:(n-1)) {
    if (abs(W[i]) < min(abs(W[i - 1]), abs(W[i + 1]))) {
      ind_min = which(c(abs(W[i - 1]), abs(W[i + 1])) == 
                        min(abs(W[i - 1]), abs(W[i + 1])))[1]
      W[i] = c(W[i - 1], W[i + 1])[ind_min]
    }
  }
  return(W)
}


## Run -------------------------------------------------------------------------
data$date = as.Date(data$date, format = "%d %b %y")
exchange$date = as.Date(exchange$date, format = "%d %b %y")

# df = NULL
# for (i in 1:nrow(iphone4)) {
#   target = iphone4[i,]
#   models = strsplit(target$Model, split = ", ")[[1]]
#   res = data %>% filter(id %in% models)
#   res$a.number = target$A.Number
#   res$identifier = target$Identifier
#   res$finish = target$Finish
#   res$storage = target$Storage
#   name = paste0(target$Identifier, "-", target$Finish, "-", target$Storage)
#   df[[name]] = res
#   print(name)
# }
# 
# df_8gb = rbind(df[[1]], df[[4]]) %>% 
#   group_by(country, date) %>% 
#   summarise(price = mean(price, na.rm = TRUE))
# 
# df_8gb %>% 
#   ggplot(aes(x = date, y = price, color = country)) +
#   geom_line()

exchange = exchange %>% 
  group_by(country_string) %>% 
  mutate(index_ratio = usdx/first(usdx)) %>% 
  ungroup()

df = left_join(data, exchange, by = c("country" = "country_string", 
                                    "date" = "date"))

df = df %>% 
  group_by(id, country, idc) %>% 
  mutate(price_adj = price*usdx,
         relative_price_adj = price_adj/first(price_adj)) %>% 
  ungroup() %>% 
  group_by(country) %>% 
  mutate(diff_relative_price_adj = relative_price_adj - lag(relative_price_adj)) %>% 
  ungroup()

df = df %>% 
  group_by(country, date) %>% 
  summarise(avg_delta = mean(diff_relative_price_adj, na.rm = TRUE))

# df = df %>% 
#   mutate(index_adj = index/index_ratio)

df_na = df %>% 
  group_by(country) %>% 
  summarise(count_na = sum(is.na(delta_price)))

df = df %>% 
  filter(!(country %in% c("ph", "cd", "pl", "xf")))

df %>%
  filter(date <= as.Date("2012-01-01")) %>% 
  ggplot(aes(x = date, y = avg_delta, color = country)) +
  geom_line()

countries = unique(df$country)
df_a = df %>% filter(country == "jp")
df_b = df %>% filter(country == "cn")


align = dtw::dtw(diff(t.normalize(df_a$avg_delta)),
                 diff(t.normalize(df_b$avg_delta)),
                 keep = TRUE,
                 step.pattern = dtw::symmetricP2)
dtw::dtwPlotThreeWay(align) 
P = Matrix::sparseMatrix(align$index1,
                         align$index2)
W = warp2weight(P)
a = fitted(forecast::ets(W))
plot(ts(a))

