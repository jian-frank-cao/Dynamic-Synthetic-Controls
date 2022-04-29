for (i in (1:10)*2+3) {
  for (j in 2:6) {
    for (z in 1:5) {
      if (z > j | i <= j) {
        next
      }
      gdp = reshape2::dcast(data %>% select(c("country", "year", "gdp")), year ~ country, value.var = "gdp")
      # row.names(gdp) = gdp$year
      # gdp = gdp %>% select(-c("Norway", "Portugal"))
      
      ## minmax
      # gdp = gdp %>% mutate_at(setdiff(colnames(gdp), "year"), ~normalize(., "minmax"))
      
      ## t
      gdp = gdp %>% mutate_at(setdiff(colnames(gdp), "year"), ~normalize(., "t"))
      
      ## adding buffer
      add_buffer = function(TS, n){
        model_right = forecast::auto.arima(TS)
        right <- as.numeric(forecast::forecast(model_right, h = n)$mean)
        model_left = forecast::auto.arima(rev(TS))
        left <- rev(as.numeric(forecast::forecast(model_left, h = n)$mean))
        
        return(c(left, TS, right))
      }
      
      width = i
      gdp_b = sapply(gdp %>% select(-year), add_buffer, n = (width - 1)/2) %>% 
        data.frame(.)
      
      ## derivative
      gdp_b = gdp_b %>%
        mutate_all(~signal::sgolayfilt(., j, width, z)) %>%
        .[((width - 1)/2 + 1):((width - 1)/2 + nrow(gdp)),]
      gdp[-1] = gdp_b
      
      res = mean(mean(var(gdp_b)))
      
      print(paste0(c(i, j, z, res), collapse = "-"))
    }
  }
}

gdp = reshape2::dcast(data %>% select(c("country", "year", "gdp")), year ~ country, value.var = "gdp")
# row.names(gdp) = gdp$year
# gdp = gdp %>% select(-c("Norway", "Portugal"))

## minmax
# gdp = gdp %>% mutate_at(setdiff(colnames(gdp), "year"), ~normalize(., "minmax"))

## t
gdp = gdp %>% mutate_at(setdiff(colnames(gdp), "year"), ~normalize(., "t"))

## adding buffer
add_buffer = function(TS, n){
  model_right = forecast::auto.arima(TS)
  right <- as.numeric(forecast::forecast(model_right, h = n)$mean)
  model_left = forecast::auto.arima(rev(TS))
  left <- rev(as.numeric(forecast::forecast(model_left, h = n)$mean))
  
  return(c(left, TS, right))
}

width = 9
gdp_b = sapply(gdp %>% select(-year), add_buffer, n = (width - 1)/2) %>% 
  data.frame(.)

## derivative
gdp_b = gdp_b %>%
  mutate_all(~signal::sgolayfilt(., 3, width, 2)) %>%
  .[((width - 1)/2 + 1):((width - 1)/2 + nrow(gdp)),]
gdp[-1] = gdp_b

## detrend
# library(remotes)
# remotes::install_github("jrevenaugh/TSAUMN")
# gdp = gdp %>% mutate_at(setdiff(colnames(gdp), "year"), ~TSAUMN::detrend(., 3))


df <- reshape2::melt(gdp ,  id.vars = 'year', variable.name = 'country')


# df = df %>% filter(!(country %in% c("Portugal", "New Zealand", "Norway", "West Germany",
#                                     "UK", "USA")))
# plot on same grid, each series colored differently -- 
# good if the series have same scale
ggplot(df, aes(year,value)) + geom_line(aes(colour = country))

# ------------------------------------------------------------------------------

load("./data/smoking.rda")
prop99_raw = read.csv("./data/prop99.csv")

prop99 = prop99_raw %>% 
  filter(SubMeasureDesc == 'Cigarette Consumption (Pack Sales Per Capita)') %>% 
  reshape2::dcast(., Year ~ LocationDesc, value.var = "Data_Value")



ggplot(smoking, aes(year,cigsale)) + geom_line(aes(colour = as.factor(state)))



