# setwd("F:/Trinity/Dynamic-Synthetic-Control")
source("./R/TwoStepDTW.R")
source("./R/synthetic_control.R")
set.seed(20220407)

## Data ------------------------------------------------------------------------
data = foreign::read.dta("data/repgermany.dta")

# # adding sin curves
# frequency = data %>% select(country) %>% distinct
# frequency$freq = round(runif(nrow(frequency)) * 5 + 6, 1)
# 
# data = left_join(data, frequency, by = "country")
# data = data %>% group_by(country) %>% 
#   mutate(x = seq(0, 3.14*mean(freq), len = length(freq)),
#          sinx = sin(x),
#          gdp_bak = gdp,
#          gdp = gdp + sinx*2000) %>% 
#   ungroup %>% 
#   mutate_at(vars(index, year, gdp, infrate, trade, schooling,
#                  invest60, invest70, invest80,
#                  industry),
#             as.numeric) %>%
#   mutate_at(vars(country), as.character)

# data = data %>% filter(!(country %in% c("Austria", "Japan")))


# # quarterly gdp
# qtrly_gdp = read.csv("./data/quarterly_gdp.csv")
# population = read.csv("./data/population.csv")
# 
# population = reshape2::melt(population[1:17, -c(1,2,4)],
#                             id.vars = "Country.Name")
# population = population %>% 
#   mutate(
#     country = Country.Name,
#     year = as.numeric(substr(variable, 2, 5)),
#     population = value
#   ) %>% 
#   select(c("country", "year", "population"))
#   
# data_qtrly = qtrly_gdp %>% 
#   mutate(
#     country = Country,
#     year = as.numeric(substr(TIME, 1, 4)),
#     quarter = as.numeric(substr(TIME, 7, 7)),
#     gdp_total = Value * 1000000
#   ) %>% 
#   select(c("country", "year", "quarter", "gdp_total"))
# 
# data_qtrly = left_join(data_qtrly, population,
#                        by = c("country", "year")) %>% 
#   mutate(
#     gdp = round(gdp_total/population, 2)
#   )
# 
# data_qtrly = left_join(data_qtrly,
#                        data %>% mutate(country = case_when(country == "USA" ~ "United States",
#                                                            country == "UK" ~ "United Kingdom",
#                                                            country == "West Germany" ~ "Germany",
#                                                            TRUE ~ country)) %>% 
#                                          select(-gdp),
#                        by = c("country", "year"))


# preprocessing
gdp = reshape2::dcast(data %>% select(c("country", "year", "gdp")),
                      year ~ country, value.var = "gdp")
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
colnames(df)[3] = "gdp_ready"

data = left_join(data, df, by = c("country", "year"))
data = data %>% mutate(gdp_original = gdp,
                       gdp = gdp_ready)


## Function --------------------------------------------------------------------
stretch_var = function(Country, col_name, stretch){
  target = data %>%
    filter(country == Country) %>% 
    .[[col_name]]
  len = length(target)
  if (stretch == 1) {   # fix this
    res = target
  }else{
    res = approx(c(0, target, 0),
                 n = (len + 1) * stretch + 1,
                 # na.rm = FALSE,
                 method = "linear"
    ) %>% 
      .[["y"]] %>% 
      .[(stretch + 1):(len * stretch + 1)]
  }
  
  return(res)
}

compare_methods = function(data,
                           last_year = 2003,
                           dependent = "West Germany",
                           dependent_id = 7,
                           stretch = 5,
                           t_treat = (21 - 1) * stretch + 1,
                           n = (44 - 1) * stretch + 1,
                           n_1st = (30 - 1) * stretch + 1,
                           k = 4 * stretch,
                           n_q = 1,
                           n_r = 1,
                           normalize_method = "minmax",
                           dtw_method = "dtw",
                           margin = 10, ...){
  # prepare data
  y_raw = data %>% 
    filter(country == dependent & year <= last_year) %>%
    .[["gdp"]]
  
  y_original = data %>% 
    filter(country == dependent & year <= last_year) %>%
    .[["gdp_original"]]
  
  x_list = data %>% 
    filter(country != dependent & year <= last_year) %>% 
    select(c("country", "gdp", "gdp_original")) %>% 
    group_by(country) %>% 
    group_split(.keep = TRUE)
  
  # Two Step DTW
  # results = x_list %>% 
  #   future_map(
  #     ~{
  #       item = .
  #       country = item[["country"]][1]
  #       x = item[["gdp"]]
  #       
  #       y = approx(y_raw, n = (length(x) - 1) * stretch + 1, method = "linear")$y
  #       x = approx(x, n = (length(x) - 1) * stretch + 1, method = "linear")$y
  #       
  #       res = TwoStepDTW(y, x, t_treat, k, n_1st, dtw_method = dtw_method,
  #                        normalize_method = normalize_method, ...)
  #       
  #       df = data.frame(time = 1:length(x), y = res$y[1:length(x)], x = res$x[1:length(x)], warped = res$x_w[1:length(x)]) %>% 
  #         `colnames<-`(c("time", dependent, country, paste0(country, "-Warped"))) %>% 
  #         reshape2::melt(., id.vars = "time") %>% 
  #         `colnames<-`(c("Time", "Country", "GDP"))
  #       
  #       fig = df %>% 
  #         ggplot(aes(x = Time, y = GDP, color = Country)) +
  #         geom_line() +
  #         scale_color_manual(values = c("#2a4d69", "#ee4035", "#7bc043")) +
  #         geom_vline(xintercept = t_treat, linetype="dashed", 
  #                    color = "grey30", size = 0.3) +
  #         theme_bw() +
  #         ggtitle(country) +
  #         theme(plot.title = element_text(hjust = 0.5),
  #               legend.position = c(0.3, 0.8))
  #       
  #       list(country = country,
  #            x = item[["gdp"]],
  #            res = res,
  #            df = df,
  #            fig = fig)
  #     }
  #   )
  
  results = NULL
  for (z in 1:length(x_list)) {
    item = x_list[[z]]
    country = item[["country"]][1]
    x = item[["gdp"]]
    x_original = item[["gdp_original"]]

    y = approx(y_raw, n = (length(x) - 1) * stretch + 1, method = "linear")$y
    x = approx(x, n = (length(x) - 1) * stretch + 1, method = "linear")$y
    
    y_original = approx(y_original, n = (length(x_original) - 1) * stretch + 1, method = "linear")$y
    x_original = approx(x_original, n = (length(x_original) - 1) * stretch + 1, method = "linear")$y

    res = TwoStepDTW(y, x, x_original, t_treat, k, n_1st, dtw_method = dtw_method,
                     normalize_method = normalize_method, ...)

    df = data.frame(time = 1:length(x), y = y_original, x = x_original, warped = res$x_w[1:length(x)]) %>%
      `colnames<-`(c("time", dependent, country, paste0(country, "-Warped"))) %>%
      reshape2::melt(., id.vars = "time") %>%
      `colnames<-`(c("Time", "Country", "GDP"))

    fig = df %>%
      ggplot(aes(x = Time, y = GDP, color = Country)) +
      geom_line() +
      scale_color_manual(values = c("#2a4d69", "#ee4035", "#7bc043")) +
      geom_vline(xintercept = t_treat, linetype="dashed",
                 color = "grey30", size = 0.3) +
      theme_bw() +
      ggtitle(country) +
      theme(plot.title = element_text(hjust = 0.5),
            legend.position = c(0.3, 0.8))

    results[[country]] = list(country = country,
                              x = item[["gdp"]],
                              res = res,
                              df = df,
                              fig = fig)
  }
  
  # plot warped data
  plot_warped(results, dependent, stretch, k)
  
  # prepare data for synthetic control
  df = results %>% 
    future_map(
      ~{
        item = .
        Country = item$country
        x_w = item$res$x_w
        res = data.frame(time = 1:((length(y_raw)-1)*stretch + 1),
                         country = Country,
                         gdp_warped = NA,
                         gdp = stretch_var(Country, "gdp_original", stretch),
                         infrate = stretch_var(Country, "infrate", stretch),
                         trade = stretch_var(Country, "trade", stretch),
                         schooling = stretch_var(Country, "schooling", stretch),
                         invest60 = stretch_var(Country, "invest60", stretch),
                         invest70 = stretch_var(Country, "invest70", stretch),
                         invest80 = stretch_var(Country, "invest80", stretch),
                         industry = stretch_var(Country, "industry", stretch))
        res$gdp_warped = x_w[1:((length(y_raw)-1)*stretch + 1)]                  
        res
      }
    ) %>% 
    do.call("rbind", .) %>% 
    `row.names<-`(NULL)
  
    df = rbind(df,
             data.frame(time = 1:((length(y_raw)-1)*stretch + 1),
                        country = dependent,
                        gdp_warped = stretch_var(dependent, "gdp_original", stretch),
                        gdp = stretch_var(dependent, "gdp_original", stretch),
                        infrate = stretch_var(dependent, "infrate", stretch),
                        trade = stretch_var(dependent, "trade", stretch),
                        schooling = stretch_var(dependent, "schooling", stretch),
                        invest60 = stretch_var(dependent, "invest60", stretch),
                        invest70 = stretch_var(dependent, "invest70", stretch),
                        invest80 = stretch_var(dependent, "invest80", stretch),
                        industry = stretch_var(dependent, "industry", stretch)))
  
  df = right_join(data[c("index", "country")] %>% distinct, df, by = "country")
  df = data.frame(df)
  
  # Original method
  synth_origin = do_synth_80(df, "gdp", dependent_id, n, stretch)
  plot_synth(synth_origin, "gdp", dependent, t_treat, stretch, k)
  
  # New method
  synth_new = do_synth_80(df, "gdp_warped", dependent_id, n, stretch)
  plot_synth(synth_new, "gdp_warped", dependent, t_treat, stretch, k)
  
  # mean absolute difference: mean(abs(synthetic - dependent))
  abs_diff1 = abs(synth_origin$synthetic - synth_origin$gdp_dependent)
  abs_diff2 = abs(synth_new$synthetic - synth_new$gdp_dependent)
  
  mad1_1 = mean(abs_diff1[1:(t_treat - 1)], na.rm = T)
  mad2_1 = mean(abs_diff2[1:(t_treat - 1)], na.rm = T)
  
  mad1_2 = mean(abs_diff1[t_treat:((1987 - 1960)*stretch)], na.rm = T)
  mad2_2 = mean(abs_diff2[t_treat:((1987 - 1960)*stretch)], na.rm = T)
  
  # area between dependent and synthetic (sum of difference)
  diff1 = synth_origin$synthetic - synth_origin$gdp_dependent
  diff2 = synth_new$synthetic - synth_new$gdp_dependent
  
  sd1_1 = sum(diff1[1:(t_treat - 1)], na.rm = T)
  sd2_1 = sum(diff2[1:(t_treat - 1)], na.rm = T)
  
  sd1_2 = sum(diff1[t_treat:((1987 - 1960)*stretch)], na.rm = T)
  sd2_2 = sum(diff2[t_treat:((1987 - 1960)*stretch)], na.rm = T)
  
  return(list(
    dtw_results = results,
    df = df,
    synth_origin = synth_origin,
    synth_new = synth_new,
    mad = data.frame(mad1_pre = mad1_1,
                     mad2_pre = mad2_1,
                     mad1_post = mad1_2,
                     mad2_post = mad2_2),
    sd = data.frame(sd1_pre = sd1_1,
                    sd2_pre = sd2_1,
                    sd1_post = sd1_2,
                    sd2_post = sd2_2)
  ))
}


## Run -------------------------------------------------------------------------
countries = data[c("index", "country")] %>% distinct
k = 6
stretch = 1
result = NULL
for (i in 1:nrow(countries)) {
  dependent = countries$country[i]
  dependent_id = countries$index[i]
  print(paste0(dependent, ":", i, "-", k, " start..."))
  res = compare_methods(data = data,
                        dependent = dependent,
                        dependent_id = dependent_id,
                        stretch = stretch,
                        normalize_method = "minmax",
                        k = k,
                        step.pattern = dtw::symmetricP2)
  result = rbind(result,
                 res$mad %>% mutate(dependent = dependent, k = k))
  print(paste0(dependent, ":", i, "-", k, " start...Done."))
}
result = result %>% mutate(ratio = mad2_post/mad1_post)


# countries = data_qtrly[c("index", "country")] %>% distinct
# k = 24
# stretch = 1
# result = NULL
# for (i in 1:nrow(countries)) {
#   dependent = countries$country[i]
#   dependent_id = countries$index[i]
#   print(paste0(dependent, ":", i, "-", k, " start..."))
#   res = compare_methods(data = data_qtrly,
#                         dependent = dependent,
#                         dependent_id = dependent_id,
#                         stretch = stretch,
#                         t_treat = (81 - 1) * stretch + 1,
#                         n = (176 - 1) * stretch + 1,
#                         n_1st = (120 - 1) * stretch + 1,
#                         normalize_method = "minmax",
#                         k = k,
#                         step.pattern = dtw::symmetric2)
#   result = rbind(result,
#                  res$mad %>% mutate(dependent = dependent, k = k))
#   print(paste0(dependent, ":", i, "-", k, " start...Done."))
# }
# result = result %>% mutate(ratio = mad2_post/mad1_post)




