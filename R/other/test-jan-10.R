mse = mse %>% 
  mutate(log.ratio.pre = log(mse.preT.TFDTW/mse.preT.raw))

search.grid = search.grid %>% 
  mutate(id = as.character(1:nrow(search.grid)),
         filter.width = as.character(filter.width),
         k = as.character(k))

mse = left_join(mse, search.grid, by = "id")

model = lm(log.ratio ~ unit + filter.width + k + step.pattern, data = mse)
summary(model)


model2 = lm(log.ratio.pre ~ unit + filter.width + k + step.pattern, data = mse)
summary(model2)
