# basque
pre.start = 7
pre.end = 16
post.start = 17
post.end = 26

# tobacco
pre.start = 11
pre.end = 20
post.start = 21
post.end = 28

# germany
pre.start = 22
pre.end = 31
post.start = 32
post.end = 41

mse = map2(
  results,
  as.list(1:length(results)),
  ~{
    item = .x
    index = .y
    mse = future_map2(
      item,
      names(item),
      ~{
        data = .x[["results.TFDTW.synth"]]
        id = .y
        mse = data %>% 
          map(
            ~{
              task = .
              unit = task$dependent
              gap.raw = task$gap.raw
              gap.TFDTW = task$gap.TFDTW
              data.frame(unit = unit,
                         mse.preT.raw = mean(gap.raw[pre.start:pre.end]^2, na.rm = T),
                         mse.preT.TFDTW = mean(gap.TFDTW[pre.start:pre.end]^2, na.rm = T),
                         mse.postT.raw = mean(gap.raw[post.start:post.end]^2, na.rm = T),
                         mse.postT.TFDTW = mean(gap.TFDTW[post.start:post.end]^2, na.rm = T)
              )
            }
          ) %>% 
          do.call("rbind", .)
        mse %>% mutate(id = id)
      }
    ) %>% 
      do.call("rbind", .)
    mse %>% mutate(data.id = index)
  }
) %>% 
  do.call("rbind", .)

filter.width.range = (1:9)*2+3
k.range = 4:9
step.pattern.range = list(
  # symmetricP0 = dtw::symmetricP0, # too bumpy
  # symmetricP05 = dtw::symmetricP05,
  symmetricP1 = dtw::symmetricP1,
  symmetricP2 = dtw::symmetricP2,
  # asymmetricP0 = dtw::asymmetricP0, # too bumpy
  # asymmetricP05 = dtw::asymmetricP05,
  asymmetricP1 = dtw::asymmetricP1,
  asymmetricP2 = dtw::asymmetricP2,
  typeIc = dtw::typeIc,
  # typeIcs = dtw::typeIcs,
  # typeIIc = dtw::typeIIc,  # jumps
  # typeIIIc = dtw::typeIIIc, # jumps
  # typeIVc = dtw::typeIVc,  # jumps
  typeId = dtw::typeId,
  # typeIds = dtw::typeIds,
  # typeIId = dtw::typeIId, # jumps
  mori2006 = dtw::mori2006
)

search.grid = expand.grid(filter.width.range, k.range,
                          names(step.pattern.range)) %>% 
  `colnames<-`(c("filter.width", "k", "step.pattern"))

search.grid = search.grid %>% 
  mutate(id = as.character(1:nrow(search.grid)),
         filter.width = as.character(filter.width),
         k = as.character(k))

mse = left_join(mse, search.grid, by = "id")

mse = mse %>% 
  mutate(log.ratio = log(mse.postT.TFDTW/mse.postT.raw),
         log.ratio.pre = log(mse.preT.TFDTW/mse.preT.raw))

model = lm(log.ratio ~ unit + filter.width + k + step.pattern, data = mse)
summary(model)

model2 = lm(log.ratio.pre ~ unit + filter.width + k + step.pattern, data = mse)
summary(model2)

res.test = data.frame(i = 1:378, t = rep(NA, 378))
for (i in 1:378) {
  res = mse %>% 
    filter(id == i)
  res.test$t[i] = t.test(res$log.ratio)$statistic
}

# --------------------------------------------------
res = mse.basque %>% 
  group_by(data.id, unit) %>% 
  top_n(-1, mse.preT.TFDTW) %>% 
  top_n(-1, id)

t.test(res$log.ratio)




units = unique(mse$unit)
res2 = NULL
for (i in units) {
  temp = mse %>% filter(unit == i)
  ind = which(temp$mse.preT.TFDTW == min(temp$mse.preT.TFDTW, na.rm = T))[1]
  res2[[i]] = temp[ind,]
}
res2 = res2 %>% do.call("rbind",.)
res = mse %>% 
  group_by(unit) %>% 
  top_n(-1, mse.preT.TFDTW) %>% 
  top_n(-1, id)
