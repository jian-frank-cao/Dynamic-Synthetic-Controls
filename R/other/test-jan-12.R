start1 = 1971
start2 = 1975
end = 1980

t.value = data.frame(start = start1, end = start2:end, t.value = NA, p.value = NA)
wilcox = data.frame(start = start1, end = start2:end, V = NA, p.value = NA)
for (i in 1:length(start2:end)) {
  year = t.value$end[i]
  res = df.gap %>%
    filter(time %in% start1:year) %>%
    group_by(unit, data.id, grid.id) %>%
    summarise(mse.sc = mean(gap.sc^2, na.rm = T),
              mse.dsc = mean(gap.dsc^2, na.rm = T))
  res = res %>% mutate(log.ratio = log(mse.dsc/mse.sc))
  t.value$t.value[i] = t.test(res$log.ratio)$statistic
  t.value$p.value[i] = t.test(res$log.ratio)$p.value
  wilcox$V[i] = wilcox.test(res$log.ratio)$statistic
  wilcox$p.value[i] = wilcox.test(res$log.ratio)$p.value
}

rbind(data.frame(method = "sc", mse = res$mse.sc),
      data.frame(method = "dsc", mse = res$mse.dsc)) %>% 
  ggplot(aes(x = method, y = mse)) +
  geom_boxplot()

sum(res$log.ratio<0)/nrow(res)

hist(res$log.ratio, n = 100)
