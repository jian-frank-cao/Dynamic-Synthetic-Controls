my_data <- ToothGrowth
sample_n(my_data, 10)
res.ftest <- var.test(len ~ supp, data = my_data)
res.ftest

df2 = df %>% filter(time %in% 60:75)
original = df2$diff_original
new = df2$diff_new

var.test(new, original, alternative = "less")
