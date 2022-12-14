Input =("
 Tech  Rat Protein
 Janet 1   1.119
 Janet 1   1.2996
 Janet 1   1.5407
 Janet 1   1.5084
 Janet 1   1.6181
 Janet 1   1.5962
 Janet 1   1.2617
 Janet 1   1.2288
 Janet 1   1.3471
 Janet 1   1.0206
 Janet 2   1.045
 Janet 2   1.1418
 Janet 2   1.2569
 Janet 2   0.6191
 Janet 2   1.4823
 Janet 2   0.8991
 Janet 2   0.8365
 Janet 2   1.2898
 Janet 2   1.1821
 Janet 2   0.9177
 Janet 3   0.9873
 Janet 3   0.9873
 Janet 3   0.8714
 Janet 3   0.9452
 Janet 3   1.1186
 Janet 3   1.2909
 Janet 3   1.1502
 Janet 3   1.1635
 Janet 3   1.151
 Janet 3   0.9367
 Brad  5   1.3883
 Brad  5   1.104
 Brad  5   1.1581
 Brad  5   1.319
 Brad  5   1.1803
 Brad  5   0.8738
 Brad  5   1.387
 Brad  5   1.301
 Brad  5   1.3925
 Brad  5   1.0832
 Brad  6   1.3952
 Brad  6   0.9714
 Brad  6   1.3972
 Brad  6   1.5369
 Brad  6   1.3727
 Brad  6   1.2909
 Brad  6   1.1874
 Brad  6   1.1374
 Brad  6   1.0647
 Brad  6   0.9486
 Brad  7   1.2574
 Brad  7   1.0295
 Brad  7   1.1941
 Brad  7   1.0759
 Brad  7   1.3249
 Brad  7   0.9494
 Brad  7   1.1041
 Brad  7   1.1575
 Brad  7   1.294
 Brad  7   1.4543
")

Data = read.table(textConnection(Input),header=TRUE)

### Since Rat is read in as an integer variable, convert it to factor

Data$Rat = as.factor(Data$Rat)

library(nlme)

fit = aov(Protein ~ Tech + Error(Rat), data=Data)
summary(fit)

df = df %>% 
  mutate(data.id = substr(unit, 1, 3))

fit = aov(gap_original ~ data.id + Error(unit), data=df)
summary(fit)

pf(q=939.8/33.92,
   df1=1405,
   df2=10108,
   lower.tail=F)


npk.aov <- aov(yield ~ block + N*P*K, npk)
summary(npk.aov)
coefficients(npk.aov)

aov(yield ~ block + N * P + K, npk)
aov(terms(yield ~ block + N * P + K, keep.order = TRUE), npk)

npk.aovE <- aov(yield ~  N*P*K + Error(block), npk)
npk.aovE
summary(npk.aovE)


df = df %>%
  mutate(
    data.id = str_split(unit, "-", simplify = TRUE)[,1],
    unit = str_split(unit, "-", simplify = TRUE)[,2]
  )

res.aov.sc = aov(gap_original ~ data.id*unit, df)
summary.aov.sc = summary(res.aov.sc)

DF.sc = summary.aov.sc[[1]]$Df %>% sum
sum.sq.sc = summary.aov.sc[[1]]$`Sum Sq` %>% sum
var.sc = sum.sq.sc/DF.sc


res.aov.dsc = aov(gap_new ~ data.id*unit, df)
summary.aov.dsc = summary(res.aov.dsc)

DF.dsc = summary.aov.dsc[[1]]$Df %>% sum
sum.sq.dsc = summary.aov.dsc[[1]]$`Sum Sq` %>% sum
var.dsc = sum.sq.dsc/DF.dsc

pf(var.dsc/var.sc, DF.dsc, DF.sc, lower.tail = TRUE)*2

# -------------------------------------------
# res = df %>%
#   group_by(id) %>%
#   summarise(mse.sc = mean(diff_original^2, na.rm = T),
#             mse.dsc = mean(diff_new^2, na.rm = T))
# res = res %>%
#   mutate(log.ratio = log(mse.dsc/mse.sc))
# t.test(res$log.ratio)

# -------------------------------------------


# df = df %>% filter(time %in% 1990:1997)
res = df.basque %>% 
  group_by(unit) %>% 
  summarise(mse.sc = mean(gap_original^2, na.rm = T),
            mse.dsc = mean(gap_new^2, na.rm = T))

res = res %>% 
  mutate(
    data.id = str_split(unit, "-", simplify = TRUE)[,1],
    unit = str_split(unit, "-", simplify = TRUE)[,2]
  )

res = res %>% 
  mutate(log.ratio = log(mse.dsc/mse.sc))
t.test(res$log.ratio)

mse.aov.sc = aov(mse.sc ~ data.id*unit, res)
summary.aov.sc = summary(mse.aov.sc)

DF.sc = summary.aov.sc[[1]]$Df %>% sum
sum.sq.sc = summary.aov.sc[[1]]$`Sum Sq` %>% sum
var.sc = sum.sq.sc/DF.sc


mse.aov.dsc = aov(mse.dsc ~ data.id*unit, res)
summary.aov.dsc = summary(mse.aov.dsc)

DF.dsc = summary.aov.dsc[[1]]$Df %>% sum
sum.sq.dsc = summary.aov.dsc[[1]]$`Sum Sq` %>% sum
var.dsc = sum.sq.dsc/DF.dsc

pf(var.dsc/var.sc, DF.dsc, DF.sc, lower.tail = TRUE)*2

rbind(data.frame(value = res$mse.sc,
                 method = "sc"),
      data.frame(value = res$mse.dsc,
                 method = "dsc")) %>% 
  ggplot(aes(x = method, y = value)) +
  geom_boxplot()
