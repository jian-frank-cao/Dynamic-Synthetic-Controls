library(psych)
library(reshape)
library(tidyverse)

sf <- matrix(c(
  9,    2,   5,    8,
  6,    1,   3,    2,
  8,    4,   6,    8,
  7,    1,   2,    6,
  10,   5,   6,    9,
  6,   2,   4,    7),ncol=4,byrow=TRUE)
colnames(sf) <- paste("J",1:4,sep="")
rownames(sf) <- paste("S",1:6,sep="")
ICC(sf,lmer=FALSE)

sf.df <- melt(sf, varnames=c("Subject", "Rater"))
sf.df2 = rbind(sf.df %>% mutate(id = 1),
               sf.df %>% mutate(value = value + rnorm(nrow(sf)), id = 2),
               sf.df %>% mutate(value = value + 2*rnorm(nrow(sf)), id = 3),
               sf.df %>% mutate(value = value + 2*rnorm(nrow(sf)), id = 4),
               sf.df %>% mutate(value = value + 2*rnorm(nrow(sf)), id = 5),
               sf.df %>% mutate(value = value + 2*rnorm(nrow(sf)), id = 6),
               sf.df %>% mutate(value = value + 2*rnorm(nrow(sf)), id = 7))
sf.df2 = sf.df2 %>% mutate(unit = paste0(Subject, "-", Rater))
res.aov = aov(value ~ Subject*Rater, sf.df2)
res.aov = aov(value ~ unit, sf.df2)
summary(res.aov)

model = lm(value ~ Subject*Rater, data = sf.df2)
model = lm(value ~ unit, data = sf.df2)
regclass::VIF(model)

sf.df3 = reshape2::dcast(sf.df2[c("id", "unit", "value")],
                              id ~ unit, value.var = "value")
res.icc = irr::icc(
  sf.df3[,-1], model = "oneway",
  type = "c", unit = "single"
)$value
(nrow(sf.df3) - 1)*res.icc + 1

res.aov = aov(gap_original ~ unit*data.id, df)
summary(res.aov)

df$group = paste0(df$unit, "-", df$data.id)
df0 = reshape2::dcast(df[c("group", "time", "value")],
                         time ~ group, value.var = "value")

units = unique(df$unit)
df0 = rep(0,304)
for (item in units) {
  temp = df %>% filter(unit == item) %>% select(gap_original)
  df0 = cbind(df0, temp)
}
res.icc = irr::icc(
  df0[,-1], model = "twoway",
  type = "a", unit = "single"
)$value
(nrow(df0) - 1)*res.icc + 1

MST = 34728
MSJ = 28
MSI = 26
MSE = 34
L = 10
K = 39
N = 38

target = ((MST-MSI)/(L*K) + (MSJ-MSI)/(L*N) + (MSI-MSE)/(L)) 
target/(target + MSE)
