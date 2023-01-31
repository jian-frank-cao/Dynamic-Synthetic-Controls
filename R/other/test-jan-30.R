library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
library(furrr)
plan(multisession, workers = 8)
options(future.rng.onMisuse="ignore")
options(stringsAsFactors = FALSE)
source("./R/misc.R")
source("./R/TFDTW.R")
source("./R/synth.R")
source("./R/implement.R")
source("./R/simulate.R")
source("./R/grid.search.R")
set.seed(20220407)


x = c(0,1,2,1,0,-1,-2,-1,0)

xp = c(1,1,-1,-1,-1,-1,1,1,1)
s = 1.25
s1 = s
s2 = exp(-log(s))
phi1 = rep(0, length(x))
phi2 = rep(0, length(x))
phi1[xp > 0] = s1
phi1[xp < 0] = s2
phi2[xp > 0] = s2
phi2[xp < 0] = s1

y3 = warpWITHweight(x, phi1)[1:length(x)]
y2 = x
y1 = warpWITHweight(x, phi2)[1:length(x)]

plot(ts(y1))
lines(1:length(y1), y2, col = "green")
lines(1:length(y1), y3, col = "blue")

df = data.frame(y = y1,
                x1 = y2,
                x2 = y3,
                x1sq = y2^2,
                x2sq = y3^2)
model1 = lm(y ~ x1 + x2, data = df[1:5,])
y.pred1 = predict(model1, df[c("x1", "x2")])

model2 = lm(y ~ x1 + x2 + x1sq + x2sq, data = df[1:5,])
y.pred2 = predict(model2, df[c("x1", "x2", "x1sq", "x2sq")])

plot(ts(y1))
lines(1:length(y2), y2, col = "green")
lines(1:length(y3), y3, col = "blue")
lines(1:length(y.pred1), y.pred1, col = "orange")
lines(1:length(y.pred2), y.pred2, col = "red")

