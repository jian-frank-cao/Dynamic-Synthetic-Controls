
install.packages("irr")

data("anxiety", package = "irr")
head(anxiety, 4)

library(tidyverse)
library("irr")
icc(
  df, model = "twoway", 
  type = "agreement", unit = "single"
)

data = haven::read_sav("~/Downloads/adult literacy.sav")

data %>% 
  group_by(group) %>% 
  summarise(hours = mean(hours),
            postscore = mean(postscore))

df = cbind(data %>% filter(group == 0) %>% .[["postscore"]],
           c(data %>% filter(group == 1) %>% .[["postscore"]], rep(NA, 12)))

classid = unique(data$classid)
df = NULL
for (id in classid) {
  temp = data %>% filter(classid == id) %>% .[["postscore"]]
  n = length(temp)
  df = cbind(df, c(temp, rep(NA, 9-n)))
}

df = reshape2::dcast(data[c("learnid", "classid", "postscore")], learnid ~ classid, value.var = "postscore")

library(ICC)
data(ChickWeight)
ICCest(Chick, weight, data = ChickWeight, CI.type = "S")
ICCest(classid, hours, data = data, CI.type = "S")

