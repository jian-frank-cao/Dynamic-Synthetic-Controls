## Setup -----------------------------------------------------------------------
library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)


## Data ------------------------------------------------------------------------
symbols = c("TIP", "INXG.L", "XRB.TO", "IBCI.DE") #, "GTIP" , "ILB.XA", "XSTH.TO"

data_list = NULL

# Loop over the symbols and download the data
for(symbol in symbols) {
  response = quantmod::getSymbols(symbol, auto.assign = FALSE,
                             from = "2010-02-01", to = "2022-12-31")
  
  data_list[[symbol]] = response %>% 
    data.frame(.) %>% 
    mutate(date = as.Date(rownames(.))) %>% 
    `rownames<-`(NULL) %>% 
    .[,c(7:6)] %>% 
    `colnames<-`(c("date", symbol))
}

data = Reduce(function(...) merge(..., by = "date", all = TRUE), data_list)
cols_to_fill <- names(data)[-1]
data = data %>% fill(all_of(cols_to_fill))

## Functions -------------------------------------------------------------------
# normalization
t.normalize = function(data, reference = NULL){
  if (is.null(reference)) {
    mu = mean(data)
    sigma = sd(data)
  }else{
    mu = mean(reference)
    sigma = sd(reference)
  }
  res = (data - mu)/sigma
  return(res)
}

# transform warping path W to weight
warp2weight = function(W){
  w = as.matrix(W)
  count = rep(1/colSums(w), nrow(w)) %>% 
    matrix(.,
           nrow = ncol(w),
           ncol = nrow(w)) %>% 
    t(.)
  weight = rowSums(w * count)
  
  return(weight)
}


## Run -------------------------------------------------------------------------
df = reshape2::melt(data, id.vars = 1)

# plot
df %>%
  ggplot(aes(x = date, y = value, color = variable)) +
  geom_line()

# dtw
df_a = data$TIP
df_b = data$INXG.L

align = dtw::dtw(diff(t.normalize(df_a)),
                 diff(t.normalize(df_b)),
                 keep = TRUE,
                 step.pattern = dtw::symmetricP2)
dtw::dtwPlotThreeWay(align)
P = Matrix::sparseMatrix(align$index1,
                         align$index2)
W = warp2weight(P)
a = fitted(forecast::ets(W))
plot(ts(a))



