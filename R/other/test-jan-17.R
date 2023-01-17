library(checkpoint)
checkpoint("2022-04-01")

library(tidyverse)
Wa = matrix(rep(0, 5*5), ncol = 5)

Wb[1,1] = 1
Wb[2,2:3] = 1
Wa[3:4,3] = 1
Wa[4,4:5] = 1
Wa[5,5] = 1

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

weight.a = warp2weight(Wa)


Wb = matrix(rep(0, 3*3), ncol = 3)

Wb[1,1:2] = 1
Wb[2:3,2] = 1
Wb[3,3] = 1

W.pp.i = Wb
weight.a.Rs = weight.a[2:4]
weight.b = (W.pp.i %*% weight.a.Rs)/rowSums(as.matrix(W.pp.i))
weight.b = as.numeric(weight.b)

res = Wa[,2:4] %*% W.pp.i
warp2weight(res)
