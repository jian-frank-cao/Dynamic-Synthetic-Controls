set.seed(1)

# x - arima
x = arima.sim(list(order = c(1,1,0), ar = c(0.8)),
              n = 109)
plot(x)

# smooth
xSMA10 = ts(SMA(x, n=10)[-(1:9)])
plot(xSMA10)

# diff
xDiff = diff(xSMA10, difference = 1)
plot(xDiff)

# phi
Diff2PhiUnif = function(xDiff,
                    method = c("uniform"),
                    speed_upper = 1,
                    speed_lower = 0.5,
                    rnd = 0.3
){
  pos_deriv = xDiff >= 0
  pos_ratio = sum(pos_deriv)/length(pos_deriv)
  speed = (rnd - 0.5)*2*(speed_upper - speed_lower) +
      sign(rnd - 0.5)*speed_lower
  pos_speed = 1 + speed*(1 - pos_ratio)
  neg_speed = 1 - speed*(pos_ratio)
  phi = rep(0, length(xDiff))
  phi[pos_deriv] = pos_speed
  phi[!pos_deriv] = neg_speed
  phi = cumsum(phi)
  return(phi)
}

ApplyPhi = function(x, phi){
  output = NULL
  for (i in 1:length(phi)) {
    ind = floor(phi[i])
    weight = phi[i] - ind
    value = weight*(x[ind + 1] - x[ind]) + x[ind]
    output = c(output, value)
  }
  return(output)
}

rnd = runif(1)
phi = Diff2PhiUnif(xDiff, speed_upper = 1,
                   speed_lower = 0.5, rnd = rnd)
plot(ts(phi))
lines(1:100, col = "red")
x2 = ApplyPhi(x[-(1:9)], phi)
plot(ts(x[-(1:9)]))
lines(x2, col = "red")

# decompose
x <- ts(x, frequency = 20, start = c(1946, 1))

xcomp = decompose(x)
plot(xcomp)
