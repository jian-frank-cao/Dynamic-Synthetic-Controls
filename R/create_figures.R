## 1. Speed Difference ---------------------------------------------------------
set.seed(4) # 3

sim.data = function(n = 3, nCycles = 6, length = 100, extra.x = round(0.2*length),
                    t.treat = 60, shock = 0, ar.x = 0.6, n.SMA = 1,
                    n.diff = 1, speed.upper = 2, speed.lower = 0.5,
                    reweight = TRUE, rescale = TRUE, 
                    rescale.multiplier = 20, beta = 1){
  # common exogenous shocks
  x = sin(seq(0, nCycles * pi, length.out = length + extra.x + n.SMA + n.diff - 1))
  
  # smoothing
  x.SMA = ts(TTR::SMA(x, n = n.SMA)[-(1:(n.SMA - 1))])
  
  # difference
  x.diff = diff(x.SMA, difference = n.diff)
  pos.diff = x.diff > 0
  if (reweight) {
    pos.ratio = sum(pos.diff)/sum(!pos.diff)
  }
  
  # speeds
  log.speeds = seq(log(speed.lower), log(speed.upper), length.out = n)
  rnd.ind = 2
  log.speeds = c(log.speeds[rnd.ind], log.speeds[-rnd.ind])
  
  # simulate
  data = NULL
  for (i in 1:n) {
    # speed profile
    log.speed = log.speeds[i]
    if (reweight) {
      if (pos.ratio > 1) {
        pos.speed = exp(log.speed*(1/pos.ratio))
        neg.speed = exp(-log.speed)
      }else{
        pos.speed = exp(log.speed)
        neg.speed = exp(-log.speed*pos.ratio)
      }
    }else{
      pos.speed = exp(log.speed)
      neg.speed = exp(-log.speed)
    }
    
    phi.shape = rep(NA, length.out = length + extra.x)
    phi.shape[pos.diff] = pos.speed
    phi.shape[!pos.diff] = neg.speed
    
    log.phi.mean = mean(log(phi.shape), na.rm = T)
    log.phi.sd = sd(log(phi.shape), na.rm = T)
    
    phi.random = exp(rnorm(n = length + extra.x,
                           mean = log.phi.mean,
                           sd = log.phi.sd))
    
    # treatment
    if (i == 1) {
      treatment = c(rep(0, t.treat),
                    seq(0, shock, length.out = round(0.1*length)),
                    rep(shock, round(0.9*length - t.treat)))
    }else{
      treatment = 0
    }
    
    phi = beta*phi.shape + (1 - beta)*phi.random
    
    y = warpWITHweight(x[1:(length + extra.x)], phi)[1:length]
    
    y = y + rnorm(n = length(y), mean = 0, sd = 0.02)
    
    if (rescale) {
      y = minmax.normalize(y, reference = y[1:t.treat])*rescale.multiplier
    }
    
    y = y + treatment
    
    # y = cumsum(y)
    
    data = rbind(data,
                 data.frame(id = i,
                            unit = LETTERS[i],
                            time = 1:length,
                            value = y,
                            value_raw = y))
  }
  return(data)
}

data = sim.data(n = 3, nCycles = 6, length = 100,
                t.treat = 60, shock = 0, ar.x = 0.6,
                n.SMA = 1, n.diff = 1,
                speed.upper = 2,
                speed.lower = 0.5,
                reweight = TRUE,
                rescale = TRUE,
                rescale.multiplier = 20,
                beta = 1)

data %>% ggplot(aes(x = time, y = value, color = unit)) + geom_line()

## 3. Placebo Test -------------------------------------------------------------
fig_all = ggpubr::ggarrange(fig_basque + ggpubr::rremove("ylab") +
                              ggpubr::rremove("xlab"),
                            fig_tobacco + ggpubr::rremove("ylab") +
                              ggpubr::rremove("xlab"),
                            fig_germany + ggpubr::rremove("ylab") +
                              ggpubr::rremove("xlab"),
                            labels = c("Basque Terrorism", 
                                       "California Tobacco",
                                       "German Reunification"),
                            label.x = 0.1,
                            label.y = 0.9,
                            hjust = c(0, 0, 0),
                            font.label = list(size = 18, color = "grey20",
                                              face = "bold"),
                            ncol = 1, nrow = 3,
                            # common.legend = TRUE, legend = "bottom",
                            align = "hv")

fig_all = ggpubr::annotate_figure(fig_all,
                        left = grid::textGrob("Gap (y - Synthetic Control)",
                                        rot = 90, vjust = 1, 
                                        gp = grid::gpar(cex = 1.3)),
                bottom = grid::textGrob("Time", 
                                  gp = grid::gpar(cex = 1.3)))

ggsave("./figures/placebo_all_1108.pdf",
       fig_all, width = 8, height = 8,
       units = "in", limitsize = FALSE)
