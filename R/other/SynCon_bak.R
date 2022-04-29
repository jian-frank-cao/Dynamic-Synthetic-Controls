library(foreign)
library(Synth)
library(tidyverse)
# library(haven)
library(xtable)

## Functions -------------------------------------------------------------------
do_synth_80 = function(df, dep_var, dependent_id, n, stretch){
  # find v
  dataprep.out <-
    dataprep(
      foo = df,
      predictors    = c(dep_var,"trade","infrate"),
      dependent     = dep_var,
      unit.variable = 1,
      time.variable = 3,
      special.predictors = list(
        list("industry", ((1971-1960)*stretch+1):((1973-1960)*stretch+1), c("mean")),
        list("schooling",c(((1965-1960)*stretch+1),((1970-1960)*stretch+1)), c("mean")),
        list("invest60" ,((1980-1960)*stretch+1), c("mean"))
      ),
      treatment.identifier = dependent_id,
      controls.identifier = setdiff(unique(df$index), dependent_id),
      time.predictors.prior = ((1963-1960)*stretch+1):((1973-1960)*stretch),
      time.optimize.ssr = ((1973-1960)*stretch+1):((1980-1960)*stretch), 
      unit.names.variable = 2,
      time.plot = 1:n
    )
  
  # fit training model
  synth.out <- 
    synth(
      data.prep.obj=dataprep.out,
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  
  # main model
  # prepare data
  dataprep.out <-
    dataprep(
      foo = df,
      predictors    = c(dep_var,"trade","infrate"),
      dependent     = dep_var,
      unit.variable = 1,
      time.variable = 3,
      special.predictors = list(
        list("industry", ((1973-1960)*stretch+1):((1980-1960)*stretch+1), c("mean")),
        list("schooling",c(((1975-1960)*stretch+1)), c("mean")),
        list("invest70" ,((1980-1960)*stretch+1), c("mean"))
      ),
      treatment.identifier = dependent_id,
      controls.identifier = setdiff(unique(df$index), dependent_id),
      time.predictors.prior = ((1973-1960)*stretch+1):((1980-1960)*stretch+1),
      time.optimize.ssr = 1:((1980-1960)*stretch),
      unit.names.variable = 2,
      time.plot = 1:n
    )
  
  # fit training model
  synth.out <- 
    synth(
      data.prep.obj=dataprep.out,
      custom.v=as.numeric(synth.out$solution.v)
    )
  
  gdp <- df %>% group_by(time) %>%
    summarize(OECD = mean(gdp, na.rm = TRUE)) %>%
    mutate(target = as.numeric(dataprep.out$Y1plot)) %>%
    mutate(Synthetic = dataprep.out$Y0%*%synth.out$solution.w %>% as.numeric) %>%
    gather(key = "country", value="gdp", OECD, target, Synthetic)
  
  gdp_dependent <- ts(gdp %>% filter(country == "target") %>% `$`(gdp), start=1)
  oecd <- ts(gdp %>% filter(country == "OECD") %>% `$`(gdp), start=1)
  synthetic <- ts(gdp %>% filter(country == "Synthetic") %>% `$`(gdp), start=1)
  
  return(list(gdp_dependent = gdp_dependent,
              oecd = oecd,
              synthetic = synthetic))
}


do_synth_80_qtrly = function(df, dep_var, dependent_id, n, stretch){
  # find v
  dataprep.out <-
    dataprep(
      foo = df,
      predictors    = c(dep_var,"trade","infrate"),
      dependent     = dep_var,
      unit.variable = 1,
      time.variable = 3,
      special.predictors = list(
        list("industry", ((1971-1960)*4*stretch+1):((1973-1960)*4*stretch+1), c("mean")),
        list("schooling",c(((1965-1960)*4*stretch+1),((1970-1960)*4*stretch+1)), c("mean")),
        list("invest60" ,((1980-1960)*4*stretch+1), c("mean"))
      ),
      treatment.identifier = dependent_id,
      controls.identifier = setdiff(unique(df$index), dependent_id),
      time.predictors.prior = ((1963-1960)*4*stretch+1):((1973-1960)*4*stretch),
      time.optimize.ssr = ((1973-1960)*4*stretch+1):((1980-1960)*4*stretch), 
      unit.names.variable = 2,
      time.plot = 1:n
    )
  
  # fit training model
  synth.out <- 
    synth(
      data.prep.obj=dataprep.out,
      Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
    )
  
  # main model
  # prepare data
  dataprep.out <-
    dataprep(
      foo = df,
      predictors    = c(dep_var,"trade","infrate"),
      dependent     = dep_var,
      unit.variable = 1,
      time.variable = 3,
      special.predictors = list(
        list("industry", ((1973-1960)*4*stretch+1):((1980-1960)*4*stretch+1), c("mean")),
        list("schooling",c(((1975-1960)*4*stretch+1)), c("mean")),
        list("invest70" ,((1980-1960)*4*stretch+1), c("mean"))
      ),
      treatment.identifier = dependent_id,
      controls.identifier = setdiff(unique(df$index), dependent_id),
      time.predictors.prior = ((1973-1960)*4*stretch+1):((1980-1960)*4*stretch+1),
      time.optimize.ssr = 1:((1980-1960)*4*stretch),
      unit.names.variable = 2,
      time.plot = 1:n
    )
  
  # fit training model
  synth.out <- 
    synth(
      data.prep.obj=dataprep.out,
      custom.v=as.numeric(synth.out$solution.v)
    )
  
  gdp <- df %>% group_by(time) %>%
    summarize(OECD = mean(gdp, na.rm = TRUE)) %>%
    mutate(target = as.numeric(dataprep.out$Y1plot)) %>%
    mutate(Synthetic = dataprep.out$Y0%*%synth.out$solution.w %>% as.numeric) %>%
    gather(key = "country", value="gdp", OECD, target, Synthetic)
  
  gdp_dependent <- ts(gdp %>% filter(country == "target") %>% `$`(gdp), start=1)
  oecd <- ts(gdp %>% filter(country == "OECD") %>% `$`(gdp), start=1)
  synthetic <- ts(gdp %>% filter(country == "Synthetic") %>% `$`(gdp), start=1)
  
  return(list(gdp_dependent = gdp_dependent,
              oecd = oecd,
              synthetic = synthetic))
}


plot_synth = function(df, dep_var, dependent, t_treat, stretch, k){
  n = length(df$gdp_dependent)
  pdf(paste0("./figures/synth_control_",
             paste0(c(dependent, dep_var, stretch, k), collapse = "_"),
             ".pdf"),
      width=10,height=6,paper='special') 
  layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
  par(mai=c(0.5, 0.75, 1, 0.4))
  plot.ts(df$gdp_dependent,
          xlab="", ylab="",
          axes=F, xaxs="i", yaxs="i",
          xlim=c(1, n), ylim=c(0,35000),
          lwd=1.5)
  lines(df$oecd, lwd=1.5, lty=4)
  abline(v=t_treat, lty=2, lwd=0.5)
  axis(side=1, at=seq(1, n, length.out = 5), lwd.ticks=0, labels=rep("", 5))
  axis(side=1, lwd=0, lwd.ticks = 1)
  axis(side=2)
  title("(a)")
  title(xlab="Year", line=2.1)
  title(ylab="Per Capita GDP (PPP 2002 USD)",line=2.5)
  
  plot.ts(df$gdp_dependent,
          xlab="", ylab="",
          axes=F, xaxs="i", yaxs="i",
          xlim=c(1, n), ylim=c(0,35000),
          lwd=1.5)
  lines(df$synthetic, lwd=1.5, lty=2)
  abline(v=t_treat, lty=2, lwd=0.5)
  axis(side=1, at=seq(1, n, length.out = 5), lwd.ticks=0, labels=rep("", 5))
  axis(side=1, lwd=0, lwd.ticks = 1)
  axis(side=2)
  title("(b)")
  title(xlab="Year", line=2.1)
  title(ylab="Per Capita GDP (PPP 2002 USD)",line=2.5)
  
  par(mai=c(0,0,0,0))
  plot.new()
  legend(x="center", ncol=1, 
         lty=c(1, 4, 2), 
         legend = c(dependent, "OECD", paste0("Synthetic ", dependent)))
  
  graphics.off()
}

