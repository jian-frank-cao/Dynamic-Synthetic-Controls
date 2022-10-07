target = c(8,10,14,19)

data %>% filter(id %in% target) %>% ggplot(aes(x = year, y = gdp, color = country)) +
  geom_line()

data %>% ggplot(aes(x = year, y = gdp, color = country)) +
  geom_line()

dataprep.out <-
  dataprep(
    foo = d,
    predictors    = c("gdp","trade","infrate"),
    dependent     = "gdp",
    unit.variable = 1,
    time.variable = 3,
    special.predictors = list(
      list("industry", 1971:1980, c("mean")),
      list("schooling",c(1970,1975), c("mean")),
      list("invest70" ,1980, c("mean"))
    ),
    treatment.identifier = 12,
    controls.identifier = unique(d$index)[-11],
    time.predictors.prior = 1971:1980,
    time.optimize.ssr = 1981:1990,
    unit.names.variable = 2,
    time.plot = 1960:2003
  )

# fit training model
synth.out <- 
  synth(
    data.prep.obj=dataprep.out,
    Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
  )

#----------------------------------------

# data prep for main model
dataprep.out <-
  dataprep(
    foo = d,
    predictors    = c("gdp","trade","infrate"),
    dependent     = "gdp",
    unit.variable = 1,
    time.variable = 3,
    special.predictors = list(
      list("industry" ,1981:1990, c("mean")),
      list("schooling",c(1980,1985), c("mean")),
      list("invest80" ,1980, c("mean"))
    ),
    treatment.identifier = 12,
    controls.identifier = unique(d$index)[-11],
    time.predictors.prior = 1981:1990,
    time.optimize.ssr = 1960:1989,
    unit.names.variable = 2,
    time.plot = 1960:2003
  )

# fit main model with v from training model
synth.out <- synth(
  data.prep.obj=dataprep.out,
  custom.v=as.numeric(synth.out$solution.v)
)

colnames(data)[2:4] = c("country", "year", "gdp")
gdp <- data %>% group_by(year) %>%
  summarize(OECD = mean(gdp, na.rm = TRUE)) %>%
  mutate(WGer = as.numeric(dataprep.out$Y1plot)) %>%
  mutate(Synthetic = dataprep.out$Y0%*%synth.out$solution.w %>% as.numeric) %>%
  gather(key = "country", value="gdp", OECD, WGer, Synthetic)


wgermany <- ts(gdp %>% filter(country == "WGer") %>% `$`(gdp), start=1960)
oecd <- ts(gdp %>% filter(country == "OECD") %>% `$`(gdp), start=1960)
synthetic <- ts(gdp %>% filter(country == "Synthetic") %>% `$`(gdp), start=1960)


layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
par(mai=c(0.5, 0.75, 1, 0.4))
plot.ts(wgermany,
        xlab="", ylab="",
        axes=F, xaxs="i", yaxs="i",
        xlim=c(1960, 2003), ylim=c(0,35000),
        lwd=1.5)
lines(oecd, lwd=1.5, lty=4)
abline(v=1990, lty=2, lwd=0.5)
axis(side=1, at=c(1960, 1970, 1980, 1990, 2000, 2003), lwd.ticks=0, labels=rep("", 6))
axis(side=1, lwd=0, lwd.ticks = 1)
axis(side=2)
title("(a)")
title(xlab="Year", line=2.1)
title(ylab="Per Capita GDP (PPP 2002 USD)",line=2.5)

plot.ts(wgermany,
        xlab="", ylab="",
        axes=F, xaxs="i", yaxs="i",
        xlim=c(1960, 2003), ylim=c(0,35000),
        lwd=1.5)
lines(synthetic, lwd=1.5, lty=2)
abline(v=1990, lty=2, lwd=0.5)
axis(side=1, at=c(1960, 1970, 1980, 1990, 2000, 2003), lwd.ticks=0, labels=rep("", 6))
axis(side=1, lwd=0, lwd.ticks = 1)
axis(side=2)
title("(b)")
title(xlab="Year", line=2.1)
title(ylab="Per Capita GDP (PPP 2002 USD)",line=2.5)

par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=1, 
       lty=c(1, 4, 2), 
       legend = c("West Germany", "OECD", "Synthetic West Germany"))

