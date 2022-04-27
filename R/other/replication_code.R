## Note! You need to install the checkpoint package first using

###############################################################
             install.packages("checkpoint")                
###############################################################


library(checkpoint)
checkpoint("2020-05-19")

library(foreign)
library(Synth)
library(tidyverse)
library(haven)
library(xtable)


# +----------+
# | Figure 1 |
# +----------+

d <- read.dta("./data/repgermany.dta")
data <- read_dta("./data/repgermany.dta") %>%
    mutate_at(vars(year, gdp, infrate, trade, schooling,
                   invest60, invest70, invest80,
                   industry),
               as.numeric) %>%
    mutate_at(vars(index, country), as.factor)


## pick v by cross-validation
# data setup for training model
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
           treatment.identifier = 7,
           controls.identifier = unique(d$index)[-7],
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
  treatment.identifier = 7,
  controls.identifier = unique(d$index)[-7],
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


gdp <- data %>% group_by(year) %>%
    summarize(OECD = mean(gdp, na.rm = TRUE)) %>%
    mutate(WGer = as.numeric(dataprep.out$Y1plot)) %>%
    mutate(Synthetic = dataprep.out$Y0%*%synth.out$solution.w %>% as.numeric) %>%
    gather(key = "country", value="gdp", OECD, WGer, Synthetic)


wgermany <- ts(gdp %>% filter(country == "WGer") %>% `$`(gdp), start=1960)
oecd <- ts(gdp %>% filter(country == "OECD") %>% `$`(gdp), start=1960)
synthetic <- ts(gdp %>% filter(country == "Synthetic") %>% `$`(gdp), start=1960)


pdf("./figures/figure1.pdf",width=10,height=6,paper='special') 
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

graphics.off()



# +----------+
# | Figure 3 |
# +----------+
###  Placebo Reunification 1975 - Trends in Per-Capita GDP: West Germany vs.
###  Synthetic West Germany


d <- read.dta("repgermany.dta")


# data prep for training model
dataprep.out <-
  dataprep(
    foo = d,
    predictors    = c("gdp","trade","infrate"),
    dependent     = "gdp",
    unit.variable = 1,
    time.variable = 3,
    special.predictors = list(
      list("industry",1971, c("mean")),
      list("schooling",c(1960,1965), c("mean")),
      list("invest60" ,1980, c("mean"))
    ),
    treatment.identifier = 7,
    controls.identifier = unique(d$index)[-7],
    time.predictors.prior = 1960:1970,
    time.optimize.ssr = 1970:1980,
    unit.names.variable = 2,
    time.plot = 1960:1990
  )

# fit training model
synth.out <- synth(
  data.prep.obj=dataprep.out,
  Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
)


# data prep for main model
dataprep.out <-
  dataprep(
    foo = d,
    predictors    = c("gdp","trade","infrate"),
    dependent     = "gdp",
    unit.variable = 1,
    time.variable = 3,
    special.predictors = list(
      list("industry" ,1971:1975, c("mean")),
      list("schooling",c(1965,1980), c("mean")),
      list("invest80" ,1980, c("mean"))
    ),
    treatment.identifier = 7,
    controls.identifier = unique(d$index)[-7],
    time.predictors.prior = 1965:1980,
    time.optimize.ssr = 1960:1980,
    unit.names.variable = 2,
    time.plot = 1960:2003
  )

# fit main model
synth.out <- synth(
  data.prep.obj=dataprep.out,
  custom.v=as.numeric(synth.out$solution.v)
)


par(mai=c(0.5, 0.75, 1, 0.4))
pdf("figure3.pdf",width=6,height=5.5,paper='special') 

plot(1960:2003,dataprep.out$Y1plot,
        xlab="", ylab="",
        axes=F, xaxs="i", yaxs="i",
        xlim=c(1960, 2003), ylim=c(0,35000),
        lwd=1.5, lty=1, type='l'
     )
axis(side=1, at=c(1960, 1970, 1980, 1990, 2000, 2003), lwd.ticks=0, labels=rep("", 6))
axis(side=1, lwd=0, lwd.ticks = 1)
axis(side=2)

title(xlab="Year", line=2.1)
title(ylab="Per Capita GDP (PPP 2002 USD)",line=2.5)



lines(1960:2003,(dataprep.out$Y0%*%synth.out$solution.w),col="black",lty="dashed",lwd=2)
abline(v=1980, lty=2, lwd=0.5)
legend(x="bottomright",
       legend=c("West Germany","Synthetic West Germany")
      ,lty=c("solid","dashed"),col=c("black","black")
      ,cex=.8,bg="white",lwd=c(2,2))
arrows(1974,18000,1979.5,18000,col="black",length=.1)
text(1963,20000,"Reunification backdated", cex=1, adj = c(0, NA))
text(1963,18000,"to 1980", cex=1, adj = c(0, NA))

graphics.off()




# +----------+
# | Figure 4 |
# +----------+

d <- read.dta("repgermany.dta")

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
           treatment.identifier = 7,
           controls.identifier = unique(d$index)[-7],
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
  treatment.identifier = 7,
  controls.identifier = unique(d$index)[-7],
  time.predictors.prior = 1981:1990,
  time.optimize.ssr = 1960:1989,
  unit.names.variable = 2,
  time.plot = 1960:2003
)

synth.out <- synth(
  data.prep.obj=dataprep.out,
  custom.v=as.numeric(synth.out$solution.v)
  )

synthY0 <- (dataprep.out$Y0%*%synth.out$solution.w)


storegaps <- 
  matrix(NA,
        length(1960:2003),
        5)
colnames(storegaps) <- c(1,3,9,12,14)
co <- unique(d$index)[-7]

for(k in 1:5){

# data prep for training model
omit <- c(1,3,9,12,14)[k]  
  dataprep.out <-
    dataprep(
      foo = d,
      predictors    = c("gdp","trade","infrate"),
      dependent     = "gdp",
      unit.variable = 1,
      time.variable = 3,
      special.predictors = list(
        list("industry",1971:1980, c("mean")),
        list("schooling"   ,c(1970,1975), c("mean")),
        list("invest70" ,1980, c("mean"))
      ),
      treatment.identifier = 7,
      controls.identifier = co[-which(co==omit)],
      time.predictors.prior = 1971:1980,
      time.optimize.ssr = 1981:1990,
      unit.names.variable = 2,
      time.plot = 1960:2003
    )
  
  # fit training model
  synth.out <- synth(
    data.prep.obj=dataprep.out,
    Margin.ipop=.005,Sigf.ipop=7,Bound.ipop=6
  )
  
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
    treatment.identifier = 7,
    controls.identifier = co[-which(co==omit)],
    time.predictors.prior = 1981:1990,
    time.optimize.ssr = 1960:1989,
    unit.names.variable = 2,
    time.plot = 1960:2003
  )
  
  # fit main model 
  synth.out <- synth(
    data.prep.obj=dataprep.out,
    custom.v=as.numeric(synth.out$solution.v)
  )
  storegaps[,k] <- (dataprep.out$Y0%*%synth.out$solution.w)
} # close loop over leave one outs

pdf(file = "figure4.pdf", width = 6, height = 5.5, paper='special')
plot(1960:2003,dataprep.out$Y1plot,
        xlab="", ylab="",
        axes=F, xaxs="i", yaxs="i",
        xlim=c(1960, 2003), ylim=c(0,35000),
        lwd=1.5, lty=1, type='l'
     )

axis(side=1, at=c(1960, 1970, 1980, 1990, 2000, 2003), lwd.ticks=0, labels=rep("", 6))
axis(side=1, lwd=0, lwd.ticks = 1)
axis(side=2)

title(xlab="Year", line=2.1)
title(ylab="Per Capita GDP (PPP 2002 USD)",line=2.5)

abline(v=1990, lty=2, lwd=0.5)
for(i in 1:5) {
  lines(1960:2003,storegaps[,i],col="darkgrey",lty="solid")
}
lines(1960:2003,synthY0,col="black",lty="dashed",lwd=2)
lines(1960:2003,dataprep.out$Y1plot,col="black",lty="solid",lwd=2)

legend(x="bottomright",
       legend=c("West Germany",
                "synthetic West Germany",
                "synthetic West Germany (leave-one-out)")
      ,lty=c("solid","dashed","solid"),
      col=c("black","black","darkgrey")
      ,cex=.8,bg="white",lwd=c(2,2,1))

graphics.off()


# +--------+
# | Tables |
# +--------+

d <- read.dta("repgermany.dta")

# ## Table 1 & 2, Figure 1, 2, & 3

## pick v by cross-validation
# data setup for training model
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
           treatment.identifier = 7,
           controls.identifier = unique(d$index)[-7],
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
  treatment.identifier = 7,
  controls.identifier = unique(d$index)[-7],
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

#### Table 1
synth.tables <- synth.tab(
                          dataprep.res = dataprep.out,
                          synth.res = synth.out
                          ); synth.tables

colnames(synth.tables$tab.pred)[3] <- "Rest of OECD Sample"
rownames(synth.tables$tab.pred) <- c("GDP per-capita","Trade openness",
                                     "Inflation rate","Industry share",
                                     "Schooling","Investment rate")

table1 <- as_tibble(synth.tables$tab.pred, rownames = "characteristic")

# Find the closest country
data <- read_dta("repgermany.dta") %>%
    mutate_at(vars(year, gdp, infrate, trade, schooling, invest60, invest70, invest80, industry), as.numeric) %>%
    mutate_at(vars(index, country), as.factor) %>%
    mutate(westger = country == "West Germany")
piece.1 <- data %>% filter(1981 <= year, year <= 1990) %>%
    gather(key="characteristic", value="value", gdp, trade, industry, infrate) %>%
    group_by(country, characteristic) %>% summarize(m = mean(value, na.rm = T))

piece.2 <- data %>% 
    filter(1980 <= year, year <= 1985) %>%
    gather(key="characteristic", value="value", invest80, schooling) %>%
    group_by(country, characteristic) %>% summarize(m = mean(value, na.rm = T))

piece <- bind_rows(piece.1, piece.2) %>%
    mutate(country.type = factor(country == "West Germany", levels = c(TRUE, FALSE), labels = c("West Germany", "OECD"))) %>%
    ungroup()
piece %>% group_by(country.type, characteristic) %>% summarize(m = mean(m)) %>% spread(country.type, m)


sigma <- piece %>% group_by(characteristic) %>% summarize(sigma = sd(m))
standardized <- inner_join(piece, sigma) %>% ungroup() %>% mutate(m = m / sigma) %>% select(-sigma, -country.type)
standardized %>% group_by(country) %>% spread(characteristic, m)

norm <- standardized %>%
    filter(country == "West Germany") %>%
    select(-country) %>%
    rename(west.germany = m) %>%
    inner_join(standardized) %>%
    filter(country != "West Germany") %>%
    mutate(m = m - west.germany) %>% select(-west.germany) %>%
    group_by(country) %>% summarize(norm=sqrt(sum(m * m))) %>% arrange(norm)

closest_country <- norm[[1, 1]] %>% as.character
closest_chars <- piece %>% filter(country == closest_country) %>% select(characteristic, m)
names(closest_chars)[2] <- closest_country
closest_chars$characteristic <- c("GDP per-capita","Industry share","Inflation rate",
                                "Trade openness","Investment rate", "Schooling")

write_csv(full_join(table1, closest_chars), "table1.csv")




#### Table 2 and 3

# synth weights
tab1 <- data.frame(synth.tables$tab.w)
# regression weights
X0 <- cbind(1,t(dataprep.out$X0))
X1 <- as.matrix(c(1,dataprep.out$X1))
W     <- X0%*%solve(t(X0)%*%X0)%*%X1
Wdat  <- data.frame(unit.numbers=as.numeric(rownames(X0)),
                    regression.w=W)
tab1  <- merge(tab1,Wdat,by="unit.numbers")

table23 <- as_tibble(tab1) %>%
    select(country = unit.names,
           synthetic.control.weights = w.weights,
           regression.weights = regression.w) %>% 
    mutate(country = as.factor(country)) %>%
    mutate_if(is.numeric, ~ round(., 2))

write_csv(table23 %>% select(country, synthetic.control.weights), "table2.csv")
write_csv(table23 %>% select(country, regression.weights), "table3.csv")

