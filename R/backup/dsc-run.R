# Load necessary packages
library(dsc)

# # Set parallel processing plan
# n_cores <- parallel::detectCores()
# future::plan(future::multisession, workers = n_cores)

# Load the Basque dataset from the Synth package
data(basque, package = "Synth")
data <- basque

# Rename relevant columns for clarity
colnames(data)[1:4] <- c("id", "unit", "time", "value")

# Compute additional variables
data$invest_ratio <- data$invest / data$value
data$value_raw <- data$value

# Define special predictors for the model
special_preds <- expression(list(
  list(dep.var, 1960:1969, c("mean")),
  list("invest_ratio", 1964:1969, c("mean")),
  list("popdens", 1969, c("mean")),
  list("sec.agriculture", 1961:1969, c("mean")),
  list("sec.energy", 1961:1969, c("mean")),
  list("sec.industry", 1961:1969, c("mean")),
  list("sec.construction", 1961:1969, c("mean")),
  list("sec.services.venta", 1961:1969, c("mean")),
  list("sec.services.nonventa", 1961:1969, c("mean")),
  list("school.illit", 1964:1969, c("mean")),
  list("school.prim", 1964:1969, c("mean")),
  list("school.med", 1964:1969, c("mean")),
  list("school.high", 1964:1969, c("mean")),
  list("school.post.high", 1964:1969, c("mean"))
))

# Execute the DSC analysis
system.time({
  result <- dsc(
    data = data,
    start.time = 1955,
    end.time = 1997,
    treat.time = 1970,
    dependent = "Basque Country (Pais Vasco)",
    predictors = NULL,
    parallel = TRUE,
    special.predictors = special_preds,
    time.predictors.prior = 1955:1969,
    time.optimize.ssr = 1955:1969
  )
})["elapsed"]



