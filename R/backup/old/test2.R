folder = "./data/pred/tobacco/"
file.list = as.list(list.files(folder))

pre.start = 22
pre.end = 31
post.start = 32
post.end = 41

value_raw_list = data.list %>% 
  future_map(
    ~{
      data = .[["data"]]
      data %>% arrange(id, time) %>% select(id, unit, time, value_raw)
    }
  )

# df.mse
df.mse = file.list %>% 
  future_map(
    ~{
      file.name = .
      result.list = readRDS(paste0(folder, file.name))
      data = result.list[[1]]$args.TFDTW.synth$data
      data = data %>% arrange(id, time) %>% select(id, unit, time, value_raw)
      data_in_list = sapply(value_raw_list,
                            function(x) identical(x, data))
      data.id = which(data_in_list)
      
      mse = future_map2(
        result.list,
        as.list(names(result.list)),
        ~{
          result = .x
          result.synth = result[["results.TFDTW.synth"]]
          grid.id = .y
          filter.width = result$filter.width
          k = result$args.TFDTW.synth$args.TFDTW$k
          step.pattern = result$args.TFDTW.synth$args.TFDTW$step.pattern1
          step.pattern_in_list = sapply(step.pattern.range,
                               function(x) identical(x, step.pattern))
          step.pattern.id = which(step.pattern_in_list)

          mse = result.synth %>% 
            map(
              ~{
                task = .
                unit = task$dependent
                ###
                scales = df.rescale %>% filter(unit == task$dependent)
                # value.min = scales$value.min
                multiplier = scales$multiplier
                gap.raw = task$gap.raw/multiplier
                gap.TFDTW = task$gap.TFDTW/multiplier
                ###
                data.frame(unit = unit,
                           mse.preT.raw = mean(gap.raw[pre.start:pre.end]^2, na.rm = T),
                           mse.preT.TFDTW = mean(gap.TFDTW[pre.start:pre.end]^2, na.rm = T),
                           mse.postT.raw = mean(gap.raw[post.start:post.end]^2, na.rm = T),
                           mse.postT.TFDTW = mean(gap.TFDTW[post.start:post.end]^2, na.rm = T))
              }
            ) %>% do.call("rbind", .)
          mse %>% mutate(grid.id = grid.id,
                         filter.width = filter.width,
                         k = k,
                         step.pattern.id)
        }
      ) %>% do.call("rbind", .)
      if (length(data.id) > 0) {
        mse$data.id = data.id
      }else{
        mse$data.id = NA
      }

      mse %>% 
        group_by(unit) %>% 
        top_n(-1, mse.preT.TFDTW) %>% 
        top_n(-1, grid.id)
    }
  ) %>% do.call("rbind", .)

df.mse = df.mse %>% mutate(step.pattern = names(step.pattern.range)[step.pattern.id])
grid.opt = left_join(df.mse[,c(1,6,7,8,10,11)], distinct(data[, c(1,2)]), by = c("unit"))

grid.opt = grid.opt %>% filter(data.id != "1") %>% 
  select(id, unit, data.id, grid.id, filter.width, k, step.pattern) %>% 
  ungroup

a = grid.opt %>% select(filter.width, k, step.pattern)
a = data.frame(a)
b = search.grid[grid.opt$grid.id,]
b$step.pattern = as.character(b$step.pattern)
rownames(b) = NULL
colSums(a == b)
print(nrow(a))

a = grid.opt %>% select(unit, grid.id, data.id, id)
b = grid.opt2

saveRDS(grid.opt, "./data/Figure_A3_3_gridOpt.Rds")


# Set the working directory to the folder containing your R scripts
setwd("/Users/jiancao/Documents/GitHub/Dynamic-Synthetic-Controls/R/rep/script")  # Change this to your folder path

# List all R scripts in the folder
r_scripts <- list.files(pattern = "\\.R$")

# Function to extract package names from a script
extract_packages <- function(script) {
  content <- readLines(script, warn = FALSE)
  # Extract packages used in library() and require()
  library_packages <- unlist(regmatches(content, gregexpr("(?<=library\\()\\s*\\w+", content, perl = TRUE)))
  require_packages <- unlist(regmatches(content, gregexpr("(?<=require\\()\\s*\\w+", content, perl = TRUE)))
  # Extract packages used with :: or :::
  colon_packages <- unlist(regmatches(content, gregexpr("\\b\\w+(?=:::?)", content, perl = TRUE)))
  
  unique(c(library_packages, require_packages, colon_packages))
}

# Initialize an empty vector to store all package names
all_packages <- c()

# Loop through each script and extract packages
for (script in r_scripts) {
  script_packages <- extract_packages(script)
  all_packages <- c(all_packages, script_packages)
}

# Get unique package names
unique_packages <- unique(all_packages)

# Display the unique packages
print(unique_packages)








