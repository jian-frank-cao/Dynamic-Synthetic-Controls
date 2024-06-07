set.seed(20220407)
i = 1
args.TFDTW.synth.target.only[["data"]] = data.list[[i]]
res = SimDesign::quiet(
  grid.search.opt(filter.width.range = filter.width.range,
                  k.range = k.range,
                  step.pattern.range = step.pattern.range,
                  args.TFDTW.synth.target.only = args.TFDTW.synth.target.only,
                  grid.search.parallel = grid.search.parallel)
)

ind = 50
res = SimDesign::quiet(
  grid.search.opt(filter.width.range = search.grid[ind,]$filter.width,
                  k.range = search.grid[ind,]$k,
                  step.pattern.range = 
                    step.pattern.range[search.grid[ind,]$step.pattern],
                  args.TFDTW.synth.target.only = args.TFDTW.synth.target.only,
                  grid.search.parallel = grid.search.parallel)
)


filenames = list.files(folder)

a = strsplit(filenames, "_")
b = NULL
for (item in a) {
  target = item[max(length(item))]
  b = c(b, as.numeric(strsplit(target, ".Rds")[[1]]))
}

grid.opt = left_join(df.mse[,c(1,6,7)], distinct(data[, c(1,2)]), by = c("unit"))

grid.opt = df.mse[, c(6,1)] %>% arrange(data.id) %>% `rownames<-`(NULL)

grid.opt = grid.opt %>% filter(data.id != "1")
saveRDS(grid.opt, "./data/Figure_A3_3_gridOpt.Rds")

i = 1
set.seed(20220407)
args.TFDTW.synth.target.only[["data"]] = data.list[[grid.opt$data.id[i]]]
aa = grid.search.opt(filter.width.range = 
                  search.grid[grid.opt$grid.id[i],]$filter.width,
                k.range = search.grid[grid.opt$grid.id[i],]$k,
                step.pattern.range =
                  step.pattern.range[search.grid[grid.opt$grid.id[i],]$step.pattern],
                args.TFDTW.synth.target.only = args.TFDTW.synth.target.only,
                grid.search.parallel = grid.search.parallel)


set.seed(20220407)
args.TFDTW.synth.all.units = list(target = "A",
                                  # data = data, 
                                  args.TFDTW.synth = args.TFDTW.synth,
                                  ## 2nd
                                  all.units.parallel = TRUE)
args.TFDTW.synth.all.units[["data"]] = data.list[[grid.opt$data.id[i]]]
bb = SimDesign::quiet(
  grid.search(filter.width.range = 
                search.grid[grid.opt$grid.id[i],]$filter.width,
              k.range = search.grid[grid.opt$grid.id[i],]$k,
              step.pattern.range =
                step.pattern.range[search.grid[grid.opt$grid.id[i],]$step.pattern],
              args.TFDTW.synth.all.units = args.TFDTW.synth.all.units,
              grid.search.parallel = grid.search.parallel)
)

df.mse.now2 = cbind(grid.opt[,2], df.mse, grid.opt[,1])
colnames(df.mse.now2) = colnames(df.mse.new)

grid.opt = df.mse %>% select(data.id, grid.id) %>% arrange(data.id)
saveRDS(grid.opt, "./data/Figure_4_gridOpt.Rds")


a = data.list.new[[99]][["data"]] == data.list.old[[99]][["data"]]
colSums(a)


a = data.frame(table(df.mse.old$grid.id))
b = data.frame(table(df.mse.old2$grid.id))


grid.opt = df.mse %>% select(data.id, grid.id) %>% arrange(data.id)
saveRDS(grid.opt, "./data/Figure_4_gridOpt.Rds")

df.mse$data.id = as.character(df.mse$data.id)
df2 = left_join(df.mse, df.mse.old, by = c("unit", "grid.id", "data.id"))
df = df %>% 
  mutate(ratio.preT.raw = mse.preT.raw.x/mse.preT.raw.y,
         ratio.preT.TFDTW = mse.preT.TFDTW.x/mse.preT.TFDTW.y,
         ratio.postT.raw = mse.postT.raw.x/mse.postT.raw.y,
         ratio.postT.TFDTW = mse.postT.TFDTW.x/mse.postT.TFDTW.y)

a = data.list.old[[107]]$data == data.list.new[[107]]$data
colSums(a)




set.seed(20220407)
data.id = 107
grid.id = 219
unit = "Italy"
id = 8
args.TFDTW.synth.target.only = list(target = unit, id = id,
                                    data = data.list[[data.id]][["data"]],
                                    args.TFDTW.synth = args.TFDTW.synth)
result_target2 = SimDesign::quiet(
  grid.search.opt(filter.width.range = search.grid[grid.id,]$filter.width,
                  k.range = search.grid[grid.id,]$k,
                  step.pattern.range =
                    step.pattern.range[search.grid[grid.id,]$step.pattern],
                  args.TFDTW.synth.target.only = args.TFDTW.synth.target.only,
                  grid.search.parallel = grid.search.parallel)
)

pre.start = 22
pre.end = 31
post.start = 32
post.end = 41

task = result_target2[[1]][["results.TFDTW.synth"]]
unit = task$dependent
scales = df.rescale %>% filter(unit == task$dependent)
multiplier = scales$multiplier
gap.raw = task$gap.raw/multiplier
gap.TFDTW = task$gap.TFDTW/multiplier
a = data.frame(unit = unit,
           mse.preT.raw = mean(gap.raw[pre.start:pre.end]^2, na.rm = T),
           mse.preT.TFDTW = mean(gap.TFDTW[pre.start:pre.end]^2, na.rm = T),
           mse.postT.raw = mean(gap.raw[post.start:post.end]^2, na.rm = T),
           mse.postT.TFDTW = mean(gap.TFDTW[post.start:post.end]^2, na.rm = T))








is_in_list <- sapply(step.pattern.range, function(x) identical(x, target))
step.pattern.range[which(is_in_list)]
