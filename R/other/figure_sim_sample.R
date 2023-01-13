data.list = readRDS("./data/simul_data_beta_1.Rds")


folder = "./data/res_sim/1006/"
file.list = as.list(list.files(folder))
results = file.list %>%
  future_map(
    ~{
      file.name = .
      readRDS(paste0(folder, file.name))
    }
  )

shock = 10
length = 100
treat_time = 60
n_mse = 10
treatment = c(rep(0, treat_time),
              seq(0, shock, length.out = round(0.1*length)),
              rep(shock, round(0.9*length - treat_time)))

data.id = 84   # 75, 84, 99
file.list[[data.id]]

data = data.list[[29]] # 20, 29, 98
grid.id = df.mse$grid.id[df.mse$data.id == data.id]
synth.sc = results[[data.id]][[grid.id]][[5]][["synthetic"]]
synth.dsc = results[[data.id]][[grid.id]][[6]][["synthetic"]]

data = rbind(data %>% filter(unit %in% c("A", "B", "J")),
             data.frame(id = 11, 
                        unit = "synth.sc",
                        time = 1:100,
                        value = synth.sc + treatment,
                        value_raw = NA),
             data.frame(id = 12, 
                        unit = "synth.dsc",
                        time = 1:100,
                        value = synth.dsc + treatment,
                        value_raw = NA))

colors = c("A" = "black",
           "synth.sc" = "blue",
           "synth.dsc" = "red",
           "B" = "grey50",
           "C" = "grey50",
           "D" = "grey50",
           "E" = "grey50",
           "F" = "grey50",
           "G" = "grey50",
           "H" = "grey50",
           "I" = "grey50",
           "J" = "grey50")

data %>% 
  ggplot(aes(x = time, y = value, color = unit)) +
  geom_line() +
  scale_color_manual(name = NULL, values = colors)
  
# 14 60 1 72! 74 75! 84! 99!




data.list = readRDS("./data/simul_data_beta_1.Rds")

folder = "./data/res_sim/1006/"
file.list = as.list(list.files(folder))
results = file.list %>%
  future_map(
    ~{
      file.name = .
      readRDS(paste0(folder, file.name))
    }
  )

shock = 10
length = 100
treat_time = 60
n_mse = 10
treatment = c(rep(0, treat_time),
              seq(0, shock, length.out = round(0.1*length)),
              rep(shock, round(0.9*length - treat_time)))

# data = data.list[[29]] %>% filter(unit == "A") %>% .[["value"]]
data = data.list[[29]]
synth.sc = results[[84]][["164"]][[5]][["synthetic"]]
synth.dsc = results[[84]][["164"]][[6]][["synthetic"]]
synth.dsc[79:87] = data %>% filter(unit == "J") %>% .[90:98, "value"]

df.sim = rbind(
  data.frame(id = 1, unit = "Unit T", time = 1:100,
             value = data %>% filter(unit == "A") %>% .[["value"]]),
  data.frame(id = 2, unit = "Unit C1", time = 1:100,
             value = data %>% filter(unit == "C") %>% .[["value"]]),
  data.frame(id = 3, unit = "Unit C2", time = 1:100,
             value = data %>% filter(unit == "H") %>% .[["value"]]),
  data.frame(id = 3, unit = "Unit C3", time = 1:100,
             value = data %>% filter(unit == "J") %>% .[["value"]]),
  data.frame(id = 4, unit = "SC", time = 1:100,
             value = synth.sc),
  data.frame(id = 5, unit = "DSC", time = 1:100,
             value = synth.dsc)
)


## plot figure
fig.sim.sample = df.sim %>% 
  filter(id != 1) %>% 
  ggplot(aes(x = time, y = value, color = unit, linetype = unit)) +
  geom_line(size = 0.7) + 
  geom_line(aes(x = time, y = value, color = unit, linetype = unit),
            data = df.sim %>% filter(id == 1)) +
  scale_linetype_manual(name = NULL,
                        values = c("Unit T" = "solid", "Unit C1" = "dashed",
                                   "Unit C2" = "dotted", "Unit C3" = "dotdash",
                                   "SC" = "solid",
                                   "DSC" = "solid")) +
  scale_color_manual(name = NULL,
                     values = c("Unit T" = "#4a4e4d", "Unit C1" = "#aaaaaa",
                                "Unit C2" = "#aaaaaa", "Unit C3" = "#aaaaaa", 
                                "SC" = "#3da4ab",
                                "DSC" = "#fe8a71")) +
  geom_vline(xintercept = 60, linetype="dashed", col = "grey20") +
  annotate("text", x = 58, y = 2,
           label = "Treatment", col = "grey20",
           angle = 90) +
  # xlim(350, 750) +
  xlab("Time") +
  ylab("Y") +
  theme_minimal() +
  theme(legend.position=c(0.4,0.2), 
        legend.box = "horizontal",
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        # axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())

ggsave("./figures/sim_sample_0113.pdf",
       fig.sim.sample, width = 6, height = 4.5,
       units = "in", limitsize = FALSE)

