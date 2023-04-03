# Dynamic Synthetic Controls: Accounting for Varying Speeds in Comparative Case Studies

## 1. Main Results

### 1.1 Speed Problem
* R/figure_speed_problem.R --> plot Figure 1.
* R/figure_speed_problem_v2.R --> plot Figure A1.

### 1.2 Two Fold Dynamic Time Warping
* R/figure_TFDTW_method.R --> plot Figure 2.

### 1.3 Monte Carlo Study
* R/figure_sim_sample.R --> plot Figure 3.
* R/run_sim_beta_0.R --> Set beta=0 and run simulation study.
* R/run_sim_beta_05.R --> Set beta=0.5 and run simulation study.
* R/run_sim_beta_1.R --> Set beta=1, run simulation study, and plot Figure 4.

### 1.4 Re-evaluating Empirical Findings
* R/run_basque.R --> run placebo test on Basque data.
* R/run_tobacco.R --> run placebo test on tobacco data.
* R/run_germany.R --> run placebo test on Germany data.
* R/figure_placebo.R --> plot Figure 5.

## 2. Sample Run
* R/run_sample_fast.R --> compare SC and DSC on Basque Country data (20 seconds)
* R/run_sample_slow.R --> compare SC and DSC on Basque Country data (full process)
