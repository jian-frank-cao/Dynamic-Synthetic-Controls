---
contributors:
  - Lars Vilhuber
  - Miklós Koren
  - Joan Llull
  - Marie Connolly
  - Peter Morrow
---

## Overview

The code in this replication package constructs the analysis file from the three data sources (Abadie and Gardeazabal, 2003; Abadie, Diamond, and Hainmueller, 2010; Abadie, Diamond, and Hainmueller, 2015) using R. Two main files run all of the code to generate the data for the 5 figures in the main paper and 3 figures in the Appendix. The replicator should expect the code to run for about 10 hours.

## Data Availability and Provenance Statements

- [ ] This paper does not involve analysis of external data (i.e., no data are used or the only data are generated by the authors via simulation in their code).

### Statement about Rights

- [x] I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript. 
- [x] I certify that the author(s) of the manuscript have documented permission to redistribute/publish the data contained within this replication package.


### Summary of Availability

- [x] All data **are** publicly available.
- [ ] Some data **cannot be made** publicly available.
- [ ] **No data can be made** publicly available.


- [ ] Confidential data used in this paper and not provided as part of the public replication package will be preserved for ___ years after publication, in accordance with journal policies. 

### Details on each Data Source

| Data.Name  | Data.Files | Location | Provided | Citation |
| -- | -- | -- | -- | -- | 
| “Panel Data from Spanish Regions 1955-1997” | "basque" in R package "Synth" | None | FALSE | Abadie and Gardeazabal, 2003 |
| “Cigarette Consumption in United States 1970-2000” | smoking.rda; prop99.csv | ./data/ | TRUE | Abadie, Diamond, and Hainmueller, 2010 |
| “OECD Country Panel Data 1960-2003” | repgermany.Rds | ./data/ | TRUE | Abadie, Diamond, and Hainmueller, 2015 |

The data set "basque" is publicly available in R package "Synth". The R package "Synth" can be obtained from [CRAN](https://cran.r-project.org/web/packages/Synth/index.html).

The data files "smoking.rda" and "prop99.csv" are publicly available at [Github](https://github.com/tom-beer/Synthetic-Control-Examples?tab=readme-ov-file). The files are included in data archive.

The data file "repgermany.tab" is available for replication at [Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/24714). The data file was downloaded and an R version "repgermany.Rds" is included in data archive.

## Computational requirements

### Software Requirements

- [x] The replication package contains one or more programs to install all dependencies and set up the necessary directory structure. [HIGHLY RECOMMENDED]

- R 4.2.0
  - `dtw` (1.23-1)
  - `emojifont` (0.5.5)
  - `forecast` (8.21.1)
  - `furrr` (0.3.1)
  - `ggpubr` (0.6.0)
  - `gridExtra` (2.3)
  - `purrr` (1.0.2)
  - `reshape2` (1.4.4)
  - `signal` (0.7-7)
  - `SimDesign` (2.15.1)
  - `Synth` (1.1-8)
  - `tidyverse` (2.0.0)
  - `TTR` (0.24.3)
  - `zoo` (1.8-12)
  - the file "`00_setup.R`" will install all dependencies (latest version), and should be run once prior to running other programs.

### Controlled Randomness

- [x] Random seed is set at line 2 of program "`00_setup.R`"
- [ ] No Pseudo random generator is used in the analysis described here.

### Memory, Runtime, Storage Requirements

#### Summary

Approximate time needed to reproduce the analyses on a standard (CURRENT YEAR) desktop machine:

- [ ] <10 minutes
- [ ] 10-60 minutes
- [ ] 1-2 hours
- [ ] 2-8 hours
- [x] 8-24 hours
- [ ] 1-3 days
- [ ] 3-14 days
- [ ] > 14 days

Approximate storage space needed:

- [ ] < 25 MBytes
- [ ] 25 MB - 250 MB
- [x] 250 MB - 2 GB
- [ ] 2 GB - 25 GB
- [ ] 25 GB - 250 GB
- [ ] > 250 GB

- [ ] Not feasible to run on a desktop machine, as described below.

#### Details

The code was last run on a **8-core Apple M1 Pro laptop with MacOS version 14.5 with 300GB of free space**.  Computation took **11 hours**.

The code was run on a **12-core Intel desktop with Windows version 11 with 24GB of RAM, 200 GB free space**.  Computation took **20 hours**.

## Description of programs/code

- Programs in `00_setup.R` install necessary packages for the other program files in the package.
- Programs in `01_main_figures.R` source and run program files `Figure_1.R`, `Figure_2.R`, `Figure_3.R`, `Figure_4.R`, `Figure_5_1.R`, `Figure_5_2.R`, `Figure_5_3.R`, `Figure_5_all.R`, output `Figure_1.pdf`, `Figure_2.pdf`, `Figure_3.pdf`, `Figure_4.pdf`, `Figure_5.pdf`, and save them in folder `./results/`.
- Programs in `02_appendix_figures.R` source and run program files `Figure_A1.R`, `Figure_A2.R`, `Figure_A3_1.R`, `Figure_A3_2.R`, `Figure_A3_3.R`, `Figure_A3_all.R`, output `Figure_A1.pdf`, `Figure_A2.pdf`, `Figure_A3.pdf`, and save them in folder `./results/`.
- Programs in `Figure_1.R` create an illustration figure `Figure_1.pdf` for the speed problem, and save it in folder `./results/`.
- Programs in `Figure_2.R` create an illustration figure `Figure_2.pdf` for the Dynamic Synthetic Control method, and save it in folder `./results/`.
- Programs in `Figure_3.R` create an illustration figure `Figure_3.pdf` for the Dynamic Time Warping path, and save it in folder `./results/`.
- Programs in `Figure_4.R` simulate 100 data sets, apply the standard Synthetic Control method and the Dynamic Synthetic Control method, create a figure `Figure_4.pdf` for a comparison of errors, and save it in folder `./results/`.
- Programs in `Figure_5_1.R` replicate the application in Abadie and Gardeazabal (2003) using the standard Synthetic Control and the Dynamic Synthetic Control methods, create a sub-figure for Figure 5 for a comparison of errors in placebo tests, and save it into data file `Figure_5_1.Rds` in folder `./data/`.
- Programs in `Figure_5_2.R` replicate the application in Abadie, Diamond, and Hainmueller (2010) using the standard Synthetic Control and the Dynamic Synthetic Control methods, create a sub-figure for Figure 5 for a comparison of errors in placebo tests, and save it into data file `Figure_5_2.Rds` in folder `./data/`.
- Programs in `Figure_5_3.R` replicate the application in Abadie, Diamond, and Hainmueller (2015) using the standard Synthetic Control and the Dynamic Synthetic Control methods, create a sub-figure for Figure 5 for a comparison of errors in placebo tests, and save it into data file `Figure_5_3.Rds` in folder `./data/`.
- Programs in `Figure_5_all.R` load sub-figures from data files `Figure_5_1.Rds`, `Figure_5_2.Rds`, and `Figure_5_3.Rds`, create a figure `Figure_5.pdf`, and save it in folder `./results/`.
- Programs in `Figure_A1.R` create an illustration figure `Figure_A1.pdf` for the speed problem with lags and polynomials, and save it in folder `./results/`.
- Programs in `Figure_A2.R` create an illustration figure `Figure_A2.pdf` for the Monte-Carlo analysis, and save it in folder `./results/`.
- Programs in `Figure_A3_1.R` replicate the application in Abadie and Gardeazabal (2003) using the standard Synthetic Control and the Dynamic Synthetic Control methods with warped covariates, create a sub-figure for Figure A3 for a comparison of errors in placebo tests, and save it into data file `Figure_A3_1.Rds` in folder `./data/`.
- Programs in `Figure_A3_2.R` replicate the application in Abadie, Diamond, and Hainmueller (2010) using the standard Synthetic Control and the Dynamic Synthetic Control methods with warped covariates, create a sub-figure for Figure A3 for a comparison of errors in placebo tests, and save it into data file `Figure_A3_2.Rds` in folder `./data/`.
- Programs in `Figure_A3_3.R` replicate the application in Abadie, Diamond, and Hainmueller (2015) using the standard Synthetic Control and the Dynamic Synthetic Control methods with warped covariates, create a sub-figure for Figure A3 for a comparison of errors in placebo tests, and save it into data file `Figure_A3_3.Rds` in folder `./data/`.
- Programs in `Figure_A3_all.R` load sub-figures from data files `Figure_A3_1.Rds`, `Figure_A3_2.Rds`, and `Figure_A3_3.Rds`, create a figure `Figure_A3.pdf`, and save it in folder `./results/`.
- Programs in `grid.search.R` conduct grid search to find out the optimal parameters.
- Programs in `implement.R` implement the Dynamic Synthetic Control method.
- Programs in `misc.R` conduct data pre-processing, data transformation, and plot figures.
- Programs in `simulate.R` simulate artificial data sets for Monte-Carlo analysis.
- Programs in `synth.R` implement the standard Synthetic Control method.
- programs in `TFDTW.R` implement the two-fold Dynamic Time Warping method, which is the first step of the Dynamic Synthetic Control method.

## Instructions to Replicators

- Edit and run `./code/00_setup.R` to adjust the default path and install necessary packages.
- Run `./code/01_main_figures.R` to produce main paper figures 1-5. The figures are saved in folder `./results/`.
- Run `./code/02_appendix_figures.R` to produce appendix figures A1-A3. The figures are saved in folder `./results/`.


## List of tables and programs

The provided code reproduces:

- [ ] All numbers provided in text in the paper
- [x] All tables and figures in the paper
- [ ] Selected tables and figures in the paper, as explained and justified below.


| Figure/Table #    | Program                  | Line Number | Output file                              | Note                            |
|-------------------|--------------------------|-------------|------------------------------------------|---------------------------------|
| Figure 1          | ./code/Figure_1.R             |             | ./results/Figure_1.pdf                     |                                 |
| Figure 2          | ./code/Figure_2.R             |             | ./results/Figure_2.pdf                     |                                 |
| Figure 3          | ./code/Figure_3.R             |             | ./results/Figure_3.pdf                     |                                 |
| Figure 4          | ./code/Figure_4.R             |             | ./results/Figure_4.pdf                     |                                 |
| Figure 5          | ./code/Figure_5_all.R         |             | ./results/Figure_5.pdf                     |                                 |
| Figure A1         | ./code/Figure_A1.R            |             | ./results/Figure_A1.pdf                    |                                 |
| Figure A2         | ./code/Figure_A1.R            |             | ./results/Figure_A2.pdf                    |                                 |
| Figure A3         | ./code/Figure_A3_all.R        |             | ./results/Figure_A3.pdf                    |                                 |


## References

Abadie, Alberto, and Javier Gardeazabal. 2003. The economic costs of conflict: a case study of the basque country. American Economic Review 93 (1): 113–132.

Abadie, Alberto, Alexis Diamond, and Jens Hainmueller. 2010. Synthetic control methods for comparative case studies: estimating the effect of California’s tobacco control program. Journal of the American Statistical Association 105 (490): 493–505.

Abadie, Alberto, Alexis Diamond, and Jens Hainmueller. 2015. Comparative politics and the synthetic control method. American Journal of Political Science 59 (2): 495–510.

---

## Acknowledgements

Some content on this page was copied from [Hindawi](https://www.hindawi.com/research.data/#statement.templates). Other content was adapted  from [Fort (2016)](https://doi.org/10.1093/restud/rdw057), Supplementary data, with the author's permission.