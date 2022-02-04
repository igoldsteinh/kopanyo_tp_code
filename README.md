# kopanyo_tp_code

This repository has all code and data required to generate the results reported in [Investigating methods for using genetic data to estimate the association between host factors and being an infection source with application to probing the association between HIV infection and tuberculosis transmission](https://www.medrxiv.org/content/10.1101/2021.12.12.21267687v1). 

## Navigation

### Simulation Study
Files used to generate simulated TB-like outbreaks are located in the [simulation folder](https://github.com/igoldsteinh/kopanyo_tp_code/tree/main/R/simulation_accuracy_and_performance/simulation_files). These files were written to be run on a computing cluster using the [Slurm scheduler](https://slurm.schedmd.com/slurm.html). Pre-processed results assessing odds ratio estimation pipeline performance used in the manuscript are stored in the [simulation results folder](https://github.com/igoldsteinh/kopanyo_tp_code/tree/main/R/simulation_accuracy_and_performance/simulation_results), and can be visualized using the available [R Code](https://github.com/igoldsteinh/kopanyo_tp_code/blob/main/R/simulation_accuracy_and_performance/visualize_simulation_results.R). 

Code used to calculate the true odds ratios of interest using Monte Carlo simulation are available in the [monte carlo folder](https://github.com/igoldsteinh/kopanyo_tp_code/tree/main/R/monte_carlo_OR_estimation). 

Code used to summarise the characteristics of simulated TB-like outbreaks are located in the [simulation summary folder](https://github.com/igoldsteinh/kopanyo_tp_code/tree/main/R/simulation_summary_table)

### Real Data Analysis
SNP Data and associated metadata used in the real data results portion of the manuscript are in the [data folder](https://github.com/igoldsteinh/kopanyo_tp_code/tree/main/data). BEAST2 .xml files and maximum clade credibility trees are located in the [BEAST2 folder](https://github.com/igoldsteinh/kopanyo_tp_code/tree/main/BEAST2). R Code for the creation of clusters by SNP cutoff, as well as analysis and visualization of real data are located in [real_data_analysis folder](https://github.com/igoldsteinh/kopanyo_tp_code/tree/main/R/real_data_analysis). 