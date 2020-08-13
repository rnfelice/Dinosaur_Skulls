# Dinosaur_Skulls
Code from: 

Felice RN, Watanabe A, Cuff AR, Hanson M, Bhullar B-AS, Rayfield ER, et al. (2020) Decelerated dinosaur skull evolution with the origin of birds. PLoS Biol 18(8): e3000801. <https://doi.org/10.1371/journal.pbio.3000801> 



## Contents:
### /Scripts:

Analyses_for_9_Module_data.R - calculate rate and dispaity in non-avian dinosaur 9 module dataset, generate figs 3 and S42-44

Analyses_for_11_Module_data.R - calculate rate and dispaity in non-avian dinosaur 11 module dataset, generate figs 3, 4 and S42-44

BayesTraits_Script.cmd- control file for BayesTraits Variable Rates models

BT_execution.txt - commands for running Bayestraits to generate Fig 1 and Figs S2-S36

Plotting_bayestraits_results.R - code for plotting trees with branches colored by rate (Fig 1 and S2-S36)

Plotting_global_GPA_results.R - code for calculating mean rate scalars per group (Fig 2 and S37-S40)

### /Data:
/Coordinate_data/ - contains Procrustes aligned coordinate data for 9 module and 11 module datasets, mean shape for 11 module dataset, and tables defining module identities of each landmark

/For_Figures_1_and_S2-S22/ - Phylogenetic trees and Phylogenetic PC scores for each module for birds and non-avian dinosaurs together

/For_Figures_S23-S23/ - Phylogenetic trees and Phylogenetic PC scores for each module for non-avian dinosaurs only

/Rate_Comparison/ - output files and input trees for BayesTraits RJMCMC analyses, used with "Plotting_global_GPA_results.R"

/Raw_data/ - contains full phylogenic topologies, raw data for figures 2 and s38, and a csv table with group identity (ie, bird, non avian theropod, etc)

