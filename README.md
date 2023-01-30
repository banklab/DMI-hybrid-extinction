# DMI-hybrid-extinction
## Simulations
This repository contains the code for the simulations and the analusis relative of the hard-selection DMI purging project.

The simulation program was written in C++ using the GSL library (see the hardselection_readme.txt document for compilation details). The source code is provided in the hybrid_speciation_hard_selection.cpp file.

An example of a bash script used to run the program and explore the parameter space is provided (run_fourloci_example.sh).

The outcome is written in a folder specify by the user, and the aprameter used are always saved in the output file.

## Analysis
The analysis was conducted using Rstudio and R. The main file (analysis.Rmd) contains the code for the analysis and the output cna be seen in the corresponding html file (analysis.html). R functions used to format the data are provided in their own file (functions.R).

In particular, for efficiency, the output files are read once and their result compiles in csv file that can be readily read in R. The later are provided here (whole_dataset.csv for 1,000 replicates per parameter combination, whole_dataset_extra.csv for 10,000 replicates per parameter combination and ER_dataset.csv for the data related to the evolutionary rescue scenario).

