# DMI-hybrid-extinction
This repository contains the code for the simulations and the analysis relative to the "In search of the Goldilocks zone for hybrid speciation II: hard times for hybrid speciation?" manuscript.
## Simulations

The simulation program was written in C++ using the GSL library. The source code is provided in the hybrid_speciation_hard_selection.cpp file. On Unix system, one can compile the source code using the following command: g++ -Wall -O3 -o fourloci_hs hybrid_speciation_hard_selection.cpp -lgsl

An example of a bash script used to run the program and explore the parameter space is provided (run_fourloci_example.sh).

The outcome is written in a folder specified by the user, and the parameters used are always saved in the output file. Importantly, it is strongly recommended to always write the output in a new folder. An example of the output of the simulations is provided here (example_output.txt). The first line corresponds to the parameters and seed used in the simulations. The following lines corresponds to independent replicates of the simulations. There is a total of 14 numbers per line. The first 8 numbers correspond to the allele that fixed (0 or 1) and the time of fixation of the alleles at the different loci (the first two numbers correspond to locus 1, 3 and 4 to locus 2, 5 and 6 to locus 3 and 7 and 8 to locus 4). The 9th number is the maximum population size and the 10th number the smallest population size reached during the simulation. The last 4 numbers correspond to the population size when the polymorphism is lost at each locus.

## Analysis

The analysis was conducted using RStudio (v2022.07.2+576) and R (v4.2.0). The main file (analysis.Rmd) contains the code for the analysis and the output can be seen in the corresponding html file (analysis.html). R functions used to format the data are provided in their own file (functions.R).

For efficiency, the (simulation) output files are read once and their result compiled in a csv file that can be readily read in R. The later are provided here (whole_dataset.csv for 1,000 replicates per parameter combination, whole_dataset_extra.csv for 10,000 replicates per parameter combination and ER_dataset.csv for the data related to the evolutionary rescue scenario). The known_folder.txt,  known_folder_extra.txt and know_folder_f.txt keep track of the folders already processed. Finally, the deterministic trajectories are provided in the "deterministic_output" archive.

## Usage

Below, we present the steps to follow to investigate a specific parameter combination.

  - The first step is to compile the source code and make the resulting file executable.
- The second step is to edit the bash script and choose the parameters to investigate.
-   The third step is to run the "add_dataframe()" function in R. This function computes many relevant metrics over the independent replicates generate by the simulation program, such as the survival probability, the hybrid speciation probability, the mean time of resolution of the genetic conflict(s), the mean minimum population size. This function takes three arguments, the folder name, the dataframe to store the newly computed metrics, and the number of replicates accepted.
