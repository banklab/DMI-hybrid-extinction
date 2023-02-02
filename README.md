# DMI-hybrid-extinction
## Simulations
This repository contains the code for the simulations and the analyis relative of the hard-selection DMI purging project.

The simulation program was written in C++ using the GSL library. The source code is provided in the hybrid_speciation_hard_selection.cpp file.
One can compile the source code using the following command:
g++ -Wall -lgsl -O3 -o fourloci_hs hybrid_speciation_hard_selection.cpp

An example of a bash script used to run the program, and explore the parameter space is provided (run_fourloci_example.sh).

The outcome is written in a folder specified by the user, and the parameters used are always saved in the output file. An example of the output of the simualtions is provided here (example_output.txt).
The first line corresponds to the parameters and seed used in the simulations. The following lines corresponds to independent replciates of the simulations. There is a total of 14 numbers per line. The first 8 numbers correspond to the allele that fixed (0 or 1) and the time of fixation of the alleles at the different loci (the first two numbers corresponds to locus 1, 3 and 4 to locus 2, 5 and 6 to locus 3 and 7 and 8 to locus 4). The 9th number is the maximum population size and the 10th number the smallest population size reached during the simulation. The last 4 numbers correspond to the population size when the polymorphism is lost at eahc locus.

## Analysis
The analysis was conducted using Rstudio (v2022.07.2+576) and R(v4.2.0). The main file (analysis.Rmd) contains the code for the analysis and the output can be seen in the corresponding html file (analysis.html). R functions used to format the data are provided in their own file (functions.R).

In particular, for efficiency, the output files are read once and their result compiles in csv file that can be readily read in R. The later are provided here (whole_dataset.csv for 1,000 replicates per parameter combination, whole_dataset_extra.csv for 10,000 replicates per parameter combination and ER_dataset.csv for the data related to the evolutionary rescue scenario).

## Usage
Below, we present the steps to follow to investigate a specific parameter combinations.
- The first step is to compile the source code and make the resulting file executable.
- The second step is to edit the bash script and choose the parameters to investigate.
- The third step is to run the "add_dataframe()" function in R. This fuctions compute many relevant metrics over the independent replicates generate by the simulation program, such as the survival probability, the hybrid speciation probability, the mean time of resolution of the genetic conflict(s), the mean minimum population size. This function takes three arguments, the folder name, the dataframe to store the newly computed metrics, and the number of replicates accepted.
