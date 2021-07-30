Readme.txt

—————————————————————————————————————————————
    Compilation and Running the Source Code                                        
—————————————————————————————————————————————

—Install the GNU Scientific Library. On macOS you could do this via homebrew, by running 
>>brew install gel

—To generate an executable, compile and run the C++ code, as given below,
>> g++ -Wall -lgsl -O3 -o fourloci_hs hybrid_speciation_hard_selection.cpp

—This generates an executable version of the code - fourloci_hs

—To make the bash script executable, run, 
>> chmod +x run_fourloci.sh

—To run a simulation 
>> ./run_fourloci.sh



—————————————————————————————————————————————
       Changing simulation parameters                    						
—————————————————————————————————————————————

— The simulation parameters can be easily edited in the bash script run_fourloci.sh. 

— The seed, selection coefficient, epistasis (codominant or recessive), recombination, number of generations, number of individuals (at t=0), maximum population size, haplotype contribution and the number of iterations per parameter combination can be controlled. 

— We are considering a diallelic four loci model there are a total of 16 possible haplotypes that can be obtained. We describe these 16 haplotypes using a binary conversion system, where ‘0’ corresponds to an ancestral allele and ‘1’ for a derived allele

ex. The haplotype a1B1a2B2 is recast as 0101, which in binary is equal to 6. In the bash script this corresponds to f6. 

— For each linkage architecture, ensure the correct parental haplotypes are chosen and for equal contribution of parental haplotypes, divide n_gen equally between them.

ex. For the Adjacent “ABAB” linkage architecture, there are two parental haplotypes - a1B1a2B2 and A1b1A2b2, which correspond to f6 and f11 using the method described above.



—————————————————————————————————————————————
         .txt file interpretation          						
—————————————————————————————————————————————

— The name of the output folder (for the simulation) and individual .txt files within each folder can be customised to the users preferences by making changes in the bash script. The default setting will append the recombination rate (varied between 0.0001 to 0.5) to the file name of the .txt files 

— The bash script also generates a file called parameter.txt within the output folder that contains a line with the parameters used for each simulation (for reference)

— The output in the .txt file is similar to what is given below:

Parameters : (…)
73 0 107 0 234 0 83 1001433
0 53 0 111 0 65 1 11515 1005137
0 146 0 95 0 559 0 122 0
0 53 0 124 0 125 0 63 1001137
0 95 0 95 0 36 0 280 1000752
0 50 0 388 0 36 1 5408 1003960
0 95 0 51 0 49 0 154 0
0 121 0 42 0 75 1 4095 1002783
0 146 0 44 0 53 0 72 1003358
…

— Each .txt file contains a line with the simulation parameters, and each line of the file corresponds to a single simulation. 

—Each column represents:
Allele that fixes at locus 1 (here A1), Time of fixation at locus 1, Allele that fixes at locus 2 (here B1) Time of fixation at locus 2,  Allele that fixes at locus 3 (here B2), Time of fixation at locus 3, Allele that fixes at locus 4 (here A2), Time of fixation at locus 4, n_ind in population with the fixed haplotype

— ‘0’ indicates to the fixation of an ancestral allele and ‘1’ indicates the fixation of a derived allele. 




