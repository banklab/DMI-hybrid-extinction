#!/bin/bash

# This an example of the bash script file to run the simulation program. It allows to choose all parameters of the model.

# all parameters value are defines here

OFFSPRING=1.01 # default expected growth of the ancestral genotype
S1=-0.001 # selection coefficient locus 1
S2=-0.001 # selection coefficient locus 2
S3=-0.001 # selection coefficient locus 3
S4=-0.001 # selection coefficient locus 4

E12=-0.99 # epistasis between loci 1 and 2
E13=0    # epistasis between loci 1 and 3
E14=0    # epistasis between loci 1 and 4
E23=0    # epistasis between loci 2 and 3
E24=0    # epistasis between loci 2 and 4
E34=-0.99 # epistasis between loci 3 and 4

R12=0.5 # recombination between loci 1 and 2
R23=0.5 # recombination between loci 2 and 3
R34=0.5 # recombination between loci 3 and 4

NGEN=10000000 # maximum number of generation to stop the simulation. 100 times the number of haplotype is usually safe. Ideally, this value should never be reached. All fixation should happen before.
NIND=500 # Number of haplotype in the populations. The population size is half this value.
MAX_POP_SIZE=1000000 # Maximum population size

H12=1 # 1 for codominance between loci 1 and 2 and 0 for recessivity
H13=1 # 1 for codominance between loci 1 and 3 and 0 for recessivity
H14=1 # 1 for codominance between loci 1 and 4 and 0 for recessivity
H23=1 # 1 for codominance between loci 2 and 3 and 0 for recessivity
H24=1 # 1 for codominance between loci 2 and 4 and 0 for recessivity
H34=1 # 1 for codominance between loci 3 and 4 and 0 for recessivity

f1=0 #hap0000
f2=0 #hap0001
f3=0 #hap0010
f4=0 #hap0011
f5=0 #hap0100
f6=250 #hap0101
f7=0 #hap0110
f8=0 #hap0111
f9=0 #hap1000
f10=0 #hap1001
f11=250 #hap1010
f12=0 #hap1011
f13=0 #hap1100
f14=0 #hap1101
f15=0 #hap1110
f16=0 #hap1111

nb_ite=1000 # nb_ite per parameter combination



name="output_sym_cod_rfunction_pi_0d5_lethal_small" # name of the output folder

# check that the folder does not exisit to avoid any overwriting issue
if test -d "$name"/ 
then
  echo "Please rename or remove the directory $name/ first, then rerun"
  exit 1
fi

mkdir $name # create the folder
cd $name # move into the folder

# Here we choose to make recombination varies and have all loci equidistant. Generate a set of simulations for all values of r between 0.1 to 0.5 with an increment of 0.01
for (( i=10; i<=50; i+=1 ))
  do
  # remove any preexisting parameter file
  test -f parameter.txt &&  rm -f parameter.txt
  # write the parameter file
  echo "$OFFSPRING $S1 $S2 $S3 $S4 $E12 $E13 $E14 $E23 $E24 $E34 .$i .$i .$i $H12 $H13 $H14 $H23 $H24 $H34 $NGEN $NIND $MAX_POP_SIZE $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 $f16" >>parameter.txt
  # print in the console the current stage for the simulation
  echo ".$i"
  # call for the program; output file has the same name than the folder plus an extra part corresponding to the recombination rate for easy identification. Note that since the parameters value are written in the file itself a mistake here, while not desirable is not critical.
  ./../fourloci_hs parameter.txt "${name}_rdot${i}.txt" $nb_ite -1
  done

# we choose to make recombination varies and have all loci equidistant. Generate a set of simulations for all values of r between 0.01 to 0.09 with an increment of 0.01

for (( i=1; i<=9; i+=1 ))
  do
   test -f parameter.txt && rm -f parameter.txt
   echo "$OFFSPRING $S1 $S2 $S3 $S4 $E12 $E13 $E14 $E23 $E24 $E34 .0$i .0$i .0$i $H12 $H13 $H14 $H23 $H24 $H34 $NGEN $NIND $MAX_POP_SIZE $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 $f16" >>parameter.txt
   echo ".$i"
   ./../fourloci_hs parameter.txt "${name}_rdot0${i}.txt" $nb_ite -1
  done

# we choose to make recombination varies and have all loci equidistant. Generate a set of simulations for all values of r between 0.001 to 0.009 with an increment of 0.001
for (( i=1; i<=9; i+=1 ))
  do
   test -f parameter.txt && rm -f parameter.txt
   echo "$OFFSPRING $S1 $S2 $S3 $S4 $E12 $E13 $E14 $E23 $E24 $E34 .00$i .00$i .00$i $H12 $H13 $H14 $H23 $H24 $H34 $NGEN $NIND $MAX_POP_SIZE $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 $f16" >>parameter.txt
   echo ".$i"
   ./../fourloci_hs parameter.txt "${name}_rdot00${i}.txt" $nb_ite -1
  done


# we choose to make recombination varies and have all loci equidistant. Generate a set of simulations for all values of r between 0.0001 to 0.0009 with an increment of 0.0001 
for (( i=1; i<=9; i+=1 ))
  do
   test -f parameter.txt && rm -f parameter.txt
   echo "$OFFSPRING $S1 $S2 $S3 $S4 $E12 $E13 $E14 $E23 $E24 $E34 .000$i .000$i .000$i $H12 $H13 $H14 $H23 $H24 $H34 $NGEN $NIND $MAX_POP_SIZE $f1 $f2 $f3 $f4 $f5 $f6 $f7 $f8 $f9 $f10 $f11 $f12 $f13 $f14 $f15 $f16" >>parameter.txt
   echo ".$i"
   ./../fourloci_hs parameter.txt "${name}_rdot000${i}.txt" $nb_ite -1
  done

# The whole set up was chosen because bash does not do proper division.
