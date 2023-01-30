#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <array>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <time.h>

//to compile  g++ -Wall -O3 hybrid_speciation_hard_selection.cpp -o fourloci_hs -lgsl -lgslcblas
using namespace std;

int read_parameter(const char* input_file, double& mean_off, double& s1, double& s2, double& s3, double& s4, double& e12, double& e13, double& e14, double& e23, double& e24, double& e34, double& r12, double& r23, double& r34, int& h12, int& h13, int& h14, int& h23, int& h24, int& h34, int& n_gen,  int& n_Ind, int& K, array<int, 16>& n_ind_hap);

void run(gsl_rng* r, int& T_locus1, int& T_locus2, int& T_locus3, int& T_locus4, bool& a1, bool& a2, bool& a3, bool& a4,const int& n_gen, const int& n_ind, const int& K, const double& r12, const double& r23, const double& r34, const array<double, 256>& fitness_table, const array<int, 16>& initial_hap_nb, int& n_final, int& min_pop, int& N1, int& N2, int& N3, int& N4);

void run(gsl_rng* r, int& T_locus1, int& T_locus2, int& T_locus3, int& T_locus4, bool& a1, bool& a2, bool& a3, bool& a4,const int& n_gen, const int& n_ind, const int& K, const double& r12, const double& r23, const double& r34, const array<double, 256>& fitness_table, const array<int, 16>& initial_hap_nb, int& n_final, const char* freq_file, const int& iterator, int& min_pop, int& N1, int& N2, int& N3, int& N4);

void gametes_prod(array<double, 16>& hap, array<double, 256>& gen, const double& rho12, const double& rho23, const double& rho34);

void track_alleles(const array<int, 256>& gen, const int& n_gen, const int& n_ind, int& T1, int& T2, int& T3, int& T4, bool& a1, bool& a2, bool& a3, bool& a4, int& N1, int& N2, int& N3, int& N4);

void track_alleles(const array<int, 256>& gen, const int& n_gen, const int& n_ind, int& T1, int& T2, int& T3, int& T4, bool& a1, bool& a2, bool& a3, bool& a4, int& N1, int& N2, int& N3, int& N4, int& locus1, int& locus2, int& locus3, int& locus4);

void test_fix(const int& n_gen, const int& n_ind, const int & locus, int & T, bool & a, int& N, const string& name);

void run(int& T_locus1,const int& n_gen, const int& n_ind, const int& K, const double& r12, const double& r23, const double& r34, const array<double, 256>& fitness_table, const array<int, 16>& initial_hap_nb, const char* freq_file);

void track_alleles(const array<double, 256>& gen, double& locus1, double& locus2, double& locus3, double& locus4);

int main(int argc, char *argv[])
{
    // check that the call of the program has the right number of parameters: 4, 5 or 6
    if (argc!=5 && argc!=6)
    {
        cout << "You need four mandatory parameters\n the parameter file name\t the output file, number of iterations and the seed\n the parameter file contains the different parameters in this order s1, s2, s3, s4, e12, e13, e14, e23, e24, e34, r12,r23,r34, h12, h13, h14, h23, h24, h34 \n There is one optional parameter, a file to save the frequency of the population at each generation. Currently, " << argc << "parameters are given.";
        return -1;
    }

    // declare variables:
    // selection: s1,s2,s3,s4 ; epistasis e12, e13, e14 e23 e24 e34; dominance of the epistasis h12, h13, h14, h23, h24, h34; recombination r12, r23, r34
    // initial frequency f1,...f16.

	//printf ("seed = %lu\n", gsl_rng_default_seed);


    double mean_off, s1, s2, s3, s4, e12, e13, e14, e23, e24, e34, r12,r23,r34, epis_12, epis_13, epis_14, epis_23, epis_24, epis_34;
    int h12, h13, h14, h23, h24, h34, indice, indicem, indicef, i1, i2,i3,i4, n_gen, n_ind, value, n_final, K, min_pop;

    array<double, 256> fitness_table; // create fitness array
	array<int, 16> initial_hap_nb; // initial frequencies array

    int T_locus1=0, T_locus2=0, T_locus3=0, T_locus4=0;
    bool allele1=0, allele2=0, allele3=0,allele4=0;
	int N1=0, N2=0, N3=0, N4=0;
    // read the parameter file
    read_parameter(argv[1], mean_off, s1, s2, s3, s4, e12, e13, e14, e23, e24, e34, r12,r23,r34,h12, h13, h14, h23, h24, h34,n_gen, n_ind, K, initial_hap_nb);

    // number of iterations of the same parameter set
	int nb_int = atoi(argv[3]);



	//    cout << "input read" << endl;
	// 	   cout << "simulating 4 loci with the following parameter"<<endl;
	//    cout << "selection: s1="<<s1 <<", s2"<<s2 <<", s3="<<s3 <<", s4="<<s4<<endl;
	//    cout << "epistasis: e12="<<e12 <<", e13="<< e13 <<", e14="<<e14 <<", e23="<<e23 <<", e24="<<e24 <<", e34="<< e34 <<endl;
	//    cout << "recombination: r12="<<r12 <<", r23="<<r23<<", r34="<<r34<<endl;
	//    cout << "dominance: h12="<<h12 <<", h13="<< h13 <<", h14="<<h14 <<", h23="<<h23 <<", h24="<<h24 <<", h34="<< h34 <<endl;
	//    cout << "number gen.= "<< n_gen << ", n_ind= " << n_ind <<", seed="<<seed<< endl;
	//    cout << "initial frequencies:"<<endl;
	//    for (int i1m=0; i1m<=1; i1m++){
	//        for (int i2m=0; i2m<=1; i2m++){
	//            for (int i3m=0; i3m<=1; i3m++ ){
	//                for (int i4m=0; i4m<=1; i4m++){
	//                    indicem =8*i1m+4*i2m+2*i3m+i4m;
	//                    cout << "hap"<<i1m <<i2m <<i3m<<i4m<< ": f=" <<initial_hap_nb[indicem]<<endl;
	//                }
	//            }
	//        }
	//    }

    // write the fitness table
    for (int i1m=0; i1m<=1; i1m++){
        for (int i2m=0; i2m<=1; i2m++){
            for (int i3m=0; i3m<=1; i3m++ ){
                for (int i4m=0; i4m<=1; i4m++){
                    indicem =8*i1m+4*i2m+2*i3m+i4m;
                    for (int i1f=0; i1f<=1; i1f++){
                        for (int i2f=0; i2f<=1; i2f++){
                            for (int i3f=0; i3f<=1; i3f++ ){
                                for (int i4f=0; i4f<=1; i4f++){
                                    indicef =8*i1f+4*i2f+2*i3f+i4f;
                                    indice  = 16*indicem+indicef;
                                    i1=i1m+i1f;i2=i2m+i2f;i3=i3m+i3f;i4=i4m+i4f;
                                    if(i1==1&&i2==1){epis_12=1+h12*e12 *i1*i2;}else{epis_12=pow(1+e12,i1*i2);}
                                    if(i1==1&&i3==1){epis_13=1+h13*e13 *i1*i3;}else{epis_13=pow(1+e13,i1*i3);}
                                    if(i1==1&&i4==1){epis_14=1+h14*e14 *i1*i4;}else{epis_14=pow(1+e14,i1*i4);}
                                    if(i2==1&&i3==1){epis_23=1+h23*e23 *i2*i3;}else{epis_23=pow(1+e23,i2*i3);}
                                    if(i2==1&&i4==1){epis_24=1+h24*e24 *i2*i4;}else{epis_24=pow(1+e24,i2*i4);}
                                    if(i3==1&&i4==1){epis_34=1+h34*e34 *i3*i4;}else{epis_34=pow(1+e34,i3*i4);}
                                    fitness_table[indice]=mean_off*pow(1+s1,(2-i1))*pow(1+s2,(2-i2))*pow(1+s3,(2-i3))*pow(1+s4,(2-i4))*epis_12*epis_13*epis_14*epis_23*epis_24*epis_34;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

	//	cout << "start sim"<<endl;

    // if four arguments have been given use the regular version of run

	if (nb_int>0){
		int seed = atoi(argv[4]);
		if (seed==-1){
			srand (time(NULL));
			seed=rand();
		}

	// initialize gsl random generator number for the multinomial sampling
		const gsl_rng_type *T;
		gsl_rng *r;
		gsl_rng_env_setup();

		T=gsl_rng_default;
		r=gsl_rng_alloc(T);
		gsl_rng_set(r, seed);

    // open output file
		FILE * fichierO;
		fichierO= fopen(argv[2],"a");
		if (fichierO == NULL) {
            cout << "wrong output file" << endl;
			return -1;
		}

    // print header of the output file. It contains all the parameters value of the run, to ensure repeatability.
		fprintf(fichierO, "Parameters: seed=%i, mean_off_nb=%f, s1=%f, s2=%f, s3=%f, s4=%f, e12=%f, e13=%f, e14=%f, e23=%f, e24=%f, e34=%f, r12=%f, r23=%f, r34=%f, h12=%i, h13=%i, h14=%i, h23=%i, h24=%i, h34=%i , n_gen=%i, n_ind=%i, K=%i ", seed, mean_off, s1, s2, s3, s4, e12, e13, e14, e23, e24, e34, r12,r23,r34, h12, h13, h14, h23, h24, h34, n_gen, n_ind, K);
		fprintf(fichierO, "initial freq: %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i\n", initial_hap_nb[0],initial_hap_nb[1],initial_hap_nb[2],initial_hap_nb[3],initial_hap_nb[4],initial_hap_nb[5],initial_hap_nb[6],initial_hap_nb[7],initial_hap_nb[8],initial_hap_nb[9],initial_hap_nb[10],initial_hap_nb[11],initial_hap_nb[12],initial_hap_nb[13],initial_hap_nb[14],initial_hap_nb[15]);

	// check that initial frequencies do sum up to one
		value=0;
		for (int i=0; i<16; i++){
			value=value+initial_hap_nb[i];
		}
		if (value!=n_ind){
			fprintf(fichierO, "wrong haplotype frequencies summing to: %i ", value);
			fprintf(fichierO, "difference is: %i", n_ind-value);
			fclose(fichierO);
		return -1;
		}
		
		for (int j=0; j<nb_int;j++){
	//			cout << "iteration "<< j<< "/"<<nb_int<< endl;
			T_locus1=0; T_locus2=0; T_locus3=0; T_locus4=0; // initialize fixation time
			allele1=0; allele2=0; allele3=0; allele4=0; // initialize allele fixation place holder
			N1=0, N2=0, N3=0, N4=0;
			// call main function
			
			min_pop=K;
			if (argc==5){
				run(r, T_locus1, T_locus2, T_locus3, T_locus4, allele1, allele2, allele3, allele4, n_gen, n_ind, K, r12, r23, r34, fitness_table, initial_hap_nb, n_final, min_pop, N1, N2, N3, N4);
			} else {
				run(r, T_locus1, T_locus2, T_locus3, T_locus4, allele1, allele2, allele3, allele4, n_gen, n_ind, K, r12, r23, r34, fitness_table, initial_hap_nb, n_final, argv[5], j, min_pop, N1, N2, N3, N4);
			} 
			fprintf(fichierO,"%i %i %i %i %i %i %i %i %i %i %i %i %i %i\n", allele1, T_locus1, allele2, T_locus2, allele3, T_locus3, allele4, T_locus4, n_final, min_pop, N1, N2, N3, N4);
		}			
		fclose(fichierO);
	} else {
		run(T_locus1, n_gen, n_ind, K, r12, r23, r34, fitness_table, initial_hap_nb, argv[5]);
	}


    return 1;
}


int read_parameter(const char* input_file, double& mean_off, double& s1, double& s2, double& s3, double& s4, double& e12, double& e13, double& e14, double& e23, double& e24, double& e34, double& r12, double& r23, double& r34, int& h12, int& h13, int& h14, int& h23, int& h24, int& h34, int& n_gen, int& n_Ind, int& K, array<int, 16>& n_ind_hap){

    int test;
    FILE * param_file;

    // read the parameter file; check that the different entries of the parameter file matches the expected type and store this value in the proper variable.

    param_file= fopen(input_file,"r");
    if (param_file == NULL) {
            cout << "wrong param file" << endl;
        return -1;
    }

    test=fscanf(param_file,"%lf ",&mean_off);
    if (test!=1){
        cout << "error param file"<<endl;
        return -1;
    }
    test=fscanf(param_file,"%lf ",&s1);
    if (test!=1){
        cout << "error param file"<<endl;
        return -1;
    }
    test=fscanf(param_file,"%lf ",&s2);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%lf ",&s3);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%lf ",&s4);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%lf ",&e12);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%lf ",&e13);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%lf ",&e14);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%lf ",&e23);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%lf ",&e24);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%lf ",&e34);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%lf ",&r12);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%lf ",&r23);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%lf ",&r34);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%i ",&h12);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%i ",&h13);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%i ",&h14);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%i ",&h23);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%i ",&h24);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%i ",&h34);
           if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%i ",&n_gen);
        if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%i ",&n_Ind);
        if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    test=fscanf(param_file,"%i ",&K);
        if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    for (int i=0; i<16; i++){
       test=fscanf(param_file, "%i ", &n_ind_hap[i]);
       if (test!=1){
            cout << "error param file"<<endl;
            return -1;
       }
    }
    fclose(param_file);
    return 1;
}

void run(gsl_rng* r, int& T_locus1, int& T_locus2, int& T_locus3, int& T_locus4, bool& a1, bool& a2, bool& a3, bool& a4,const int& n_gen, const int& n_ind, const int& K, const double& r12, const double& r23, const double& r34, const array<double, 256>& fitness_table, const array<int, 16>& initial_hap_nb, int& N_dyn_ind, int& min_pop, int& N1, int& N2, int& N3, int& N4){

	int indicem, indicef, indice, counter;
	


	array<double, 16> haplotype_freq; // define an array to store the frequency of all 16 gametes.
    array<double, 256> genotype_freq; // define an array to store the frequency of all 256 different genotypes.
    array<int, 256> genotype_ind_count;

    N_dyn_ind=n_ind;

    // calculate gamete frequencies
    for(int i=0; i< 16; i++){
    //	cout << "indice "<< i << " nb" << initial_hap_nb[i] << "\n";
    	haplotype_freq[i]=1.0*initial_hap_nb[i]/N_dyn_ind;
    }

    // initialize genotype array
    for (int i=0; i< 256; i++){
    	genotype_freq[i]=0.0;
    	genotype_ind_count[i]=0;
    }

    // form pure species parents in generation 0 - in the parameter files, the input is in terms of haplotype and not genotypes for easier usage but it corresponds to adults indiviudals nonetheless

    for (int i1m=0; i1m<=1; i1m++){
        for (int i2m=0; i2m<=1; i2m++){
            for (int i3m=0; i3m<=1; i3m++ ){
                for (int i4m=0; i4m<=1; i4m++){
                	indicem =8*i1m+4*i2m+2*i3m+i4m;
                	indice = 16*indicem+indicem;
                	genotype_freq[indice]=haplotype_freq[indicem];
    //            	cout << "indice " << indice << "gen freq  " << genotype_freq[indice] << " hap_freq " << haplotype_freq[indicem] << "\n";
                }
            }
        }
    }




    // life cycle
    for (int gen=0; gen<n_gen; gen++){


    	// determine number of offprings per genotype (selection only on female).
        counter=0;
        for (int i=0; i<256; i++){
            genotype_ind_count[i]=gsl_ran_poisson(r, N_dyn_ind*genotype_freq[i]*fitness_table[i]);
            //if (genotype_freq[i] >0){
                //cout << "gen "<< gen << "indice " << i << "gen freq  " << genotype_freq[i] << "fitness " << fitness_table[i]<< "param poisson " << N_dyn_ind*genotype_freq[i]*fitness_table[i] << "value " << genotype_ind_count[i] << "\n";
            //}
            counter=counter+genotype_ind_count[i];
        }
        N_dyn_ind=counter;

        for (int i=0; i<256; i++){
            genotype_freq[i]=1.0*genotype_ind_count[i]/N_dyn_ind;
        }



        // check whether any locus got fixed for a specific allele
        track_alleles(genotype_ind_count, (gen-1), N_dyn_ind, T_locus1, T_locus2, T_locus3, T_locus4, a1, a2, a3, a4, N1, N2, N3, N4);
		
		if (N_dyn_ind<min_pop){
			min_pop=N_dyn_ind;
		}
        //cout<<"time " << T_locus1 << "  "<< T_locus2 << "  "<< T_locus3 << "  "<< T_locus4 << "  " <<endl;
        // if all loci are no longer polymorphic stop the iteration process
        if((T_locus1>0&& T_locus2>0&& T_locus3>0&& T_locus4>0 &&N_dyn_ind>K)||(N_dyn_ind==0)){

            break;
        }


		
        // Population regulation
        if (N_dyn_ind>K){
            N_dyn_ind=K;
        }

        // empty the gamete vector.
        for (int i=0; i<16; i++){
            haplotype_freq[i]=0;
        }

    	// calculate the male gamete frequency
    	gametes_prod(haplotype_freq,genotype_freq,r12,r23,r34);


    	// empty the genotype frequency vector.
		for (int i=0; i<256; i++){
    		genotype_freq[i]=0;
    	}

    	// calculate genotype freq of next generation.
        for (int i1m=0; i1m<=1; i1m++){
            for (int i2m=0; i2m<=1; i2m++){
                for (int i3m=0; i3m<=1; i3m++ ){
                    for (int i4m=0; i4m<=1; i4m++){
                        indicem =8*i1m+4*i2m+2*i3m+i4m;

                        for (int i1f=0; i1f<=1; i1f++){
                            for (int i2f=0; i2f<=1; i2f++){
                                for (int i3f=0; i3f<=1; i3f++ ){
                                    for (int i4f=0; i4f<=1; i4f++){
                                        indicef =8*i1f+4*i2f+2*i3f+i4f;
                                        indice = 16*indicem+indicef;
                                        genotype_freq[indice]=haplotype_freq[indicem]*haplotype_freq[indicef];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
	}
}

void gametes_prod(array<double, 16>& hap, array<double, 256>& gen, const double& rho12, const double& rho23, const double& rho34){
    // define the gamete production by considering the allele at each locus from the paternal (i.m) and maternal (i.f) origin
    int indice, indicef, indicem, i1, i2,i3, i4;
    for (int i1m=0; i1m<=1; i1m++){
        for (int i2m=0; i2m<=1; i2m++){
            for (int i3m=0; i3m<=1; i3m++ ){
                for (int i4m=0; i4m<=1; i4m++){
                    for (int i1f=0; i1f<=1; i1f++){
                        for (int i2f=0; i2f<=1; i2f++){
                            for (int i3f=0; i3f<=1; i3f++ ){
                                for (int i4f=0; i4f<=1; i4f++){
                                    indicem =8*i1m+4*i2m+2*i3m+i4m; // define the paternal haplotype
                                    indicef =8*i1f+4*i2f+2*i3f+i4f; // define the maternal haplotype
                                    indice  = 16*indicem+indicef;   // define the genotype
                                    i1=i1m+i1f;i2=i2m+i2f;i3=i3m+i3f;i4=i4m+i4f; // define the status at each locus: 0 for 0/0, 1 for 0/1 or 1/0 and 2 for 1/1
                                    // check that there is individuals of this genotype
                                    if (gen[indice]>0){
                                    // recombination matters only if there is at least two loci heterozygotes
                                    // first assume that first locus is heterozygote
                                    if(i1==1&&(i2==1||i3==1||i4==1)){
                                        // assume that the third and fourth loci are homozygote. Recombination only matter between loci 1 and 2
                                        if (i4!=1&&i3!=1){
                                            hap[indicem]=hap[indicem]+ (1-rho12)*gen[indice]/2; // add half the non recombining frequency of the genotype to the haplotype matching the paternal one
                                            hap[indicef]=hap[indicef]+ (1-rho12)*gen[indice]/2; // add half the non recombining frequency of the genotype to the haplotype matching the maternal one
                                            // calculate the indices of the two recombinant haplotypes; recombination only matter if it happens between the first and second loci
                                            indicem =8*i1m+4*i2f+2*i3f+i4f;
                                            indicef =8*i1f+4*i2m+2*i3m+i4m;
                                            hap[indicem]=hap[indicem]+ rho12*gen[indice]/2; // add half of the frequency of the recombining genotype to the "first" haplotype
                                            hap[indicef]=hap[indicef]+ rho12*gen[indice]/2; // add half of the frequency of the recombining genotype to the "second" haplotype
                                            //cout << indice  << " rec 12" << endl;
                                        // assume that the fourth locus is homozygote
                                        } else if(i4!=1){
                                            // add the frequency of non-recombining genotypes.
                                            hap[indicem]=hap[indicem]+ (1-rho12)*(1-rho23)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ (1-rho12)*(1-rho23)*gen[indice]/2;
                                            // recombination happens between loci 1 and 2 only
                                            indicem =8*i1m+4*i2f+2*i3f+i4f;
                                            indicef =8*i1f+4*i2m+2*i3m+i4m;
                                            hap[indicem]=hap[indicem]+ rho12*(1-rho23)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho12*(1-rho23)*gen[indice]/2;
                                            // recombination happens between loci 2 and 3 only
                                            indicem =8*i1m+4*i2m+2*i3f+i4f;
                                            indicef =8*i1f+4*i2f+2*i3m+i4m;
                                            hap[indicem]=hap[indicem]+ rho23*(1-rho12)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho23*(1-rho12)*gen[indice]/2;
                                            // recombination happens between loci 1 and 2 and between loci 2 and 3
                                            indicem =8*i1m+4*i2f+2*i3m+i4m;
                                            indicef =8*i1f+4*i2m+2*i3f+i4f;
                                            hap[indicem]=hap[indicem]+ rho23*rho12*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho23*rho12*gen[indice]/2;
                                            //cout << indice  << " rec 12 and 23" << endl;
                                        // assume that the second, third  and fourth loci can be heterozygote
                                        } else{
                                            // add the frequency of non-recombining genotypes.
                                            hap[indicem]=hap[indicem]+ (1-rho12)*(1-rho23)*(1-rho34)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ (1-rho12)*(1-rho23)*(1-rho34)*gen[indice]/2;
                                            // recombination happens between loci 1 and 2 only
                                            indicem =8*i1m+4*i2f+2*i3f+i4f;
                                            indicef =8*i1f+4*i2m+2*i3m+i4m;
                                            hap[indicem]=hap[indicem]+ rho12*(1-rho23)*(1-rho34)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho12*(1-rho23)*(1-rho34)*gen[indice]/2;
                                            // recombination happens between loci 2 and 3 only
                                            indicem =8*i1m+4*i2m+2*i3f+i4f;
                                            indicef =8*i1f+4*i2f+2*i3m+i4m;
                                            hap[indicem]=hap[indicem]+ rho23*(1-rho12)*(1-rho34)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho23*(1-rho12)*(1-rho34)*gen[indice]/2;
                                            // recombination happens between loci 1 and 2 and between loci 2 and 3
                                            indicem =8*i1m+4*i2f+2*i3m+i4m;
                                            indicef =8*i1f+4*i2m+2*i3f+i4f;
                                            hap[indicem]=hap[indicem]+ rho23*rho12*(1-rho34)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho23*rho12*(1-rho34)*gen[indice]/2;
                                            // recombination happens between loci 3 and 4 only
                                            indicem =8*i1m+4*i2m+2*i3m+i4f;
                                            indicef =8*i1f+4*i2f+2*i3f+i4m;
                                            hap[indicem]=hap[indicem]+ rho34*(1-rho23)*(1-rho12)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho34*(1-rho23)*(1-rho12)*gen[indice]/2;
                                            // recombination happens between loci 1 and 2 and between loci 3 and 4
                                            indicem =8*i1m+4*i2f+2*i3f+i4m;
                                            indicef =8*i1f+4*i2m+2*i3m+i4f;
                                            hap[indicem]=hap[indicem]+ rho12*(1-rho23)*rho34*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho12*(1-rho23)*rho34*gen[indice]/2;
                                            // recombination happens between loci 2 and 3 and between loci 3 and 4
                                            indicem =8*i1m+4*i2m+2*i3f+i4m;
                                            indicef =8*i1f+4*i2f+2*i3m+i4f;
                                            hap[indicem]=hap[indicem]+ rho23*(1-rho12)*rho34*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho23*(1-rho12)*rho34*gen[indice]/2;
                                            // recombination happens between loci 1 and 2, between loci 2 and 3 and between loci 3 and 4
                                            indicem =8*i1m+4*i2f+2*i3m+i4f;
                                            indicef =8*i1f+4*i2m+2*i3f+i4m;
                                            hap[indicem]=hap[indicem]+ rho12*rho23*rho34*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho12*rho23*rho34*gen[indice]/2;
                                            //cout << indice  << " rec 12, 23 and 34 " << endl;
                                        }
                                    // assume that locus 1 is homozygote and locus 2 is heterozygote
                                    } else if (i2==1&&(i3==1||i4==1)){
                                        // assume that locus 3 is heterozygote
                                        if (i4!=1){
                                            // add the frequency of non-recombining genotypes.
                                            hap[indicem]=hap[indicem]+ (1-rho23)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ (1-rho23)*gen[indice]/2;
                                            // recombination happens between loci 2 and 3
                                            indicem =8*i1m+4*i2m+2*i3f+i4f;
                                            indicef =8*i1f+4*i2f+2*i3m+i4m;
                                            hap[indicem]=hap[indicem]+ rho23*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho23*gen[indice]/2;
                                            //cout << indice  << " rec 23" << endl;
                                        // assume that both locus 3 and 4 can be heterozygote
                                        } else{
                                            // add the frequency of non-recombining genotypes.
                                            hap[indicem]=hap[indicem]+ (1-rho34)*(1-rho23)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ (1-rho34)*(1-rho23)*gen[indice]/2;
                                            // recombination happens between loci 2 and 3
                                            indicem =8*i1m+4*i2m+2*i3f+i4f;
                                            indicef =8*i1f+4*i2f+2*i3m+i4m;
                                            hap[indicem]=hap[indicem]+ rho23*(1-rho34)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho23*(1-rho34)*gen[indice]/2;
                                            // recombination happens between loci 3 and 4
                                            indicem =8*i1m+4*i2m+2*i3m+i4f;
                                            indicef =8*i1f+4*i2f+2*i3f+i4m;
                                            hap[indicem]=hap[indicem]+ rho34*(1-rho23)*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho34*(1-rho23)*gen[indice]/2;
                                            // recombination happens between loci 2 and 3 and loci 3 and 4
                                            indicem =8*i1m+4*i2m+2*i3f+i4m;
                                            indicef =8*i1f+4*i2f+2*i3m+i4f;
                                            hap[indicem]=hap[indicem]+ rho23*rho34*gen[indice]/2;
                                            hap[indicef]=hap[indicef]+ rho23*rho34*gen[indice]/2;
                                            //cout << indice  << " rec 23 and 34" << endl;
                                        }
                                    // loci 1 and 2 are homozygote. Therefore loci 3 and 4  have to be heterozygote.
                                    } else if(i3==1&&i4==1){
                                        // add the frequency of non-recombining genotypes.
                                        hap[indicem]=hap[indicem]+ (1-rho34)*gen[indice]/2;
                                        hap[indicef]=hap[indicef]+ (1-rho34)*gen[indice]/2;
                                        // recombination happens between loci 3 and 4
                                        indicem =8*i1m+4*i2m+2*i3m+i4f;
                                        indicef =8*i1f+4*i2f+2*i3f+i4m;
                                        hap[indicem]=hap[indicem]+ rho34*gen[indice]/2;
                                        hap[indicef]=hap[indicef]+ rho34*gen[indice]/2;
                                        //cout << indice  << " rec 34" << endl;
                                    } else{
                                        // // add the frequency of non-recombining genotypes.
                                        hap[indicem]=hap[indicem]+ gen[indice]/2;
                                        hap[indicef]=hap[indicef]+  gen[indice]/2;
                                        //cout << indice  << "no rec" << endl;
                                    }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void track_alleles(const array<int, 256>& gen, const int& n_gen, const int& n_ind, int& T1, int& T2, int& T3, int& T4, bool& a1, bool& a2, bool& a3, bool& a4, int& N1, int& N2, int& N3, int& N4){
    int locus1=0, locus2=0, locus3=0, locus4=0;
    int indicem, indicef, indice;
    // calculate the frequency of allele 1 at each locus
    for (int i1m=0; i1m<=1; i1m++){
        for (int i2m=0; i2m<=1; i2m++){
            for (int i3m=0; i3m<=1; i3m++ ){
                for (int i4m=0; i4m<=1; i4m++){
                    indicem =8*i1m+4*i2m+2*i3m+i4m;

                    for (int i1f=0; i1f<=1; i1f++){
                        for (int i2f=0; i2f<=1; i2f++){
                            for (int i3f=0; i3f<=1; i3f++ ){
                                for (int i4f=0; i4f<=1; i4f++){
                                    indicef =8*i1f+4*i2f+2*i3f+i4f;
                                    indice = 16*indicem+indicef;
                                    locus1=locus1+(i1m+i1f)*gen[indice];
                                    locus2=locus2+(i2m+i2f)*gen[indice];
                                    locus3=locus3+(i3m+i3f)*gen[indice];
                                    locus4=locus4+(i4m+i4f)*gen[indice];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //cout << n_ind << " " << locus1 << " " << locus2 << " " << locus3 << " " << locus4 <<endl;
    // check for fixation of said allele, and when it happened.
    test_fix(n_gen, n_ind, locus1, T1, a1, N1, "locus 1");
    test_fix(n_gen, n_ind, locus2, T2, a2, N2, "locus 2");
    test_fix(n_gen, n_ind, locus3, T3, a3, N3, "locus 3");
    test_fix(n_gen, n_ind, locus4, T4, a4, N4, "locus 4");
}

void test_fix(const int& n_gen, const int& n_ind, const int & locus, int & T, bool & a, int& N , const string& name){
    // for a given locus, t record the time to fixation of one of the two alleles 0 or 1, and a save which allele fixed.
    if (T==0){
        if (locus==0){
            T=n_gen;
			N=n_ind;
	            //cout << "fixation "<< name << " allele 0 gen. "<< n_gen << endl;
        } else if(locus==2*n_ind){
            T=n_gen;
            a=1;
			N=n_ind;
	            //cout << "fixation "<< name << " allele 1 gen. "<< n_gen <<endl;
        }
    }
}

// tracking every generation

void run(gsl_rng* r, int& T_locus1, int& T_locus2, int& T_locus3, int& T_locus4, bool& a1, bool& a2, bool& a3, bool& a4,const int& n_gen, const int& n_ind, const int& K, const double& r12, const double& r23, const double& r34, const array<double, 256>& fitness_table, const array<int, 16>& initial_hap_nb, int& N_dyn_ind, const char* freq_file, const int& iterator, int& min_pop, int& N1, int& N2, int& N3, int& N4){

    int indicem, indicef, indice, counter, locus1, locus2, locus3, locus4;
    


    array<double, 16> haplotype_freq; // define an array to store the frequency of all 16 gametes.
    array<double, 256> genotype_freq; // define an array to store the frequency of all 256 different genotypes.
    array<int, 256> genotype_ind_count;

    N_dyn_ind=n_ind;

    // calculate gamete frequencies
    for(int i=0; i< 16; i++){
    //  cout << "indice "<< i << " nb" << initial_hap_nb[i] << "\n";
        haplotype_freq[i]=1.0*initial_hap_nb[i]/N_dyn_ind;
    }

    // initialize genotype array
    for (int i=0; i< 256; i++){
        genotype_freq[i]=0.0;
        genotype_ind_count[i]=0;
    }

    // form pure species parents in generation 0 - in the parameter files, the input is in terms of haplotype and not genotypes for easier usage but it corresponds to adults indiviudals nonetheless

    for (int i1m=0; i1m<=1; i1m++){
        for (int i2m=0; i2m<=1; i2m++){
            for (int i3m=0; i3m<=1; i3m++ ){
                for (int i4m=0; i4m<=1; i4m++){
                    indicem =8*i1m+4*i2m+2*i3m+i4m;
                    indice = 16*indicem+indicem;
                    genotype_freq[indice]=haplotype_freq[indicem];
    //              cout << "indice " << indice << "gen freq  " << genotype_freq[indice] << " hap_freq " << haplotype_freq[indicem] << "\n";
                }
            }
        }
    }

    FILE * frequency_file;
	
	string name_file= freq_file;
	unsigned int find_ext=0;
	
	find_ext=name_file.find(".");
	
	if ( find_ext!= string::npos){
		name_file=name_file.substr(0, find_ext)+to_string(iterator)+name_file.substr(find_ext);
		cout << name_file <<endl;
	} else {
		cout << "Frequency file has no extension \n all outcomes will be append to the same file" << endl;;
	} 
	
	
    // read the parameter file; check that the different entries of the parameter file matches the expected type and store this value in the proper variable.

    frequency_file= fopen(name_file.c_str(),"a");


    // life cycle
    for (int gen=0; gen<n_gen; gen++){


        // determine number of offprings per genotype (selection only on female).
        counter=0;
        for (int i=0; i<256; i++){
            genotype_ind_count[i]=gsl_ran_poisson(r, N_dyn_ind*genotype_freq[i]*fitness_table[i]);
            //if (genotype_freq[i] >0){
                //cout << "gen "<< gen << "indice " << i << "gen freq  " << genotype_freq[i] << "fitness " << fitness_table[i]<< "param poisson " << N_dyn_ind*genotype_freq[i]*fitness_table[i] << "value " << genotype_ind_count[i] << "\n";
            //}
            counter=counter+genotype_ind_count[i];
        }
        N_dyn_ind=counter;

        for (int i=0; i<256; i++){
            genotype_freq[i]=1.0*genotype_ind_count[i]/N_dyn_ind;
        }

		if (N_dyn_ind<min_pop){
			min_pop=N_dyn_ind;
		}

        // check whether any locus got fixed for a specific allele
        track_alleles(genotype_ind_count, (gen-1), N_dyn_ind, T_locus1, T_locus2, T_locus3, T_locus4, a1, a2, a3, a4, N1, N2, N3, N4, locus1, locus2, locus3, locus4);
        fprintf(frequency_file, "%i %i %i %i %i\n", N_dyn_ind, locus1, locus2, locus3, locus4);
        //cout<<"time " << T_locus1 << "  "<< T_locus2 << "  "<< T_locus3 << "  "<< T_locus4 << "  " <<endl;
        // if all loci are no longer polymorphic stop the iteration process
        if((T_locus1>0&& T_locus2>0&& T_locus3>0&& T_locus4>0&&N_dyn_ind>K)||(N_dyn_ind==0)){
            fclose(frequency_file);
            break;
        }

        // Population regulation
        if (N_dyn_ind>K){
            N_dyn_ind=K;
        }

        // empty the gamete vector.
        for (int i=0; i<16; i++){
            haplotype_freq[i]=0;
        }

        // calculate the male gamete frequency
        gametes_prod(haplotype_freq,genotype_freq,r12,r23,r34);


        // empty the genotype frequency vector.
        for (int i=0; i<256; i++){
            genotype_freq[i]=0;
        }

        // calculate genotype freq of next generation.
        for (int i1m=0; i1m<=1; i1m++){
            for (int i2m=0; i2m<=1; i2m++){
                for (int i3m=0; i3m<=1; i3m++ ){
                    for (int i4m=0; i4m<=1; i4m++){
                        indicem =8*i1m+4*i2m+2*i3m+i4m;

                        for (int i1f=0; i1f<=1; i1f++){
                            for (int i2f=0; i2f<=1; i2f++){
                                for (int i3f=0; i3f<=1; i3f++ ){
                                    for (int i4f=0; i4f<=1; i4f++){
                                        indicef =8*i1f+4*i2f+2*i3f+i4f;
                                        indice = 16*indicem+indicef;
                                        genotype_freq[indice]=haplotype_freq[indicem]*haplotype_freq[indicef];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    fclose(frequency_file);
}

void track_alleles(const array<int, 256>& gen, const int& n_gen, const int& n_ind, int& T1, int& T2, int& T3, int& T4, bool& a1, bool& a2, bool& a3, bool& a4, int& N1, int& N2, int& N3, int& N4, int& locus1, int& locus2, int& locus3, int& locus4){
    locus1=0;
    locus2=0;
    locus3=0;
    locus4=0;
    int indicem, indicef, indice;
    // calculate the frequency of allele 1 at each locus
    for (int i1m=0; i1m<=1; i1m++){
        for (int i2m=0; i2m<=1; i2m++){
            for (int i3m=0; i3m<=1; i3m++ ){
                for (int i4m=0; i4m<=1; i4m++){
                    indicem =8*i1m+4*i2m+2*i3m+i4m;

                    for (int i1f=0; i1f<=1; i1f++){
                        for (int i2f=0; i2f<=1; i2f++){
                            for (int i3f=0; i3f<=1; i3f++ ){
                                for (int i4f=0; i4f<=1; i4f++){
                                    indicef =8*i1f+4*i2f+2*i3f+i4f;
                                    indice = 16*indicem+indicef;
                                    locus1=locus1+(i1m+i1f)*gen[indice];
                                    locus2=locus2+(i2m+i2f)*gen[indice];
                                    locus3=locus3+(i3m+i3f)*gen[indice];
                                    locus4=locus4+(i4m+i4f)*gen[indice];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //cout << n_ind << " " << locus1 << " " << locus2 << " " << locus3 << " " << locus4 <<endl;
    // check for fixation of said allele, and when it happened.
    test_fix(n_gen, n_ind, locus1, T1, a1, N1, "locus 1");
    test_fix(n_gen, n_ind, locus2, T2, a2, N2, "locus 2");
    test_fix(n_gen, n_ind, locus3, T3, a3, N3, "locus 3");
    test_fix(n_gen, n_ind, locus4, T4, a4, N4, "locus 4");
}

// removing drift

void run(int& T_locus1,const int& n_gen, const int& n_ind, const int& K, const double& r12, const double& r23, const double& r34, const array<double, 256>& fitness_table, const array<int, 16>& initial_hap_nb, const char* freq_file){

	int indicem, indicef, indice;
	double counter, N_dyn_ind, locus1, locus2, locus3, locus4; 


	array<double, 16> haplotype_freq; // define an array to store the frequency of all 16 gametes.
    array<double, 256> genotype_freq; // define an array to store the frequency of all 256 different genotypes.
    array<double, 256> genotype_ind_count;

    N_dyn_ind=n_ind;

    // calculate gamete frequencies
    for(int i=0; i< 16; i++){
    //	cout << "indice "<< i << " nb" << initial_hap_nb[i] << "\n";
    	haplotype_freq[i]=1.0*initial_hap_nb[i]/N_dyn_ind;
    }

    // initialize genotype array
    for (int i=0; i< 256; i++){
    	genotype_freq[i]=0.0;
    	genotype_ind_count[i]=0.0;
    }

    // form pure species parents in generation 0 - in the parameter files, the input is in terms of haplotype and not genotypes for easier usage but it corresponds to adults indiviudals nonetheless

    for (int i1m=0; i1m<=1; i1m++){
        for (int i2m=0; i2m<=1; i2m++){
            for (int i3m=0; i3m<=1; i3m++ ){
                for (int i4m=0; i4m<=1; i4m++){
                	indicem =8*i1m+4*i2m+2*i3m+i4m;
                	indice = 16*indicem+indicem;
                	genotype_freq[indice]=haplotype_freq[indicem];
                	//cout << "indice " << indice << "gen freq  " << genotype_freq[indice] << " hap_freq " << haplotype_freq[indicem] << "\n";
                }
            }
        }
    }
	
	
	
	for (int i=0; i<16;i++){
		for (int j=0; j<16; j++){
			cout << fitness_table[16*i+j] << "\t";
		}
		cout <<endl;
	}

	FILE * frequency_file;
	frequency_file= fopen(freq_file,"a");

	track_alleles(genotype_freq, locus1, locus2, locus3, locus4);
	locus1=locus1/2;
	locus2=locus2/2;
	locus3=locus3/2;
	locus4=locus4/2;
	fprintf(frequency_file, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", N_dyn_ind, locus1, locus2, locus3, locus4, haplotype_freq[0], haplotype_freq[1], haplotype_freq[2], haplotype_freq[3], haplotype_freq[4], haplotype_freq[5], haplotype_freq[6], haplotype_freq[7], haplotype_freq[8], haplotype_freq[9], haplotype_freq[10], haplotype_freq[11], haplotype_freq[12], haplotype_freq[13], haplotype_freq[14], haplotype_freq[15]);
	
    // life cycle
    for (int gen=0; gen<n_gen; gen++){


    	// determine number of offprings per genotype.
        counter=0.0;
        for (int i=0; i<256; i++){
            genotype_ind_count[i]=N_dyn_ind*genotype_freq[i]*fitness_table[i];
            if (genotype_freq[i] >0){
                cout << "gen "<< gen << "\t indice " << i << "\t gen freq  " << genotype_freq[i] << "\t fitness " << fitness_table[i]<< "\t param poisson " << N_dyn_ind*genotype_freq[i]*fitness_table[i] << "\t value " << genotype_ind_count[i] << "\n";
            }
            counter=counter+genotype_ind_count[i];
        }
        N_dyn_ind=counter;

        for (int i=0; i<256; i++){
            genotype_freq[i]=1.0*genotype_ind_count[i]/N_dyn_ind;
        }




		
		

        // empty the gamete vector.
        for (int i=0; i<16; i++){
            haplotype_freq[i]=0;
        }

    	// calculate the gamete frequency
    	gametes_prod(haplotype_freq,genotype_freq,r12,r23,r34);


    	// empty the genotype frequency vector.
		for (int i=0; i<256; i++){
    		genotype_freq[i]=0;
    	}

    	// calculate genotype freq of next generation.
        for (int i1m=0; i1m<=1; i1m++){
            for (int i2m=0; i2m<=1; i2m++){
                for (int i3m=0; i3m<=1; i3m++ ){
                    for (int i4m=0; i4m<=1; i4m++){
                        indicem =8*i1m+4*i2m+2*i3m+i4m;

                        for (int i1f=0; i1f<=1; i1f++){
                            for (int i2f=0; i2f<=1; i2f++){
                                for (int i3f=0; i3f<=1; i3f++ ){
                                    for (int i4f=0; i4f<=1; i4f++){
                                        indicef =8*i1f+4*i2f+2*i3f+i4f;
                                        indice = 16*indicem+indicef;
                                        genotype_freq[indice]=haplotype_freq[indicem]*haplotype_freq[indicef];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
		track_alleles(genotype_freq, locus1, locus2, locus3, locus4);
		locus1=locus1/2;
		locus2=locus2/2;
		locus3=locus3/2;
		locus4=locus4/2;
		fprintf(frequency_file, "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", N_dyn_ind, locus1, locus2, locus3, locus4, haplotype_freq[0], haplotype_freq[1], haplotype_freq[2], haplotype_freq[3], haplotype_freq[4], haplotype_freq[5], haplotype_freq[6], haplotype_freq[7], haplotype_freq[8], haplotype_freq[9], haplotype_freq[10], haplotype_freq[11], haplotype_freq[12], haplotype_freq[13], haplotype_freq[14], haplotype_freq[15]);

		//cout<<"time " << T_locus1 << "  "<< T_locus2 << "  "<< T_locus3 << "  "<< T_locus4 << "  " <<endl;
		// if all loci are no longer polymorphic stop the iteration process
		if((N_dyn_ind>K)||(N_dyn_ind<1)){
			break;
		}
	}
	        // check whether any locus got fixed for a specific allele
    
}

void track_alleles(const array<double, 256>& gen, double& locus1, double& locus2, double& locus3, double& locus4){
    locus1=0.0;
    locus2=0.0;
    locus3=0.0;
    locus4=0.0;
    int indicem, indicef, indice;
    // calculate the frequency of allele 1 at each locus
    for (int i1m=0; i1m<=1; i1m++){
        for (int i2m=0; i2m<=1; i2m++){
            for (int i3m=0; i3m<=1; i3m++ ){
                for (int i4m=0; i4m<=1; i4m++){
                    indicem =8*i1m+4*i2m+2*i3m+i4m;

                    for (int i1f=0; i1f<=1; i1f++){
                        for (int i2f=0; i2f<=1; i2f++){
                            for (int i3f=0; i3f<=1; i3f++ ){
                                for (int i4f=0; i4f<=1; i4f++){
                                    indicef =8*i1f+4*i2f+2*i3f+i4f;
                                    indice = 16*indicem+indicef;
                                    locus1=locus1+(i1m+i1f)*gen[indice];
                                    locus2=locus2+(i2m+i2f)*gen[indice];
                                    locus3=locus3+(i3m+i3f)*gen[indice];
                                    locus4=locus4+(i4m+i4f)*gen[indice];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

