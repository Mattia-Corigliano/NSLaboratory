/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Studente: Mattia Corigliano, Matr. 944964
_/    _/  _/_/_/  _/_/_/_/ email: mattia.corigliano@studenti.unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm> 
#include "random.h"

using namespace std;

double integrand(double );
double distribution(double ); //for importance sampling
double sd(vector<double>, vector<double>, int);
 
int main (){
	//*****initializing random numbers generation
	Random rnd;
	int seed[4];
   	int p1, p2;
   	ifstream Primes("Primes");
   	if (Primes.is_open()){
      		Primes >> p1 >> p2 ;
   	} else cerr << "PROBLEM: Unable to open Primes" << endl;
   	Primes.close();
   	ifstream input("seed.in");
   	string property;
   	if (input.is_open()){
      		while ( !input.eof() ){
         		input >> property;
         		if( property == "RANDOMSEED" ){
            		input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            		rnd.SetRandom(seed,p1,p2);
         		}
      		}
      		input.close();
   	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
	rnd.SaveSeed();
        //*******

	//initializing output files
	ofstream output("output.txt");
	ofstream output2("output2.txt");
	ofstream output3("output3.txt"); //for visualizing the walkers's paths
	ofstream output4("output4.txt"); //for the statistical analysis
	ofstream output5("output5.txt"); //for the statistical analysis
//***********************************ESERCIZIO 2.1*********************************************************
	int M = 1E6;
	int N_blocks = 1E2;
	int M_xblock = M/N_blocks;	
	
	vector<double>random_number(M, 0.);
	vector<double>integral(N_blocks, 0.);
	vector<double>mean(N_blocks, 0.);
	vector<double>mean_squared(N_blocks, 0.);
	vector<double>error(N_blocks, 0.);

	//genero i numeri random
	for(int i=0; i<M; i++)
		random_number[i]=rnd.Rannyu();

	//carico le medie in ogni blocco
	for(int block=0; block<N_blocks; block++){
		double mean = 0;
		for(int points=0; points<M_xblock; points++)
			mean = mean + integrand(random_number[points+M_xblock*block])/M_xblock;
		integral[block] = mean;
	}
	for(int block=0; block<N_blocks; block++){
		for(int i=0; i<block+1; i++){
			mean[block] = mean[block] + integral[i];
			mean_squared[block]=mean_squared[block]+integral[i]*integral[i];
		}
		mean[block]=mean[block]/(block+1);
		mean_squared[block]=mean_squared[block]/(block+1);
		error[block]=sd(mean_squared, mean, block+1);
		output << double((block+1)*M_xblock) << setw(15) << mean[block] << setw(15) << error[block] << endl;
	}
	for(int i=0; i<N_blocks; i++){
		integral[i]=0.;
		mean[i]=0.;
		mean_squared[i]=0.;
		error[i]=0.;
	}
	//carico le medie in ogni blocco
	for(int block=0; block<N_blocks; block++){
		double mean = 0;
		for(int points=0; points<M_xblock; points++){
			double x = distribution(random_number[points+M_xblock*block]);
			double temp = integrand(x)/(2*(1-x));
			mean = mean + temp/M_xblock;
		}
		integral[block] = mean;
	}
	for(int block=0; block<N_blocks; block++){
		for(int i=0; i<block+1; i++){
			mean[block] = mean[block] + integral[i];
			mean_squared[block]=mean_squared[block]+integral[i]*integral[i];
		}
		mean[block]=mean[block]/(block+1);
		mean_squared[block]=mean_squared[block]/(block+1);
		error[block]=sd(mean_squared, mean, block+1);
		output2 << double((block+1)*M_xblock) << setw(15) << mean[block] << setw(15) << error[block] << endl;
	}

//************************************3D Random walk*********************************************
	int n_walkers = 1E4; //number of different realizations of the process, i.e. walkers	
	int n_steps = 1E2; //number of steps taken by each walker
	double position_step=1.;
	vector<double> position_lattice(3*(n_steps+1)*n_walkers, 0.);
	vector<double> position_continuum(3*(n_steps+1)*n_walkers, 0.);
	for(int walker=0; walker<n_walkers; walker++){
		int x=0, y=0, z=0;
		double x2=0, y2=0, z2=0;
		for(int step=1; step<n_steps+1; step++){
			int direction = int(rnd.Rannyu(0., 3.)); //simulate uniformly a direction (x, y, z)--> 0, 1, 2
			double verse = rnd.Rannyu();
			if(direction == 0){
				if(verse<0.5)
					x=x-position_step;
				else
					x=x+position_step;
			}
			if(direction == 1){
				if(verse<0.5)
					y=y-position_step;
				else
					y=y+position_step;
			}
			if(direction == 2){
				if(verse<0.5)
					z=z-position_step;
				else
					z=z+position_step;
			}
			position_lattice[3*step+3*(n_steps+1)*walker] = x;
			position_lattice[3*step+1+3*(n_steps+1)*walker]=y;
			position_lattice[3*step+2+3*(n_steps+1)*walker]=z;
			double phi = rnd.Rannyu(0., 2.*M_PI);
			double theta = acos(1.-2.*rnd.Rannyu());
			x2 = x2 + position_step*sin(theta)*cos(phi);
			y2 = y2 + position_step*sin(theta)*sin(phi);
			z2 = z2 + position_step*cos(theta);
			position_continuum[3*step+3*(n_steps+1)*walker] = x2;
			position_continuum[3*step+1+3*(n_steps+1)*walker]=y2;
			position_continuum[3*step+2+3*(n_steps+1)*walker]=z2;
		}
	}
	for(int i=0; i<int(position_lattice.size()); i++){
		output3<<position_lattice[i]<<endl;
		output4<<position_continuum[i]<<endl;
	}

	vector<double> squared_distance(n_walkers*(n_steps+1), 0.);
	vector<double> squared_distance2(n_walkers*(n_steps+1), 0.);
	for(int i=0; i<int(position_lattice.size()); i=i+3){
		squared_distance[i/3] = position_lattice[i]*position_lattice[i]+position_lattice[i+1]*position_lattice[i+1]+position_lattice[i+2]*position_lattice[i+2];
		squared_distance2[i/3] = position_continuum[i]*position_continuum[i]+position_continuum[i+1]*position_continuum[i+1]+position_continuum[i+2]*position_continuum[i+2];
	}
	for(int i=0; i<n_steps+1; i++){
		double mean =0., mean2=0., error=0.;
		double mean_2=0., mean2_2=0., error2=0.;
		for(int walker=0; walker<n_walkers; walker++){
			mean = mean + squared_distance[i+101*walker]/n_walkers;
			mean2 = mean2+squared_distance[i+101*walker]*squared_distance[i+101*walker]/n_walkers;
			mean_2 = mean_2 + squared_distance2[i+101*walker]/n_walkers;
			mean2_2 = mean2_2+squared_distance2[i+101*walker]*squared_distance2[i+101*walker]/n_walkers;
		}
		error = sqrt((mean2-mean*mean)/(n_walkers-1));
		error2 = sqrt((mean2_2-mean_2*mean_2)/(n_walkers-1));
		if(i==0 || i==1){
			error = 0.;
			error2 = 0.;
		}
		else{
			error = error/(2*sqrt(mean));
			error2 = error2/(2*sqrt(mean_2));
		}
		output5 << i<< setw(15) << sqrt(mean) << setw(15) << sqrt(mean_2) << setw(15) << error<< setw(15) << error2 <<endl;
	}
	
	output.close();
	output2.close();
	output3.close();
	output4.close();
	output5.close();
return 0;
}

//implementations of the functions
double integrand(double x){
	return 0.5*M_PI*cos(0.5*M_PI*x);
}
double distribution(double x){
	return 1.-sqrt(1.-x);
}
double sd(vector <double> mean2, vector<double> mean, int n) {
	if (n==1) return 0;
	else 
		return sqrt((mean2[n-1]-mean[n-1]*mean[n-1])/(n-1));
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Studente: Mattia Corigliano, Matr. 944964
_/    _/  _/_/_/  _/_/_/_/ email: mattia.corigliano@studenti.unimi.it
*****************************************************************
*****************************************************************/
