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

double angle_uniform(Random); // this function returns theta uniformly distributed in [0, pi/2]
double sd(vector<double>, vector<double>, int);
 
int main (){
	//initializing random numbers generation
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
 
	//initializing output files
	ofstream output("output.txt");
	ofstream output2("output2.txt");
	ofstream output3("output3.txt");
	ofstream output4("output4.txt");
	ofstream output5("output5.txt");
//***********************************ESERCIZIO 1.1*********************************************************
	//variables
	int M=1000000;           //Total number of throws
	int N=100;              //Number of blocks
	int L=int(M/N);        //Number of throws in each block, please use for M a multiple of N
	vector<double> rnd_numbers(M, 0.); //vector containing M zeros
	for (int i = 0; i<M; i++)
		rnd_numbers[i]=rnd.Rannyu();
	vector<int> x(N, 0);
	for (int i =0; i<N; i++)
		x[i] = i+1;
	vector<double> ave(N, 0.);
	vector<double> ave2(N, 0.);
	vector<double> sum_prog(N, 0.);
	vector<double> sum2_prog(N, 0.);
	vector<double> err_prog(N, 0.);

	for(int i=0; i<N; i++){
		double sum = 0.;
		for(int j=0; j<L; j++)
			sum = sum + rnd_numbers[j+i*L];
		ave[i] = sum/L;
		ave2[i] = ave[i]*ave[i];
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i] = sum_prog[i] + ave[j]; //SUM_{j=0,i} r_j
        		sum2_prog[i] = sum2_prog[i] + ave2[j]; //SUM_{j=0,i} (r_j)^2
		}
		sum_prog[i]=sum_prog[i]/(i+1); //Cumulative average
    		sum2_prog[i]=sum2_prog[i]/(i+1); //Cumulative square average
    		err_prog[i] = sd(sum2_prog, sum_prog, i); //Statistical uncertainty
		output << double(x[i]*L) << setw(15) << sum_prog[i] << setw(15) << err_prog[i] << endl;
	}

//***********************************ESERCIZIO 1.2*********************************************************
	for(int i=0; i<N; i++) {
		ave[i]=0;
		ave2[i]=0;
		sum_prog[i]=0;
		sum2_prog[i]=0;
		err_prog[i]=0;
	}
	for(int i=0; i<N; i++){
		double sum = 0.;
		for(int j=0; j<L; j++)
			sum = sum + (rnd_numbers[j+i*L]-0.5)*(rnd_numbers[j+i*L]-0.5);
		ave[i] = sum/L;
		ave2[i] = ave[i]*ave[i];
	}
	for(int i=0; i<N; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i] = sum_prog[i] + ave[j]; //SUM_{j=0,i} r_j
        		sum2_prog[i] = sum2_prog[i] + ave2[j]; //SUM_{j=0,i} (r_j)^2
		}
		sum_prog[i]=sum_prog[i]/(i+1); //Cumulative average
    		sum2_prog[i]=sum2_prog[i]/(i+1); //Cumulative square average
    		err_prog[i] = sd(sum2_prog, sum_prog, i); //Statistical uncertainty
		output2 << double(x[i]*L) << setw(15) << sum_prog[i] << setw(15) << err_prog[i] << endl;
	}

//***********************************ESERCIZIO 1.3*********************************************************
	// N --> M = 100 is the number of subintervals in [0,1]
	// L --> n the number of throws 
	// n/M = L/N
	double dx = double(1./N); //width of the ith interval
	double n = double(L/N);
	vector <double> chi(N, 0.);
	double chi2 = 0;

	for(int i=0; i<N; i++){
		for(int j=0; j<L; j++){
			if(rnd_numbers[i*L+j]>=double(i*dx) && rnd_numbers[i*L+j]<double((i+1)*dx))
				chi[i]++;		
		}
		chi[i] = double((chi[i]-n)*(chi[i]-n)/n);
		chi2 = chi2 + chi[i];
		output3 << i+1 << setw(15) << chi[i] << endl;
		
	}
	
//***********************************ESERCIZIO 2*********************************************************
	double realizations = 10000;
	vector <int> points = {1, 2, 10, 100};
	vector <double> Sn_binom(4, 0.);
	vector <double> Sn_exp(4, 0.);
	vector <double> Sn_lor(4, 0.);
	for (int k=0; k<realizations; k++){
		for (int i=0; i<int(points.size()); i++){
			Sn_binom[i] = 0.;
			Sn_exp[i] = 0.;
			Sn_lor[i] = 0.;
			for (int j=0; j<points[i]; j++){
				Sn_binom[i]=Sn_binom[i]+double(int(rnd.Rannyu(0, 6))+1)/points[i];
				Sn_exp[i]=Sn_exp[i]+rnd.Exp(1.)/points[i];
				Sn_lor[i]=Sn_lor[i] + rnd.Lor(0., 1.)/points[i];
			}
		}
		output4 << Sn_binom[0] << setw(15) << Sn_binom[1] << setw(15) << Sn_binom[2] << setw(15) << Sn_binom[3]<<setw(15);
		output4 << Sn_exp[0] << setw(15) << Sn_exp[1] << setw(15) << Sn_exp[2] << setw(15) << Sn_exp[3]<<setw(15);
		output4 << Sn_lor[0] << setw(15) << Sn_lor[1] << setw(15) << Sn_lor[2] << setw(15) << Sn_lor[3] << endl;
	}
//***********************************ESERCIZIO 3*********************************************************
//SIMULATION OF THE BUFFON EXPERIMENT
	double lines_width = 1.;
	double needle_lenght = lines_width/2.;
	int n_needles=1E6;
	int n_lines_crossed;
	int n_blocks = 1E3;
	int throws_per_block = int(n_needles/n_blocks);
	vector <double> pi (n_blocks, 0.);

	for(int j = 0; j<n_blocks; j++){
		n_lines_crossed=0;
		for (int i = 0; i<throws_per_block; i++){
			double x1 = rnd.Rannyu(0., lines_width/2.); //draw uniformly in [0, nd] the bottom of the i-th needle
			double theta = angle_uniform(rnd);
			if(x1<=needle_lenght/2.*sin(theta))
				n_lines_crossed++;
		}
		pi[j] = 2.*needle_lenght/lines_width*(double(throws_per_block)/double(n_lines_crossed));
	}
	x.resize(n_blocks);
	x.shrink_to_fit();
	sum_prog.resize(n_blocks);
	sum_prog.shrink_to_fit();
	sum2_prog.resize(n_blocks);
	sum2_prog.shrink_to_fit();
	err_prog.resize(n_blocks);
	err_prog.shrink_to_fit();
	for(int i=0; i<n_blocks; i++) {
		sum_prog[i]=0;
		sum2_prog[i]=0;
		err_prog[i]=0;
	}
	for (int i =0; i<n_blocks; i++)
		x[i] = i+1;
	for(int i = 0;i<n_blocks;i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i] = sum_prog[i] + pi[j]; //SUM_{j=0,i} r_j
        		sum2_prog[i] = sum2_prog[i] + pi[j]*pi[j]; //SUM_{j=0,i} (r_j)^2
		}
		sum_prog[i]=sum_prog[i]/(i+1); //Cumulative average
    		sum2_prog[i]=sum2_prog[i]/(i+1); //Cumulative square average
    		err_prog[i] = sd(sum2_prog, sum_prog, i); //Statistical uncertainty
		output5 << double(x[i]*throws_per_block) << setw(15) << sum_prog[i] << setw(15) << err_prog[i] << endl;
	}
	
	output.close();
	output2.close();
	output3.close();
	output4.close();
	output5.close();
return 0;
}

//implementations of the functions
double angle_uniform(Random rnd){
	double x, y;
	do{
		x = rnd.Rannyu(0., 1.);
		y = rnd.Rannyu(0., 1.);
		if(sqrt(x*x+y*y)<1.)
			break;
	}while(sqrt(x*x+y*y)>1.);
	double theta=acos(x/(sqrt(x*x+y*y)));
	return theta;
}
		
			
/*double find_min(vector <double> v){
	double min = v[0];
	for(int i=0; i<int(v.size()); i++){
		if(v[i]<min)
			min = v[i];
	}
	return min;
}*/
double sd(vector <double> mean2, vector<double> mean, int n) {
	if (n==0) return 0;
	else 
		return sqrt((mean2[n]-mean[n]*mean[n])/n);
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
