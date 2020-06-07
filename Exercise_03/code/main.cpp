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
//***********************************ESERCIZIO 3.1*********************************************************
// PLAIN VANILLA OPTION PRICING
// parameters
	int M = 1E6; //number of asset prices
	int N = 1E4; //number of blocks
	int n = M/N;
	int S_i = 100;
	double T = 1.; //delivery time
	int time_steps = 100;
	double dt = T/time_steps;
	double strike_price = 100.;
	double r = 0.1; //risk-free interest rate
	double volatility = 0.25;
	
	vector<double> call(N, 0.);
	vector<double> call2(N, 0.);
	vector<double> put(N, 0.);
	vector<double> put2(N, 0.);

	for(int i=0; i<N; i++){
		double S1, S2;
		for(int j=0; j<n; j++){ 
			S1 = S_i * exp((r-0.5*volatility*volatility)*T+volatility*rnd.Gauss(0., 1.));
			call[i] = call[i] + exp(-r*T)*max(0., S1-strike_price)/n;
			put[i] = put[i] + exp(-r*T)*max(0., strike_price-S1)/n;
			S2 = S_i;
			for(int k=0; k<time_steps; k++)
				S2 = S2 * exp((r-0.5*volatility*volatility)*dt+volatility*rnd.Gauss(0., 1.)*sqrt(dt));
			call2[i] = call2[i] + exp(-r*T)*max(0., S2-strike_price)/n;
			put2[i] = put2[i] + exp(-r*T)*max(0., strike_price-S2)/n;	
		}
	}
	for(int i=0; i<N; i++){
		double ave_call = 0., ave2_call=0., error_call = 0., ave_put=0., ave2_put=0., error_put=0.;
		double Ave_call = 0., Ave2_call=0., Error_call = 0., Ave_put=0., Ave2_put=0., Error_put=0.;
		for(int j=0; j<i+1; j++){
			ave_call = ave_call + call[j];
			ave2_call = ave2_call + call[j]*call[j];
			ave_put = ave_put + put[j];
			ave2_put = ave2_put + put[j]*put[j];

			Ave_call = Ave_call + call2[j];
			Ave2_call = Ave2_call + call2[j]*call2[j];
			Ave_put = Ave_put + put2[j];
			Ave2_put = Ave2_put + put2[j]*put2[j];
		}
		ave_call = ave_call/(i+1);
		ave2_call = ave2_call/(i+1);
		ave_put = ave_put/(i+1);
		ave2_put = ave2_put/(i+1);

		Ave_call = Ave_call/(i+1);
		Ave2_call = Ave2_call/(i+1);
		Ave_put = Ave_put/(i+1);
		Ave2_put = Ave2_put/(i+1);
		if(i+1==1){
			error_call = 0.;
			error_put = 0.;
			Error_call = 0.;
			Error_put = 0.;
		}
		else{
			error_call = sqrt((ave2_call-ave_call*ave_call)/(i+1));
			error_put = sqrt((ave2_put-ave_put*ave_put)/(i+1));
			Error_call = sqrt((Ave2_call-Ave_call*Ave_call)/(i+1));
			Error_put = sqrt((Ave2_put-Ave_put*Ave_put)/(i+1));
		}
		output << setw(15) << (i+1)*n << setw(15) << ave_call << setw(15) << error_call << setw(15) << ave_put << setw(15) << error_put << setw(15);
		output <<  Ave_call << setw(15) << Error_call << setw(15) << Ave_put << setw(15) << Error_put << setw(15) << endl;
	}
		
	output.close();
return 0;
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
