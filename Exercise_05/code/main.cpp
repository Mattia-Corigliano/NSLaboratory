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

double p1 (double, double, double);
double p2 (double, double, double);
bool accept_reject1(double, double, double, double, double, double, Random);
bool accept_reject2(double, double, double, double, double, double, Random);
 
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
	ofstream output("output_uniform.txt");
	ofstream output2("output2_uniform.txt");
	ofstream output3("output3_uniform.txt");
	ofstream output4("output4_uniform.txt");
	//ofstream output("output_gaussian.txt");
	//ofstream output2("output2_gaussian.txt");
	//ofstream output3("output3_gaussian.txt");
	//ofstream output4("output4_gaussian.txt");
//***********************************ESERCIZIO 5.1*********************************************************
// parameters
	int M = 1E6; //number of steps
	int N = 1E3; //number of blocks
	int n = M/N;//number of steps per blocks
	//lenght are in unit of the bohr radius
	double delta = 0.; //transition semi-amplituted (x --> x-delta, x+delta etc...)
	// small values of delta results in a large number of transitions accepted.
	// by increasing delta from a low value we can find its optimal value, which results in the 50% of the total transitions accepted
	double precision = 0.001;
	double error = 0.;
	vector <double> r1s(3*M, 0); //3d positions
	vector <double> r2p(3*M, 0);

	//finding the optimal delta value
	bool transition_accepted = false;
	cout << endl;
	cout << "start finding the optimal value for delta..." << endl;
	cout << "desired precision = " <<precision << endl << endl;
	do{
		delta = delta+0.001;
		double f = 0.;
		//setting the initial position in 3d space
		r1s[0] = 10.; //x_0
		r1s[1] = 10.; //y_0
		r1s[2] = 10.; //z_0
		
		for(int i=1; i<M; i++){
			double x = r1s[3*(i-1)] + rnd.Rannyu(-delta, delta);
			double y = r1s[3*(i-1)+1] + rnd.Rannyu(-delta, delta);
			double z =  r1s[3*(i-1)+2] + rnd.Rannyu(-delta, delta);
			//double x = rnd.Gauss(r1s[3*(i-1)], delta);
			//double y = rnd.Gauss(r1s[3*(i-1)+1], delta);
			//double z = rnd.Gauss(r1s[3*(i-1)+2], delta);

			transition_accepted = accept_reject1(x, y, z, r1s[3*(i-1)], r1s[3*(i-1)+1], r1s[3*(i-1)+2], rnd);
			if(transition_accepted == true){
				f++;
				r1s[3*i] = x;
				r1s[3*i+1] = y;
				r1s[3*i+2] = z;
			}
			else{
				r1s[3*i] = r1s[3*(i-1)];
				r1s[3*i+1] =r1s[3*(i-1)+1];
				r1s[3*i+2] = r1s[3*(i-1)+2];
			}
		}
		error = abs(f/M-0.5);
	}while(error>precision);
	cout <<"the search ended with optimal value delta = " << delta << endl;

	//finding the optimal delta value
	transition_accepted = false;
	cout << endl;
	cout << "start finding the optimal value for delta..." << endl;
	cout << "desired precision = " << precision << endl << endl;
	delta = 0.;
	error = 0.;
	do{
		delta = delta+0.005;
		double f = 0.;
		//setting the initial position in 3d space
		r2p[0] = 10.; //x_0
		r2p[1] = 10.; //y_0
		r2p[2] = 10.; //z_0
		
		for(int i=1; i<M; i++){
			double x = r2p[3*(i-1)] + rnd.Rannyu(-delta, delta);
			double y = r2p[3*(i-1)+1] + rnd.Rannyu(-delta, delta);
			double z =  r2p[3*(i-1)+2] + rnd.Rannyu(-delta, delta);
			//double x = rnd.Gauss(r2p[3*(i-1)], delta);
			//double y = rnd.Gauss(r2p[3*(i-1)+1], delta);
			//double z = rnd.Gauss(r2p[3*(i-1)+2], delta);
			transition_accepted = accept_reject2(x, y, z, r2p[3*(i-1)], r2p[3*(i-1)+1], r2p[3*(i-1)+2], rnd);
			if(transition_accepted == true){
				f++;
				r2p[3*i] = x;
				r2p[3*i+1] = y;
				r2p[3*i+2] = z;
			}
			else{
				r2p[3*i] = r2p[3*(i-1)];
				r2p[3*i+1] =r2p[3*(i-1)+1];
				r2p[3*i+2] = r2p[3*(i-1)+2];
			}
		}
		error = abs(f/M-0.5);
	}while(error>precision);
	cout <<"the search ended with optimal value delta = " << delta << endl;
	
	
	for (int i=0; i<M; i++){
		output<<setw(15) << i+1<< setw(15) << r1s[3*i] << setw(15) << r1s[3*i+1] << setw(15) << r1s[3*i+2] << endl;
		output2<<setw(15) << i+1<< setw(15) << r2p[3*i] << setw(15) << r2p[3*i+1] << setw(15) << r2p[3*i+2] << endl;
	}

	//data blocking
	vector <double> radius1s(N, 0.), radius2p(N, 0.);
	for(int i=0; i<N; i++){
		for(int j=0; j<n; j++){
			radius1s[i]=radius1s[i] + sqrt(pow(r1s[3*(j+i*n)], 2.) + pow(r1s[3*(j+i*n)+1], 2.) + pow(r1s[3*(j+i*n)+2], 2.))/n;
			radius2p[i]=radius2p[i] + sqrt(pow(r2p[3*(j+i*n)], 2.) + pow(r2p[3*(j+i*n)+1], 2.) + pow(r2p[3*(j+i*n)+2], 2.))/n;
		}
	}
	for(int i=0; i<N; i++){
		double ave1s=0., ave21s=0., error1s=0., ave2p=0., ave22p=0., error2p=0.;
		for(int j=0; j<i+1; j++){
			ave1s = ave1s + radius1s[j]/(i+1);
			ave21s = ave21s + radius1s[j]*radius1s[j]/(i+1);
			ave2p = ave2p + radius2p[j]/(i+1);
			ave22p = ave22p + radius2p[j]*radius2p[j]/(i+1);
		}
		if(i+1==1){
			error1s = 0.;
			error2p = 0.;
		}
		else{
			error1s = sqrt((ave21s-ave1s*ave1s)/i);
			error2p = sqrt((ave22p-ave2p*ave2p)/i);
		}
		output3 << setw(15) << (i+1)*n << setw(15) << ave1s << setw(15) << error1s << endl;
		output4 << setw(15) << (i+1)*n << setw(15) << ave2p << setw(15) << error2p << endl;
		
	}
	output.close();
	output2.close();
	output3.close();
	output4.close();
return 0;
}

double p1 (double x, double y, double z){
	return exp(-2.*sqrt(x*x+y*y+z*z));
}

double p2 (double x, double y, double z){
	return z*z*exp(-sqrt(x*x+y*y+z*z));
}

bool accept_reject1(double x1, double y1, double z1, double x2, double y2, double z2, Random rnd){
	double a = min(1., p1(x1, y1, z1)/p1(x2, y2, z2));
	bool accept = false;
	if(rnd.Rannyu() <= a) accept=true;
	return accept;
}
bool accept_reject2(double x1, double y1, double z1, double x2, double y2, double z2, Random rnd){
	double a = min(1., p2(x1, y1, z1)/p2(x2, y2, z2));
	bool accept = false;
	if(rnd.Rannyu() <= a) accept=true;
	return accept;
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
