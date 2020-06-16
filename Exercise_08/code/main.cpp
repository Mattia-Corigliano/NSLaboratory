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

void Init_prng(Random&);
void Simulated_annealing(double, double, double&, double&, double&,vector <double> &, Random);
void Optim_delta(double& , bool , int , double, double, Random);
double Psi_trial (double, double, double);
bool Accept_reject(double, double, double, double, Random);
double E_loc(double , double , double );
double Potential(double );

int main (){
	//Description
	cout<<endl << "---------------------------------------------------------------"<<endl;
	cout<< "Variational Monte Carlo (VMC) method to calculate ground state energy"<<endl;
	cout<< "---------------------------------------------------------------"<<endl;
	cout<<endl;
	//*****initializing random numbers generation
	Random rnd;
	Init_prng(rnd);
	//initializing output files
	//ofstream output1("sampled_x.txt");
	ofstream output2("Egs.txt");

	int N = 100; //number of blocks
	int n = 1E4;//number of steps per blocks
	double mu = 1., sigma = 0.5; //initial guess of the trial GSwavefunction	
	//trial GS-wavefunction optimization
	double delta;
	int opt_steps = 1000;
	double temp=10; //temperature of the annealing
	
	vector <double> xi(n*N, 0); //1d positions

	Simulated_annealing(opt_steps,temp, mu, sigma, delta, xi, rnd);
	
	vector <double> integral(N, 0.);
	for(int block=0; block<N; block++){
		for(int k=0; k<n; k++)
			integral[block] = integral[block] + E_loc(xi[k+block*n], mu, sigma)/n;
	}
	//data blocking
	for(int i=0; i<N; i++){
		double ave=0., ave2=0., error=0.;
		for(int j=0; j<i+1; j++){
			ave = ave + integral[j]/(i+1);
			ave2 = ave2 + integral[j]*integral[j]/(i+1);
		}
		if(i+1==1)
			error = 0.;
		else
			error = sqrt((ave2-ave*ave)/i);
		output2 << setw(15) << i+1 << setw(15) << ave << setw(15) << error << endl;
	}	
	//output1.close();
	output2.close();
return 0;
}


//******************************************************************************************************
void Init_prng(Random& rnd){
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
}
void Optim_delta(double& delta, bool optimize, int M, double mu, double sigma, Random rnd){
	if(optimize==true){
		delta = 0.;
		double delta_step = 0.02;
		double precision = 0.05, error = 0., f;
		cout << endl;
		cout << "	Performing delta optimization with desired precision = " << precision<<"..." << endl;
		do{
			delta = delta+delta_step;
			f = 0.;
			//setting the initial position in 1d space
			vector <double> r(M, 0.);

			for(int i=1; i<M; i++){
				double x = r[(i-1)] + rnd.Rannyu(-delta, delta);
				if(Accept_reject(x, r[(i-1)], mu, sigma, rnd) == true){
					f++;
					r[i] = x;
				}
				else
					r[i] = r[i-1];
			}
			error = abs(f/M-0.5);
		}while(error>precision);
		cout << "	...optimal delta value = " << delta << endl << endl;
	}
}
void Simulated_annealing(double opt_steps, double temp, double& mu, double& sigma, double& delta, vector <double> & v, Random rnd){
	cout << endl << "Performing simulated annealing optimization of the GSwavefunction parameters..." << endl << endl;
	cout << "Number of iterations =" << opt_steps << endl << endl;
	double guess_mu,guess_sigma, opt_mu, opt_sigma, opt_gs;
	opt_mu = mu;
	opt_sigma=sigma;
	opt_gs=0;
	int sample=1E4;
	double hmu=0.05;
	double hsigma=0.03;

	cout << "Initial guess: mu = " << mu << " sigma = "<< sigma << endl;

	for(int step=1; step<=opt_steps;step++){
		//adaptive learning rate
		if(step>100){
			hmu=0.02;
			hsigma=0.02;
		}
		if(step>500){
			hmu=0.01;
			hsigma=0.01;
		}
		if(step>750){
			hmu=0.005;
			hsigma=0.005;
		}
		//genero nuovi guess
		if(step==1){
			guess_mu = mu;
			guess_sigma = sigma;
		}
		else{
			guess_mu=opt_mu+rnd.Gauss(0., hmu);
			guess_sigma=opt_sigma+rnd.Gauss(0., hsigma);
		}	
		Optim_delta(delta, true, sample, guess_mu, guess_sigma, rnd);
		vector <double> r(sample, 0);
		for (int i=0; i<sample; i++){
			double xnew = r[(i-1)] + rnd.Rannyu(-delta, delta);
			if(Accept_reject(xnew, r[(i-1)], guess_mu, guess_sigma, rnd) == true)
				r[i] = xnew;
			else
				r[i] = r[i-1];
		}
		//con gli xi samplati da psi_trial calcolo E_Gs trial
		double guess_gs = 0.;
		for(int k=0; k<sample; k++)
			guess_gs = guess_gs + E_loc(r[k], guess_mu, guess_sigma)/sample;
		cout << step << " iteration tempted: " << guess_mu << " " << guess_sigma << " " << guess_gs << endl;
		if(guess_gs<opt_gs){
			opt_mu=guess_mu;
			opt_sigma=guess_sigma;
			opt_gs = guess_gs;
		}
		else{
			double r=rnd.Rannyu();
			double th=exp(-(double)step/temp);
			if(r<th){
				opt_mu=guess_mu;
				opt_sigma=guess_sigma;
				opt_gs=guess_gs;
			}
		}
		cout << step << " iteration: opt. mu = "<<opt_mu << " opt. sigma = " << opt_sigma << " opt.egs = " <<opt_gs << endl << endl; 
	}

	ofstream output1("sampled_x.txt");
	mu = opt_mu;
	sigma = opt_sigma;
	cout << "optimization ended with optimal values: mu = "<<mu << " , sigma = " << sigma << endl;
	int size = v.size();
	double f=0.;
	Optim_delta(delta, true, size, mu, sigma, rnd);
	for(int i =0; i<size; i++){
		double xnew = v[(i-1)] + rnd.Rannyu(-delta, delta);
		if(Accept_reject(xnew, v[(i-1)], mu, sigma, rnd) == true){
			v[i] = xnew;
			f++;
		}
		else
			v[i] = v[i-1];

		output1<<setw(15) << i+1<< setw(15) << v[i] << endl;		
	}

	cout << endl <<endl<< "Acceptance rate: " << f/size << endl;
	output1.close();
}

double Psi_trial (double x, double mu, double sigma){
	return exp(-1.* pow(x-mu, 2.)/pow(sigma, 2.)/2.) + exp(-1.* pow(x+mu, 2.)/pow(sigma, 2.)/2.);
}
double Potential(double x){
	return pow(x,4)-2.5*pow(x,2);
}
double E_loc(double x, double mu, double sigma){
	double kinetic = -0.5*(exp(-pow(x+mu,2)/(2*pow(sigma,2)))*((1+exp(2*mu*x/pow(sigma,2)))*(pow(mu,2)-pow(sigma,2)+pow(x,2))-2*mu*x*(-1+exp(2*mu*x/pow(sigma,2))))/pow(sigma,4));
	double e = kinetic/Psi_trial(x,mu,sigma) + Potential(x);
	return e; 
}
bool Accept_reject(double x1, double x2,  double mu, double sigma, Random rnd){
	double a = min(1., pow(Psi_trial(x1, mu, sigma), 2.)/pow(Psi_trial(x2, mu, sigma), 2.));
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
