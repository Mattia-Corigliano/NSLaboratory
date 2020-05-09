/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <algorithm>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{ 
  Input(); //Inizialization

  //equilibration of the system
  if(equilibration==0)
    Equilibrate();
    
  //measures on the equilibrated system
  if(equilibration==1){
     for(int iblk=1; iblk <= nblk; ++iblk){ //simulation
       Reset(iblk);   //Reset block averages
       for(int istep=1; istep <= nstep; ++istep){
         Move(metro);
         Measure();
         Accumulate(); //Update block averages
       }
       Averages(iblk);   //Print results for current block
     }
     ConfFinal(); //Write final configuration
  }
 
return 0;
}
//**************************************
void Input(void)
{
  ifstream ReadInput;
  cout << endl << "******************************" << endl;
  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> equilibration; //0 or 1; if 0 the code does an equilibration run, else the code compute averages on the equilibrated system
  ReadInput >> restart; //0 or 1; if 0 the code starts from the T=infty initial configuration, else the code strats from the final configuration reached in the previous run
  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility
 
  n_props = 4; //Number of observables

//initial configuration
  //infinite temperature configuration
  if(restart==0){
  	for(int i=0; i<nspin; ++i){
		if(rnd.Rannyu() >=0.5) s[i]=1;
		else s[i] = -1;
	}
  }
  //equilibrated initial configuration
  if(restart==1){
  	ifstream read("initial_eq_conf.txt");
	for(int i=0; i<nspin; ++i)
		read>>s[i];
	read.close();
  }	
  
//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}
//*********************************************************************
void Equilibrate(void){
  cout << endl;
  cout << "***This is an equilibration run***" << endl << endl;
  ofstream eq("equilibration.txt");
  //write initial energy
  eq << setw(12) << 0. << setw(12) << walker[iu]/(double)nspin << setw(12) << walker[ix]/(double)nspin << endl;
  for(int istep=1; istep <= nstep; ++istep){
      Move(metro);
      Measure();
      if(istep%100==0)
        eq << setw(12) << istep << setw(12) << walker[iu]/(double)nspin << setw(12) <<  walker[ix]/(double)nspin << endl; 
   }
   ConfFinal(); //Write final configuration (it will be the initial configuration for the measures on the equilibrated system)
   eq.close();
}
//***********************************************************************
void Move(int metro)
{
  int o;
  double a, p, dE;

  for(int i=0; i<nspin; ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    attempted++;
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1){ //metropolis
    	dE = Energy_difference(s[o], o);
	a = min(1., exp(-beta*dE));
	if(a==1.){
		s[o] = -s[o];
		accepted++;
	}
	else{
		if(rnd.Rannyu() < a){
			s[o] = -s[o];	
			accepted++;
		}
	}
    }
    else //Gibbs sampling
    {
	accepted++;
	dE = Energy_difference(1, o);
	p = pow(1 + exp(beta*dE), -1.);
	double r = rnd.Rannyu();
	if(r<p)
		s[o]=-1;	
	else
		s[o]=1;
    }
  }
}
//*********************************************************************************
double Energy_difference(int sm, int ip)
{
  double ene = 2.*J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) + 2.*h * sm;
  return ene;
}
//****************************************************************************************
void Measure()
{
  double u = 0.0, m = 0.0;//c=0.0;
//cycle over the actual configuration
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     //c += pow(beta, 2.)*pow(-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]), 2.);
     m += s[i];
  }
  walker[iu] = u;
  walker[ic] = pow(beta, 2.)*pow(u, 2.);
  walker[ix] = beta*pow(m, 2.);
  walker[im] = m;
}
//************************************************************************************
void Reset(int iblk) //Reset block averages
{
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }
   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}
//************************************************************************************
void Accumulate(void) //Update block averages
{
   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}
//*************************************************************************************
void Averages(int iblk) //Print results for current block
{
   ofstream Ene, Heat, Mag, Chi;
   const int wd=12;
    
   cout << "Block number " << iblk << endl;
   cout << "Acceptance rate " << accepted/attempted << endl << endl;
   //Energy
   Ene.open("output.ene",ios::app);
   stima_u = blk_av[iu]/blk_norm/(double)nspin;
   glob_av[iu]  += stima_u;
   glob_av2[iu] += stima_u*stima_u;
   err_u=Error(glob_av[iu],glob_av2[iu],iblk);
   Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
   Ene.close();
   //Heat capacity
   Heat.open("output.heat", ios::app);
   stima_c = (blk_av[ic]/blk_norm-pow(beta, 2.)*pow(stima_u*(double)nspin, 2.))/(double)nspin;
   glob_av[ic] += stima_c;
   glob_av2[ic] += stima_c*stima_c;
   err_c=Error(glob_av[ic], glob_av2[ic], iblk);
   Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
   Heat.close();
   //Magnetic susceptibility
   Chi.open("output.chi", ios::app);
   stima_x = blk_av[ix]/blk_norm/(double)nspin;
   glob_av[ix] += stima_x;
   glob_av2[ix] += stima_x*stima_x;
   err_x = Error(glob_av[ix], glob_av2[ix], iblk);
   Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
   Chi.close();
   //Magnetization
   Mag.open("output.mag", ios::app);
   stima_m = blk_av[im]/blk_norm/(double)nspin;
   glob_av[im] += stima_m;
   glob_av2[im] += stima_m*stima_m;
   err_m = Error(glob_av[im], glob_av2[im], iblk);
   Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
   Mag.close();
   cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    if(iblk==1) return 0.0;
    else{
	double err = (sum2/(double)iblk - pow(sum/(double)iblk,2));
	if(err < 0.)
		err = abs(err);
    	return sqrt(err/(double)(iblk-1));
    }
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
