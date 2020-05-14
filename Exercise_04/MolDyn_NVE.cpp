/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <vector>
#include <iomanip>
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
  Input(); //Inizialization         
  int nconf = 1, measurement_step=10, measurements=nstep/measurement_step, nblocks=20, i=0, j=0, measurements_per_block = measurements/nblocks;
  vector <double> epot(nblocks, 0.), ekin(nblocks, 0.), etot(nblocks, 0.), temp(nblocks, 0.);
  //initialize gdir matrix  
  double gdir[nbins][nblocks];
  for(int i=0; i<nbins; i++){
  	for(int k=0; k<nblocks; k++)
		gdir[i][k] = 0.;
  }
  ofstream output("output.gofr.0");
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%measurement_step == 0){
        Measure();     //Properties measurement
	epot[i]=epot[i] + stima_pot/measurements_per_block;
    	ekin[i]=ekin[i] + stima_kin/measurements_per_block;
        temp[i]=temp[i] + stima_temp/measurements_per_block;
        etot[i]=etot[i] + stima_etot/measurements_per_block;
	for(int bin=0; bin<nbins; bin++)	
		gdir[bin][i]=gdir[bin][i] + g_r[bin]/measurements_per_block;
	j = j+1; //update measurements counter
	if(j==measurements_per_block){
		for(int bin=0; bin<nbins; bin++){
			double r = bin*bin_size;	
			output << i+1 << setw(15) << r << setw(15) << gdir[bin][i] << endl;
		}
		i=i+1; //update block counter
		j=0;
	}
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
  }
  ConfFinal();       //Write last two configurations

  //blocking method
  cout << endl;
  cout << "Blocking Method" << endl << endl;
  ofstream outfile, outfile2, outfile3, outfile4, outfile5;
  outfile.open("ave_epot.out");
  outfile2.open("ave_ekin.out");
  outfile3.open("ave_etot.out");
  outfile4.open("ave_temp.out");
  outfile5.open("ave_gdir.out");
  for(int block = 0; block<nblocks; block++){
	cout << block+1 << endl;
	double ave_epot=0., ave2_epot=0., ave_ekin=0., ave2_ekin=0., ave_etot=0., ave2_etot=0., ave_temp=0., ave2_temp=0.;
	double error_epot=0., error_ekin=0., error_etot =0., error_temp=0.;
	for(int k=0; k<block+1; k++){
		ave_epot = ave_epot + epot[k]/(block+1);
		ave2_epot = ave2_epot + epot[k]*epot[k]/(block+1);
		ave_ekin = ave_ekin + ekin[k]/(block+1);
		ave2_ekin = ave2_ekin + ekin[k]*ekin[k]/(block+1);
		ave_etot = ave_etot + etot[k]/(block+1);
		ave2_etot = ave2_etot + etot[k]*etot[k]/(block+1);
		ave_temp = ave_temp + temp[k]/(block+1);
		ave2_temp = ave2_temp + temp[k]*temp[k]/(block+1);
	}
	if(block+1==1){
		error_epot = 0.;
		error_ekin = 0.;
		error_etot = 0.;
		error_temp= 0.;
	}
	else{
		error_epot = sqrt((ave2_epot-ave_epot*ave_epot)/block);
		error_ekin = sqrt((ave2_ekin-ave_ekin*ave_ekin)/block);
		error_etot = sqrt((ave2_etot-ave_etot*ave_etot)/block);
		error_temp = sqrt((ave2_temp-ave_temp*ave_temp)/block);
	}
	outfile  << setw(15) << (block+1)*measurements_per_block << setw(15) << ave_epot << setw(15) << error_epot << endl;
	outfile2 << setw(15) << (block+1)*measurements_per_block << setw(15) << ave_ekin << setw(15) << error_ekin << endl;
	outfile3 << setw(15) << (block+1)*measurements_per_block << setw(15) << ave_etot << setw(15) << error_etot << endl;
	outfile4 << setw(15) << (block+1)*measurements_per_block << setw(15) << ave_temp << setw(15) << error_temp << endl;
	for(int bin=0; bin < nbins; bin++){
		double r = bin*bin_size;	
		double ave_gdir=0., ave2_gdir=0.;
        	double error_gdir=0.;
		for(int k=0; k<block+1; k++){
			ave_gdir = ave_gdir + gdir[bin][k]/(block+1);
			ave2_gdir = ave2_gdir + gdir[bin][k]*gdir[bin][k]/(block+1);
  		}
		if(block+1==1)
			error_gdir = 0.;
		else
			error_gdir = sqrt((ave2_gdir-ave_gdir*ave_gdir)/block);

		outfile5 << setw(15) << (block+1)<< setw(15) << r << setw(15) << ave_gdir << setw(15) << error_gdir << endl;
	}
  }
  
  outfile.close();
  outfile2.close();
  outfile3.close();
  outfile4.close();
  outfile5.close();
  return 0;
}

//*******************************************************************************************************************************
void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf, ReadConf2;
  double ep, ek, pr, et, vir;
  int bol;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input
  
  ReadInput >> bol;
  if(bol==1) {
  	improved_mode = true;
	cout << "The program runs in the improved mode: True" << endl << endl;
  }
  else {
  	improved_mode = false;
	cout << "The program runs in the improved mode: False" << endl << endl;
  }

  ReadInput >> temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;
  bin_size = (box/2.0)/(double)nbins;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

/*
//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables
*/

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl;
  if(improved_mode==true) cout << "Read initial configuration from file old.0" << endl;
  ReadConf.open("config.0");
  ReadConf2.open("old.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
    if(improved_mode==true){
    	ReadConf2 >> xold[i] >> yold[i] >> zold[i];
    	xold[i] = xold[i]*box;
        yold[i] = yold[i]*box;
        zold[i] = zold[i]*box;
    }
  }
  ReadConf.close();
  ReadConf2.close();

//Prepare initial velocities and old configurations
   if(improved_mode==false){
   	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   	double sumv[3] = {0.0, 0.0, 0.0};
   	for (int i=0; i<npart; ++i){
     		vx[i] = rand()/double(RAND_MAX) - 0.5;
     		vy[i] = rand()/double(RAND_MAX) - 0.5;
     		vz[i] = rand()/double(RAND_MAX) - 0.5;

     		sumv[0] += vx[i];
     		sumv[1] += vy[i];
     		sumv[2] += vz[i];
   	}
   	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   	double sumv2 = 0.0, fs;
   	for (int i=0; i<npart; ++i){
     		vx[i] = vx[i] - sumv[0];
     		vy[i] = vy[i] - sumv[1];
     		vz[i] = vz[i] - sumv[2];

     		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  	}
   	sumv2 /= (double)npart;

   	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   	for (int i=0; i<npart; ++i){
    		vx[i] *= fs;
     		vy[i] *= fs;
     		vz[i] *= fs;

     		xold[i] = Pbc(x[i] - vx[i] * delta);
     		yold[i] = Pbc(y[i] - vy[i] * delta);
     		zold[i] = Pbc(z[i] - vz[i] * delta);
   	}
   }

//from x and xold compute x(t+dt) and correct velocity to obtain x_new(t) = x_old and x(t+dt) = x
   if(improved_mode==true){
   	//compute r(t+dt) with one step of the Verlet algorithm
	double v2 = 0., fs, xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];
	for(int i=0; i<npart; ++i){ //Force acting on particle i
    		fx[i] = Force(i,0);
    		fy[i] = Force(i,1);
    		fz[i] = Force(i,2);
  	}
	for(int i=0; i<npart; i++){
    		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );
		vx[i] = Pbc(xnew - x[i])/(2.0 * delta);
    		vy[i] = Pbc(ynew - y[i])/(2.0 * delta);
    		vz[i] = Pbc(znew - z[i])/(2.0 * delta);
		v2 = v2 + pow(vx[i], 2.) + pow(vy[i], 2.) + pow(vz[i], 2.);
		x[i]=xnew;
		y[i]=ynew;
		z[i]=znew;
	}
	v2 = v2/npart;
	fs = sqrt(3*temp/v2);
	for(int i=0; i<npart; i++){
		vx[i] = vx[i]*fs;
    		vy[i] = vy[i]*fs;
    		vz[i] = vz[i]*fs;
	
		xold[i] = Pbc(x[i] - vx[i] * delta);
     		yold[i] = Pbc(y[i] - vy[i] * delta);
     		zold[i] = Pbc(z[i] - vz[i] * delta);
	}
   }
   return;
}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme
    //compute r(t+dt) with one step of the Verlet algorithm
    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }

  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  double v, t, vij;
  double dx, dy, dz, dr;
 
  ofstream Epot, Ekin, Etot, Temp;
  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);
  //Gdir.open("output_gdir.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  for (int k=0; k<nbins; ++k) g_r[k]=0.0;

  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)
     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     //Potential energy
     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       v += vij;
     }
     
     //g(r) radial distribution function
     //update of the histogram of g(r)
     for(int k=0;k<nbins;k++){
        double r = k*bin_size;
    	double vol = 4.*M_PI*(pow(r+bin_size, 3.)-pow(r, 3.))/3.;
     	if(dr>=k*bin_size && dr< (k+1)*bin_size)
        	g_r[k]=g_r[k] + 2./(vol*rho*npart);
        //Gdir << r << setw(15) << g_r[k] << endl;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle

  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Temp << stima_temp << endl;
  Etot << stima_etot << endl;
  
  Epot.close();
  Ekin.close();
  Temp.close();
  Etot.close();
  //Gdir.close();

  return;
}

void ConfFinal(void){ //Write the last two configurations
  ofstream WriteConf, WriteConf2;

  cout << "Print last two configurations to file config.second_last config.last" << endl << endl;
  WriteConf.open("config.second_last");
  WriteConf2.open("config.last");

  for (int i=0; i<npart; ++i){
    WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
    WriteConf2 << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  WriteConf2.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
