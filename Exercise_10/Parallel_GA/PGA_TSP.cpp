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
#include <random>
#include "PGA_TSP.h"
#include "random.h"
#include "mpi.h" //for parallel computing

using namespace std;

int main (int argc, char* argv[]){
		
	//*******************PARALLEL*************************	
	MPI_Init(&argc, &argv);

	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &core);

	Init_prng(); //initializing parallel random number generator
	
	Input(); //initialing variables for all cores

	//Export_CitiesConfig();
	/*
	if(core == 2){
		cout << "VISUALIZING CORE " <<core << endl << endl;
		cout << "Generation "<<0<< endl;
	}
	*/
	Generate_Tours();
	Check();
	Sort_by_DecrFitness();
	/*
	if(core == 2)
		Get_OptimalSolution();
	*/
	Export_OptimalLenght(0);
	//Export_Tour();

	//cycle over generations
	for(int k=1; k<=gen; k++){
		/*
		if(core == 2)
			cout << endl << "Generation "<<k<< ":" << endl;
		*/
		Get_Elite();		
		Select_Parents();
		CrossOver();
		Reproduce();
		Mutate();
		Add_Elite();
		Check();
		Sort_by_DecrFitness();
		
		//performing migration of the best solutions between cores
		if(k%migration_step == 0)
			Migrate();
		/*
		if(core == 2)
			Get_OptimalSolution();
		*/
		Export_OptimalLenght(k);
		
		//Export_CostFunction(k);
		if(k==gen){
			Export_OptimalTour();
			for(int k=0; k<4; k++){
				if(core==k)
					cout << Eval_TourLenght(tour[0]) << endl;
			}
		}
		/*
		if(core == 2)
			cout << endl;
		*/
	}

	MPI_Finalize();

	//********************************************************
return 0;
}

//******************************************************************************************************
//functions implementations
//******************************************************************************************************

void Add_Elite(){
	//replace worst with selected elite
	Sort_by_DecrFitness();
	for(int k=0; k<Nelite; k++)
		tour[M-1-k] = elite_tour[k];
}

void Check(){
	//cout <<endl<<"assuring the goodness of the generated tours..." << endl;
	bool tours_are_okay = true;
	double s;
	double control = 0.5*N*(N+1); //if all the cities are in the tour once this is the sum of them
	for(int k=0; k<M; k++){
		s = 0.;
		if(tour[k][0] != 1)
			tours_are_okay = false;
		for(int j = 0; j<N; j++)
			s = s + tour[k][j];
		if(s != control)
			tours_are_okay = false;
	}
	if(tours_are_okay == false){
		cout << "ERROR! : tours are not okay" << endl;
		exit(EXIT_FAILURE);
	}		
}

void CrossOver(){
	if(rnd.Rannyu() < pc){
		double r;
		int cut, dim, order_found;
		for(int k=0; k<M/2; k++){
			//randomly generating the location of the cut
			r = rnd.Rannyu(N-4, N-2);
			cut = (int)r;
			dim = N-cut-1;
			order_found = 0.;
			vector <int> cutted1(dim, 0);
			vector <int> cutted2(dim, 0);
			vector <int> ordered(2*dim, 0);
			
			//generating the cutted parts
			for(int j=0; j<dim; j++){
				cutted1[j] = parents[2*k][r+1+j];
				cutted2[j] = parents[2*k+1][r+1+j];
			}
	
			//sorting the cutted parts
			for(int j=0; j<N; j++){
				double temp_city = parents[2*k+1][j];
				for(int i = 0; i<dim; i++){
					if(temp_city == cutted1[i]){
						ordered[order_found] = temp_city;
						order_found ++;
					}
				}
			}
			for(int j=0; j<N; j++){
				double temp_city = parents[2*k][j];
				for(int i = 0; i<dim; i++){
					if(temp_city == cutted2[i]){
						ordered[order_found] = temp_city;
						order_found ++;
					}
				}
			}

			//refilling with sorted orders
			for(int j=0; j<dim; j++){
				parents[2*k][r+1+j] = ordered[j];
			       	parents[2*k+1][r+1+j] =ordered[j+dim];
			}
		}
	}
}

double Eval_Distance_2D(double x1, double y1, double x2, double y2){
       return sqrt(pow(x1-x2, 2.) + pow(y1-y2, 2.));
}

double Eval_CostFunction(vector<int> v){
	double loss = 0.;
	for(int k=0; k<N; k++){
		int first_city, second_city;
		if(k==N-1){
			first_city = v[k];
			second_city = v[0];
		}
		else{
			first_city = v[k];
			second_city = v[k+1];
		}
		loss = loss + city_distance[first_city-1][second_city-1];
	}
	return loss;
}

double Eval_TourLenght(vector<int> v){
	double lenght = 0.;
	for(int k=0; k<N; k++){
		int first_city, second_city;
		if(k==N-1){
			first_city = v[k];
			second_city = v[0];
		}
		else{
			first_city = v[k];
			second_city = v[k+1];
		}
		lenght = lenght + city_distance[first_city-1][second_city-1];
	}
	return lenght;
}

void Eval_Fitness(){
	for(int k=0; k<M; k++)
		tour_fitness[k] = pow(Eval_CostFunction(tour[k]), -1.);
}

void Export_CitiesConfig(){
	ofstream config("cities_configuration.0");
	for(int k=0; k<N; k++)
		config << k+1 << setw(15) << city_position[2*k] << setw(15)<< city_position[2*k+1] << endl;

	config.close();
}

void Export_Tour(){
	ofstream out("tour_convergence.0", ios::app);
	int temp_city, which_city;
	for(int k = 0; k<N; k++){
		temp_city = tour[0][k];
		for(int j=0; j<N; j++){
			if(temp_city == city[j])
				which_city = j;
		}
		out << which_city +1 << setw(15) << city_position[2*which_city] << setw(15) << city_position[2*which_city+1] << endl;
	}
	out << 1 << setw(15) << city_position[0] << setw(15) << city_position[1] << endl;
	out.close();
}

void Export_OptimalTour(){

	int temp_city, which_city;
	for(int k=0; k<N; k++){
		temp_city = tour[0][k];
		for(int j=0; j<N; j++){
			if(temp_city == city[j])
				which_city = j;
		}
		if(core == 0){
			ofstream out1("core1/optimal_tour.0", ios::app);
			out1 << which_city+1 << setw(15) << city_position[2*which_city] << setw(15) << city_position[2*which_city+1] << endl;
			out1.close();
		}
		if(core == 1){
			ofstream out2("core2/optimal_tour.0", ios::app);	
			out2 << which_city+1 << setw(15) << city_position[2*which_city] << setw(15) << city_position[2*which_city+1] << endl;
			out2.close();
		}
		if(core == 2){
			ofstream out3("core3/optimal_tour.0", ios::app);
			out3 << which_city+1 << setw(15) << city_position[2*which_city] << setw(15) << city_position[2*which_city+1] << endl;
			out3.close();
		}
		if(core == 3){
			ofstream out4("core4/optimal_tour.0", ios::app);
			out4 << which_city+1 << setw(15) << city_position[2*which_city] << setw(15) << city_position[2*which_city+1] << endl;
			out4.close();
		}
	}
	//adding the return to the first city

	if(core == 0){
		ofstream out1("core1/optimal_tour.0", ios::app);
		out1 << 1 << setw(15) << city_position[0] << setw(15) << city_position[1] << endl;
		out1.close();
	}	
	if(core == 1){
		ofstream out2("core2/optimal_tour.0", ios::app);
		out2 << 1 << setw(15) << city_position[0] << setw(15) << city_position[1] << endl;
		out2.close();
	}
	if(core == 2){
		ofstream out3("core3/optimal_tour.0", ios::app);
		out3 << 1 << setw(15) << city_position[0] << setw(15) << city_position[1] << endl;
		out3.close();
	}
	if(core == 3){
		ofstream out4("core4/optimal_tour.0", ios::app);
		out4 << 1 << setw(15) << city_position[0] << setw(15) << city_position[1] << endl;	
		out4.close();
	}
}

void Export_OptimalLenght(int gen){
	
	double opt_lenght = Eval_TourLenght(tour[0]);

	if(core == 0){
		ofstream out1("core1/optimal_length.0", ios::app);
		out1 << gen << setw(15) << opt_lenght << endl;
		out1.close();
	}
	if(core == 1){
		ofstream out2("core2/optimal_length.0", ios::app);		
		out2 << gen << setw(15) << opt_lenght << endl;
		out2.close();
	}
	if(core == 2){
		ofstream out3("core3/optimal_length.0", ios::app);				
		out3 << gen << setw(15) << opt_lenght << endl;
		out3.close();
	}
	if(core == 3){
		ofstream out4("core4/optimal_length.0", ios::app);
		out4 << gen << setw(15) << opt_lenght << endl;
		out4.close();
	}
}

void Export_CostFunction(int gen){
	ofstream out("ave_lossfunction.0", ios::app);
	double ave_loss = 0.;
	for(int k=0; k<M/2; k++)	
		ave_loss = ave_loss + Eval_CostFunction(tour[k])/(M/2.);
	out << gen << setw(15) << ave_loss  << endl;
	out.close();
}

void Get_Elite(){
	for(int k=0; k<Nelite; k++)
		elite_tour[k] = tour[k];
}

void Get_OptimalSolution(){
	int opt = 0;
	double opt_lenght= Eval_TourLenght(tour[0]);
	for(int k=1; k<M; k++){
		double temp_lenght = Eval_TourLenght(tour[k]);
		if(temp_lenght < opt_lenght){
			opt = k;
			opt_lenght = temp_lenght;
		}
	}
	cout << "	the optimal solution is the tour: ";
	for (int k=0; k<N; k++)
		cout << tour[opt][k] << " ";
	cout << endl << "	of lenght  = " << opt_lenght << endl;
}

void Generate_Tours(){
	//cout << "Generating initial tours..." << endl;
	int swap_index, swap_value;
	double r;
    	for(int k=0; k<M; k++){
		vector <int> temp_tour = city;
		for(int j=1; j<N; j++){
			//int swap_value = temp_tour[j];
			r = rnd.Rannyu(1, N);
			swap_index = (int)r;
			if(swap_index != j){
				swap_value = temp_tour[j];
				temp_tour[j] = temp_tour[swap_index];
				temp_tour[swap_index] = swap_value;
			}
		}
		tour[k] = temp_tour;
	}
}

void Init_prng(){
	int seed[4];
   	int temp_p1, temp_p2, p1, p2;
   	
	ifstream Primes("Primes");

   	if (Primes.is_open()){
		for(int k=0; k<4; k++){
			Primes >> temp_p1 >> temp_p2;
			if(core==k){
				//ofstream out("primes_core.txt", ios::app);
				p1 = temp_p1;
				p2 = temp_p2;
				//out << core << " " << p1 << " " << p2 << endl;
				//out.close();
			}
		}
	} 
	else cerr << "PROBLEM: Unable to open Primes" << endl;
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

void Input(){
	if(core==0){
		cout<<endl << "---------------------------------------------------------------"<<endl;
        	cout<< "Traveling Salesman Problem with Parallelized Genetic Algorithm"<<endl;
		cout<< "The code performs migration between cores" << endl;
	}
	ifstream input("input.txt");
	input >> N >> M; //read number of cities and paths
	input >> gen; //read number of generations
	input >> migration_step; //read migration step	
	input >> Nelite; //read number of elite
	input >> pc; //read crossover probabilty
	input >> pm1 >> pm2 >> pm3 >> pm4 >> pm5; //read mutation operators probabilities
	input >> z; //if 0 cities on a circumference, else inside a square
	input >> L; //read radius if the cirmuference, semi-lenght of the square
	
	if(core==0){
		cout << endl << "Number of cores running in parallel = " << size << endl;
		cout << "Number of cities = " << N << endl;
		cout << "Number of paths considered = " << M << endl;
		cout << "Number of generations = " << gen << endl;
		cout << "Migration step = " << migration_step <<endl;
		cout << "Number of elite individuals = " << Nelite << endl;
	        cout<< "---------------------------------------------------------------"<<endl;
        	cout<<endl;
	}	
	
	//initialize cities
	//cout << "initializing cities...." << endl;
	for(int k=0; k<N; k++)
		city.push_back(k+1);

	//initialize cities positions
	//cout << "initializing cities positions..." << endl;
	for(int k=0; k<N; k++){
		city_position.push_back(0.);
		city_position.push_back(0.);
	}
	if(core == 0){
		if(z==0){
	
			//cout << "Placing cities on a circumference of radius r = " << L << endl;
			for (int k=0; k<N; k++){
				double theta = rnd.Rannyu(0., 2*pi);
				double x = L*cos(theta), y = L*sin(theta);
				city_position[2*k] = x;
				city_position[2*k+1] = y;	
			}
		}
		if(z==1){
			//cout << "Placing cities inside a square of edge L = " << 2*L << endl;
			for (int k=0; k<N; k++){
				double x = rnd.Rannyu(-L, L);
				double y = rnd.Rannyu(-L, L);
				city_position[2*k] = x;
				city_position[2*k+1] = y;	
			}
		}
	}
	MPI_Bcast(city_position.data(), 2*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	//computing distances between pair of cities
	//cout  << "computing distances between all pair of cities..." << endl << endl;
	for(int k=0; k<N; k++){
		vector <double> temp_row;
		for(int j=0; j<N; j++)
			temp_row.push_back(Eval_Distance_2D(city_position[2*k],city_position[2*k+1], city_position[2*j], city_position[2*j+1]));					
		city_distance.push_back(temp_row);
	}

	//initialize tours matrix and vector of fitnesses to zero
	for(int k=0; k<M; k++){
		vector <int> temp_row(N, 0);
		tour.push_back(temp_row);
		parents.push_back(temp_row);
		tour_fitness.push_back(0.);
	}

	//initialize elite tour
	for(int k=0; k<Nelite; k++){
		vector<int> temp_row(N, 0);
		elite_tour.push_back(temp_row);
	}
}

void Migrate(){
	vector<vector<int>> migrant_tour;
	for(int k=0; k<4; k++){
		vector<int> temp_migrant(N, 0);
		migrant_tour.push_back(temp_migrant);
	}
	for(int k=0; k<4; k++){
		if(core == k)
			migrant_tour[k] = tour[0];
		
		MPI_Bcast(migrant_tour[k].data(), N, MPI_INTEGER, k, MPI_COMM_WORLD);
	}

	double r = rnd.Rannyu(1., 4.);
	int how_to_swap = (int)r; //integer 1, 2, 3

	if(how_to_swap==1){
		vector<vector<int>> copy_migrant_tour = migrant_tour;
		migrant_tour[0] = copy_migrant_tour[3]; //4-->1
		migrant_tour[1] = copy_migrant_tour[0]; //1-->2
		migrant_tour[2] = copy_migrant_tour[1]; //2-->3
		migrant_tour[3] = copy_migrant_tour[2]; //3-->4
	}
	if(how_to_swap==2){
		vector<vector<int>> copy_migrant_tour = migrant_tour;
		migrant_tour[0] = copy_migrant_tour[2]; //3-->1
		migrant_tour[1] = copy_migrant_tour[3]; //4-->2
		migrant_tour[2] = copy_migrant_tour[0]; //1-->3
		migrant_tour[3] = copy_migrant_tour[1]; //3-->4
	}
	if(how_to_swap==3){
		vector<vector<int>> copy_migrant_tour = migrant_tour;
		migrant_tour[0] = copy_migrant_tour[1]; //2-->1
		migrant_tour[1] = copy_migrant_tour[2]; //3-->2
		migrant_tour[2] = copy_migrant_tour[3]; //4-->3
		migrant_tour[3] = copy_migrant_tour[0]; //1-->4
	}

	for(int k=0; k<4; k++){
		if(core == k)
			tour[0] = migrant_tour[k];
	}
	
	Sort_by_DecrFitness();
}

void Mutate(){
	//cout << "	mutations now take place...." << endl;
	
	double r1, r2;
	int shift, index1, index2, swap_value;
	
	//mutations
    	for(int k=0; k<M; k++){
		double which_mut_operator = rnd.Rannyu();
		
		//mutation at a single site
		if(which_mut_operator < pm1){
			do{
				r1 = rnd.Rannyu(1, N);
				r2 = rnd.Rannyu(1, N);
				index1 = (int)r1;
				index2 = (int)r2;
			}while(index1 == index2);
			swap_value = tour[k][index1];
			tour[k][index1] = tour[k][index2];
			tour[k][index2] = swap_value;
		}
		
		//mutations double site
		else if(which_mut_operator >= pm1 && which_mut_operator  < pm1 + pm2){
			r1 = rnd.Rannyu(1, N-3);
			index1 = (int)r1; //int rnd numbers 1, 2, ..., N-4
			
			swap_value = tour[k][index1];
			tour[k][index1] = tour[k][index1+2];
			tour[k][index1+2] = swap_value;

			swap_value = tour[k][index1+1];
			tour[k][index1+1] = tour[k][index1+3];
			tour[k][index1+3] = swap_value;
		}
		

		//inversion of a group of 3-4  genes
		else if(which_mut_operator >= pm1+pm2 && which_mut_operator < pm1+pm2+pm3){
			r1 = rnd.Rannyu(1, N-4); //selection of the LHS of the group
			index1 = (int)r1; //rnd number 1--> N-5
			r2 = rnd.Rannyu(3, 5);  //random selection of the lenght of the group to invert (3 or 4)
			index2 = (int)r2;
			index2 = index1+index2; 

			vector<int> selected_region;
			for(int j=index1; j<=index2; j++)
				selected_region.push_back(tour[k][j]);
			reverse(selected_region.begin(), selected_region.end());
			for(int j=0; j<(int)selected_region.size(); j++)
				tour[k][index1+j] = selected_region[j];
		}

		//shift of some positions of a group of 3 Genes
		else if(which_mut_operator >= pm1+pm2+pm3 && which_mut_operator < pm1+pm2+pm3+pm4){
			r1 = rnd.Rannyu(1, N-4); 
			index1 = (int)r1; //int from 1 to N-5
			index2 = index1+3; //int from 4 to N-2 
			double r = rnd.Rannyu(1, N-index2);
			shift = (int)r;
		
			for(int j=0; j<(index2-index1); j++){
				swap_value = tour[k][index1+j];
				tour[k][index1+j] = tour[k][index1+j+shift];
				tour[k][index1+j+shift] = swap_value;
			}
		}

		//push at the end first 3/4 cities (except the first city)
		else if(which_mut_operator >= pm1+pm2+pm3+pm4 && which_mut_operator < pm1+pm2+pm3+pm4+pm5){
			r1 = rnd.Rannyu(3, 5);
			index1 = (int)r1; //int rnd 3 or 4
			for(int j=1; j<=index1; j++){
				swap_value = tour[k][j];
				tour[k][j] = tour[k][N-index1-1+j];
				tour[k][N-index1-1+j] = swap_value;
			}
		}

	}
}

void Reproduce(){
	for(int k=0; k<M; k++)
		tour[k] = parents[k];
}

void Select_Parents(){
	
	//roulette wheel selection
	double tot_fitness=0.;
	for(int k=0; k<M; k++)
		tot_fitness = tot_fitness + tour_fitness[k];
	
	for(int k=0; k<M; k++){
		double partial_fitness = 0.;
		double r = rnd.Rannyu(0, tot_fitness); //random number between zero and the total of the fitnesses
		for(int j=0; j<M; j++){
			partial_fitness = partial_fitness + tour_fitness[j];
			if(partial_fitness > r){
				parents[k] = tour[j];
				break;
			}
		}
	}
}

void Show_Tours(){
	cout << "Tours: " << endl;
	for(int k=0; k<M; k++){
		for(int j=0; j<N; j++)
			cout << tour[k][j] << " ";
		cout << " " << Eval_TourLenght(tour[k]) << "   " << tour_fitness[k] << "     " << Eval_CostFunction(tour[k]) <<  endl;
	}
}

void Sort_by_DecrFitness(){
	
	Eval_Fitness();

	vector<int> swap_tour(N, 0);
	double swap_fitness;
	int index_max;
	double max;
	for(int k=0; k<M; k++){
		index_max = k;
		max = tour_fitness[k];
		//find the maximum value in k-->M-1
		for(int j=k+1; j<M; j++){
			if(tour_fitness[j] > max){
				index_max = j;
				max = tour_fitness[j];
			}
		}
		//swap k with index_max
		swap_fitness = tour_fitness[k];
		tour_fitness[k] = tour_fitness[index_max];
		tour_fitness[index_max] = swap_fitness;
		
		swap_tour = tour[k];
		tour[k] = tour[index_max];
		tour[index_max] = swap_tour;
	}
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
