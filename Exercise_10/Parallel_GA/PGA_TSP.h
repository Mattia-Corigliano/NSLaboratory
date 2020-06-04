#include "random.h"
#include <vector>

using namespace std;

void Add_Elite();
void Input();
void Init_prng();
void Export_CitiesConfig();
void Export_OptimalTour();
void Export_Tour();
void Export_CostFunction(int);
double Eval_Distance_2D(double, double, double, double);
void Generate_Tours();
double Eval_CostFunction(vector<int>);
double Eval_TourLenght(vector<int>);
void Eval_Fitness();
void Check();
void Get_OptimalSolution();
void Select_Parents();
void Reproduce();
void CrossOver();
void Migrate();
void Mutate();
void Show_Tours();
void Get_Elite();
void Export_OptimalLenght(int);
void Sort_by_DecrFitness();

Random rnd;

const int pi = 3.1415926535; 

int size, core, migration_step; //variable for parallel computing
int N, M; //number of cities and paths
int gen; //number of generations
int Nelite; //number of elite
double pc; //crossover probability
double pm1, pm2, pm3, pm4, pm5; //mutation probability
int z; //if 0 cities are placed on a circumference, if 1 inside a square
double L; //radius of the circumference, semi-lenght of the square

vector <int> city; //label for the cities
vector <double> city_position; //x-y position of the cities
vector <vector<double>> city_distance; //matrix containing distance between pair of cities
vector <vector<int>> tour; //matrix containing the M Tours around the N cities
vector <double> tour_fitness; //vector containing the fitness of each of the M tours
vector <vector<int>>  parents; //selected parents for reproduction
vector <vector<int>> elite_tour; //matrix containing the tours with the best fitness
