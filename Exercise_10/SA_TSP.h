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
void Generate_Tour();
double Eval_CostFunction(vector<int>);
double Eval_TourLenght(vector<int>);
void Check();
void CrossOver();
void Metropolis_Step(double);
void Mutate();
void Save_Tour();
void Show_Tours();
void Get_Annealing_Schedule();
void Export_OptimalLenght(int);

Random rnd;

const int pi = 3.1415926535; 

int N, M; //number of cities and paths
int gen; //number of generations
double  T0; //initial temperature for the annealing
double pc; //crossover probability
double pm1, pm2, pm3, pm4, pm5; //mutation probability
int z; //if 0 cities are placed on a circumference, if 1 inside a square
double L; //radius of the circumference, semi-lenght of the square

vector <int> city; //label for the cities
vector <double> temperature; //temperatures for the annealing
vector <double> city_position; //x-y position of the cities
vector <vector<double>> city_distance; //matrix containing distance between pair of cities
vector <vector<int>> tour; //matrix containing the M Tours around the N cities
vector <vector<int>> old_tour;
