CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex9.exe : GA_TSP.o random.o
	$(CC) random.o GA_TSP.o -o ex9.exe
GA_TSP.o : GA_TSP.cpp GA_TSP.h random.cpp random.h
	$(CC) -c GA_TSP.cpp -o GA_TSP.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
