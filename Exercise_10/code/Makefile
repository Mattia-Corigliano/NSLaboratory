CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex10.exe : SA_TSP.o random.o
	$(CC) random.o SA_TSP.o -o ex10.exe
SA_TSP.o : SA_TSP.cpp SA_TSP.h random.cpp random.h
	$(CC) -c SA_TSP.cpp -o SA_TSP.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
