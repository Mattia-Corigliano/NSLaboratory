CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

ex10.exe : PGA_TSP.o random.o
	$(CC) random.o PGA_TSP.o -o ex10.exe
PGA_TSP.o : PGA_TSP.cpp PGA_TSP.h random.cpp random.h
	$(CC) -c PGA_TSP.cpp -o PGA_TSP.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
