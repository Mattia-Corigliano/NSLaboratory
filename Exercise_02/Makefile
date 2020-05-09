CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex2.exe : main.o random.o
	$(CC) random.o main.o -o ex2.exe
main.o : main.cpp random.cpp random.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
