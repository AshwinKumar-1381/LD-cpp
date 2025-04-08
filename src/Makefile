# Makefile for the in-house Langevin Dynamics C++ code

CC = g++
CFLAGS = -w -std=c++17

msd_test: msd_test.cpp fileIO.cpp particles.cpp run.cpp compute.cpp kmc.cpp random.cpp interactions.cpp system.cpp
	$(CC) $(CFLAGS) msd_test.cpp fileIO.cpp particles.cpp run.cpp compute.cpp kmc.cpp random.cpp interactions.cpp system.cpp -o msd_test.o 

biphasic: main.cpp fileIO.cpp particles.cpp run.cpp compute.cpp kmc.cpp random.cpp interactions.cpp system.cpp
	$(CC) $(CFLAGS) main.cpp fileIO.cpp particles.cpp run.cpp compute.cpp kmc.cpp random.cpp interactions.cpp system.cpp -o biphasic.o

lane: lane.cpp fileIO.cpp particles.cpp run.cpp compute.cpp kmc.cpp random.cpp interactions.cpp system.cpp
	$(CC) $(CFLAGS) lane.cpp fileIO.cpp particles.cpp run.cpp compute.cpp kmc.cpp random.cpp interactions.cpp system.cpp -o lane.o

clean:
	rm *.o 