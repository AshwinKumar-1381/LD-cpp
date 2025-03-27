# Makefile for the in-house Langevin Dynamics C++ code

msd: msd_test.cpp fileIO.cpp particles.cpp run.cpp compute.cpp kmc.cpp random.cpp interactions.cpp system.cpp
	g++ -w msd_test.cpp fileIO.cpp particles.cpp run.cpp compute.cpp kmc.cpp random.cpp interactions.cpp system.cpp -o msd_test 
	./msd_test

lanes:

biphasic: