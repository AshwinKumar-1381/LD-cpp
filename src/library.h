// library.h

#ifndef LIBRARY_H
#define LIBRARY_H

#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string>
#include<fstream>
#include<cmath>
#include<sys/stat.h>
#include<chrono>
#include<vector>
#include<cstring>
#include<any>

using namespace std;
using namespace std::chrono;

namespace program{

    class sysInput{

    	public:
    		
        int nr;
        float dt;
        float L, S, pfrac;
        float m_str, D_str, PR, PeA, PeB;
        int N, nAtomTypes;
        char *folder = new char[50];
        bool restart = false;
        int res_step = 0;
        time_point<high_resolution_clock> begin, end;
        
    };

    void simulation(sysInput *Input);
    char* returnElapsedTime(sysInput *Input);
}

#endif /* LIBRARY_H */