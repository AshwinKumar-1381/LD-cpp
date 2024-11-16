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
#include<omp.h>
#include<chrono>

using namespace std;

namespace program{

    class sysInput{

    	public:
    		
        int nr;
        float t_run1;
        float lj_ep;
        float L, pfrac;
        float m_str, PR, PeA, PeB;
        int N;
        float bias;
        
    };
}

#endif /* LIBRARY_H */