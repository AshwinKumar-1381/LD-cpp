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
using namespace std::chrono;

namespace program{

    class sysInput{

    	public:
    		
        int nr;
        float t_run1;
        float lj_ep;
        float L, S, pfrac;
        float m_str, PR, PeA, PeB;
        int N;
        
    };

    void printElapsed(auto begin, auto end)
    {
        auto dur = duration_cast<seconds>(end - begin);

        int t_time = dur.count();
        int time[3] = {0, 0, 0}; 
        int fac1 = 3600;
        for (int i = 0; i<3; i++)
        {
            time[i] = int(t_time/fac1);
            t_time -= time[i]*fac1;
            fac1 /= 60;    
        }
    
        printf("\nTotal simulation time = %d hrs %d mins %d secs\n", time[0], time[1], time[2]);
    }
}

#endif /* LIBRARY_H */