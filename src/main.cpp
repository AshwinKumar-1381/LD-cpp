// main.cpp

#include "library.h"
#include "system.h"

#define INIT_CONFIG_FILE "../init.xyz"

using namespace program;

namespace program{}

void simulation(int nr, float t_simu, float lj_ep, float PeA, float PeB, float m_str, float pfrac, float L, int N);

int main(int argc, char *argv[])
{
    int nr = atoi(argv[1]);
    float t_simu = atof(argv[2]);
    float lj_ep = atof(argv[3]);
    float PeA = atof(argv[4]);
    float PeB = atof(argv[5]);
    float m_str = atof(argv[6]);
    float pfrac = atof(argv[7]);
    float L = atof(argv[8]);
    int N = int(pfrac*L*L); 
    
    printf("%d\n",N);
    
    simulation(nr, t_simu, lj_ep, PeA, PeB, m_str, pfrac, L, N);
    return(0);     
}

void simulation(int nr, float t_simu, float lj_ep, float PeA, float PeB, float m_str, float pfrac, float L, int N)
{
    SimBox *BOX = (SimBox*) malloc(sizeof(SimBox));
    Atoms2D *ATOMS = (Atoms2D*) malloc(N*sizeof(Atoms2D));
    
#ifdef INIT_CONFIG_FILE
    char *init_fname = {INIT_CONFIG_FILE};
    BOX -> initBox(ATOMS, pfrac, L, N, init_fname); 
    
    
    for(int i = 0; i < N; i++) printf("Atom %d - %c %f %f \n", i+1, ATOMS[i].id, ATOMS[i].rx, ATOMS[i].ry);
    printf("Hi\n");
#endif
    
    delete(BOX, ATOMS);
}
