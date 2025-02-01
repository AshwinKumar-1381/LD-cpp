// main.cpp

#include "library.h"
#include "interactions.h"
#include "system.h"
#include "fileIO.h"
//#include "run.h"
//#include "kmc.h"

#define INIT_CONFIG_FILE "../tests/frame0.dat"

using namespace program;

void simulation(sysInput *Input);

int main(int argc, char *argv[])
{
    program::sysInput *Input = new sysInput;

    Input -> nr     = 18;
    Input -> lj_ep  = 1.0;
    Input -> m_str  = 1.0;
    Input -> pfrac  = 0.5;
    Input -> L      = 100.0;
    Input -> S      = 1.0;
    Input -> PR     = 0.0;
    Input -> PeA    = 0.0;
    Input -> PeB    = 2.0;
    
    Input -> N = int(Input->pfrac*Input->L*Input->L*Input->S); 

    auto begin = high_resolution_clock::now();
    
    simulation(Input);
    
    auto end = high_resolution_clock::now();
    program::printElapsed(begin, end);

    return(0);     
}

void simulation(sysInput *Input)
{
    SimBox *BOX = new SimBox;
    atom_style *ATOMS = new atom_style[Input->N];
    interactions ***INTERACTIONS = program::createInteractions(2);
    
    INTERACTIONS[1][1] = new interactions("WCA_2P", {1.0, 1.0, 1.122462048});
    INTERACTIONS[2][2] = new interactions("WCA_2P", {1.0, 1.0, 1.122462048});
    INTERACTIONS[1][2] = new interactions("WCA_2P", {1.0, 1.0, 1.122462048});

    program::mirrorInteractions(INTERACTIONS, 2);

    //printf("%f\n", any_cast<WCA_2P>(INTERACTIONS[2][1]->pair_style).get_forces(1.2599)[0]);
    //printf("%f\n", INTERACTIONS[1][1]->getrc());

#ifdef INIT_CONFIG_FILE
    char *init_fname = {INIT_CONFIG_FILE};
    BOX -> initBox(ATOMS, BOX, INTERACTIONS, Input, init_fname); 
#else
    BOX -> initBox(ATOMS, BOX, INTERACTIONS, Input);
#endif

    char *fname = new char[50];
    
    program::makeFolder(Input);
    program::writeLog(Input);

    sprintf(fname, "../Data%d/init.xyz", Input->nr);
    program::write2xyz(ATOMS, Input, 0, fname);

    sprintf(fname, "../Data%d/frame0.dat", Input->nr);
    program::writeFrame(ATOMS, Input, fname);

    printf("%f\n", BOX->pe);

    // runNVE params - runID time dt thermo_every traj_every norm 
    // runLangevin params - runID time dt thermo_every traj_every norm zero kmc
    // runKMC params - rate bias delay verbose

    /*
    runNVE *RUN1 = new runNVE(1, 10, 5e-4, 1000, 1000, true);
    RUN1 -> integrateNVE(ATOMS, BOX, INTERACTION, Input);

    runLangevin *RUN2 = new runLangevin(2, 1000, 5e-4, 1000, 1000, true, true, false);
    RUN2 -> integrateLangevin(ATOMS, BOX, INTERACTION, Input);
    
    runLangevin *RUN3 = new runLangevin(3, 100000, 5e-4, 10000, 10000, true, true, false);
    KMC_poisson *KMC;
    
    if(RUN3 -> kmc == true)
    { 
        KMC = new KMC_poisson(5e-4, 1.0, 0, false);
        RUN3 -> integrateLangevin(ATOMS, BOX, INTERACTION, Input, KMC);
    }
    else
        RUN3 -> integrateLangevin(ATOMS, BOX, INTERACTION, Input);

    delete[] ATOMS;
    delete BOX;
    */
}