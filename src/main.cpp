// main.cpp

#include "library.h"
#include "interactions.h"
#include "system.h"
#include "fileIO.h"
#include "run.h"
#include "kmc.h"

// #define INIT_CONFIG_FILE "../test_force/init3.xyz"

using namespace program;

int main(int argc, char *argv[])
{
    program::sysInput *Input = new sysInput;

    Input -> nr         = 265;
    Input -> nAtomTypes = 2;
    Input -> m_str      = 1.0;
    Input -> dt         = 2e-4;
    Input -> pfrac      = 0.5;
    Input -> L          = 200.0;
    Input -> S          = 0.25;
    Input -> PR         = 0.0;
    Input -> PeA        = -0.0;
    Input -> PeB        = +5.0;

    program::simulation(Input);
    return(0);     
}

void program::simulation(sysInput *Input)
{
    Input -> begin = high_resolution_clock::now();

    Input -> N = int(Input->pfrac*Input->L*Input->L*Input->S); 

    SimBox *BOX = new SimBox;
    atom_style *ATOMS = new atom_style[Input->N];
    interactions ***INTERACTIONS = program::createInteractions(Input->nAtomTypes);
    
    INTERACTIONS[1][1] = new interactions("WCA_2P", {1.0, 1.0, 1.122462048});
    INTERACTIONS[2][2] = new interactions("WCA_2P", {1.0, 1.0, 1.122462048});
    INTERACTIONS[1][2] = new interactions("WCA_2P", {1.0, 1.0, 1.122462048});

    program::mirrorInteractions(INTERACTIONS, Input->nAtomTypes);

#ifdef INIT_CONFIG_FILE
    char *init_fname = {INIT_CONFIG_FILE};
    BOX -> initBox(ATOMS, BOX, INTERACTIONS, Input, init_fname); 
#else
    BOX -> initBox(ATOMS, BOX, INTERACTIONS, Input);
#endif

    char *fname = new char[50];
    
    program::makeFolder(Input);

    sprintf(fname, "../Data%d/init.xyz", Input->nr);
    program::write2xyz(ATOMS, Input, 0, fname);

    sprintf(fname, "../Data%d/frame0.dat", Input->nr);
    program::writeFrame(ATOMS, Input, fname);

    // runNVE params - runID time dt thermo_every traj_every norm 
    // runLangevin params - runID time dt thermo_every traj_every norm zero kmc
    // runBrownian params - runID time dt thermo_every traj_every norm zero kmc
    // runKMC params - rate bias delay verbose dist

    runNVE *RUN1 = new runNVE(1, 1000, Input->dt, 10000, 10000, true);
    RUN1 -> integrateNVE(ATOMS, BOX, INTERACTIONS, Input);
    
    /*
    run_style *RUN2 = new run_style(2, 50000, Input->dt, 10000, 10000, true, true, true);
    KMC_poisson *KMC;
    
    if(RUN2 -> kmc == true)
    { 
        KMC = new KMC_poisson(1e-2, 1.0, 0, false, false);
        RUN2 -> integrateLangevin(ATOMS, BOX, INTERACTIONS, Input, KMC);
    }
    else
        RUN2 -> integrateLangevin(ATOMS, BOX, INTERACTIONS, Input);
    */

    run_style *RUN2 = new run_style(2, 20000, Input->dt, 10000, 10000, true, true, true);
    KMC_poisson *KMC;

    if(RUN2 -> kmc == true)
    {
        KMC = new KMC_poisson(1e-2, 1.0, 0, false, false);
        RUN2 -> integrateBrownian(ATOMS, BOX, INTERACTIONS, Input, KMC);
    }
    else
        RUN2 -> integrateBrownian(ATOMS, BOX, INTERACTIONS, Input);

    Input -> end = high_resolution_clock::now();
    printf("\n%s", program::returnElapsedTime(Input));

    program::writeLog(Input, BOX, RUN2, KMC);

    delete[] ATOMS;
    delete BOX;
}