// main.cpp

#include "library.h"
#include "system.h"
#include "fileIO.h"
#include "run.h"
#include "kmc.h"

//#define INIT_CONFIG_FILE "../tests/frame0.dat"

using namespace program;

void simulation(sysInput *Input);

int main(int argc, char *argv[])
{
    program::sysInput *Input = new sysInput;

    Input -> nr = atoi(argv[1]);
    Input -> lj_ep = atof(argv[2]);
    Input -> m_str = atof(argv[3]);
    Input -> pfrac = atof(argv[4]);
    Input -> L = atof(argv[5]);
    Input -> PR = atof(argv[6]);
    Input -> PeA = atof(argv[7]);
    Input -> PeB = atof(argv[8]);
    
    Input -> N = int(Input->pfrac*Input->L*Input->L); 

    simulation(Input);
    return(0);     
}

void simulation(sysInput *Input)
{
    SimBox *BOX = new SimBox;
    atom_style *ATOMS = new atom_style[Input->N];
    
    pair_style *INTERACTION = new pair_style;
    INTERACTION->epsilon = Input->lj_ep;

#ifdef INIT_CONFIG_FILE
    char *init_fname = {INIT_CONFIG_FILE};
    BOX->initBox(ATOMS, BOX, INTERACTION, Input, init_fname); 
#else
    BOX->initBox(ATOMS, BOX, INTERACTION, Input);
#endif

    char *fname = new char[50];
    
    program::makeFolder(Input);
    program::writeLog(Input);

    sprintf(fname, "../Data%d/init.xyz", Input->nr);
    program::write2xyz(ATOMS, Input, 0, fname);

    sprintf(fname, "../Data%d/frame0.dat", Input->nr);
    program::writeFrame(ATOMS, Input, fname);

    // runNVE params - runID time dt thermo_every traj_every norm 
    runNVE *RUN1 = new runNVE(1, 10, 5e-4, 2000, 2000, true);
    RUN1 -> integrateNVE(ATOMS, BOX, INTERACTION, Input);

    // runLangevin params - runID time dt thermo_every traj_every norm zero kmc
    runLangevin *RUN2 = new runLangevin(2, 500, 5e-4, 10000, 500, true, true, true);
    KMC_poisson *KMC;
    
    if(RUN2->kmc == true)
    { 
        // runKMC params - rate bias verbose
        KMC = new KMC_poisson(0.05, 0.5, false);
        RUN2 -> integrateLangevin(ATOMS, BOX, INTERACTION, Input, KMC);
    }
    else
        RUN2 -> integrateLangevin(ATOMS, BOX, INTERACTION, Input);

    delete[] ATOMS;
    delete BOX;
}