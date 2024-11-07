// main.cpp

#include "library.h"
#include "system.h"
#include "fileIO.h"
#include "run.h"

#define INIT_CONFIG_FILE "../tests/frame0.dat"

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
    runNVE *run1 = new runNVE(1, 50, 1e-6, 50000, 1000, true);
    run1 -> integrateNVE(ATOMS, BOX, INTERACTION, Input);

    delete BOX;
    delete[] ATOMS;
}