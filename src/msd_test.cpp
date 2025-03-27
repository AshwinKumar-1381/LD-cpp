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

    Input -> nr         = 278;
    Input -> m_str      = 1.0;
    Input -> D_str      = 1.0;
    Input -> dt         = 5e-5;
    Input -> pfrac      = 0.1;
    Input -> L          = 100.0;
    Input -> S          = 1.0;

    program::simulation(Input);
    return(0);     
}

void program::simulation(sysInput *Input)
{
    Input -> begin = high_resolution_clock::now();

    Input -> N = int(Input->pfrac*Input->L*Input->L*Input->S); 

    SimBox *BOX = new SimBox;
    atom_style *ATOMS = new atom_style[Input->N];

#ifdef INIT_CONFIG_FILE
    char *init_fname = {INIT_CONFIG_FILE};
    BOX -> initBox(ATOMS, BOX, NULL, Input, init_fname); 
#else
    BOX -> initBox(ATOMS, BOX, NULL, Input);
#endif

    char *fname = new char[50];
    
    program::makeFolder(Input);

    sprintf(fname, "../Data%d/init.xyz", Input->nr);
    program::write2xyz(ATOMS, Input, 0, fname);

    sprintf(fname, "../Data%d/frame0.dat", Input->nr);
    program::writeFrame(ATOMS, Input, fname);

    // runNVE params - runID time dt thermo_every traj_every norm 
    // runLangevin params - runID time dt thermo_every traj_every norm zero kmc field_x

    runNVE *RUN1 = new runNVE(1, 1, Input->dt, 10000, 10000, true);
    RUN1 -> integrateNVE(ATOMS, BOX, NULL, Input);
    
    run_style *RUN2 = new run_style(2, 10000, Input->dt, 200000, 200000, true, true, false, 0.5*BOX->boxLength_x);
    RUN2 -> integrateLangevin(ATOMS, BOX, NULL, Input);

    Input -> end = high_resolution_clock::now();
    printf("\n%s", program::returnElapsedTime(Input));

    program::writeLog(Input, BOX, RUN2);

    delete[] ATOMS;
    delete BOX;
}