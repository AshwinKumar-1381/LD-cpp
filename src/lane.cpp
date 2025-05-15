// main.cpp

#include "library.h"
#include "interactions.h"
#include "system.h"
#include "fileIO.h"
#include "run.h"
#include "kmc.h"

// #define INIT_CONFIG_FILE "../lane/May2025/Data5/frames/frame200000.dat"

using namespace program;

int main(int argc, char *argv[])
{
    program::sysInput *Input = new sysInput;

    sprintf(Input -> folder, "lane/May2025");
    Input -> restart    = false;
    Input -> nr         = 5;
    Input -> nAtomTypes = 2;
    Input -> m_str      = 1.0;
    Input -> D_str      = 1.0e-3;
    Input -> dt         = 5e-4;     
    Input -> pfrac      = 0.45;
    Input -> L          = 150.0;
    Input -> S          = 0.2;
    Input -> PR         = 0.5;
    Input -> PeA        = -20.0;
    Input -> PeB        = +20.0;

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

    // runNVE params - runID time dt thermo_every traj_every norm 
    // runLangevin params - runID time dt thermo_every traj_every norm zero kmc field_x
    // runBrownian params - runID time dt thermo_every traj_every norm zero kmc
    // runKMC params - rate bias delay verbose dist

    if(Input -> restart == false)
    {
        sprintf(fname, "../%s/Data%d/init.xyz", Input->folder, Input->nr);
        program::write2xyz(ATOMS, Input, 0, fname);

        runNVE *RUN1 = new runNVE(1, "Run1", 500, Input->dt, 10000, 10000, true);
        RUN1 -> integrateNVE(ATOMS, BOX, INTERACTIONS, Input);
        
        sprintf(fname, "../%s/Data%d/frames/frame0.dat", Input->folder, Input->nr);
        program::writeFrame(ATOMS, Input, 0, fname);

        run_style *RUN2 = new run_style(2, "Run2", 5e4, Input->dt, 10000, 10000, true, true, true, 1.0*BOX->boxLength_x);
        KMC_poisson *KMC;
        
        if(RUN2 -> kmc == true)
        { 
            KMC = new KMC_poisson(1e-3, 1.0, 0);
            RUN2 -> integrateLangevin(ATOMS, BOX, INTERACTIONS, Input, KMC);
        }
        else
            RUN2 -> integrateLangevin(ATOMS, BOX, INTERACTIONS, Input);   

        sprintf(fname, "../%s/Data%d/frames/frame%ld.dat", Input->folder, Input->nr, long(RUN2->maxSteps));
        program::writeFrame(ATOMS, Input, long(RUN2->maxSteps), fname);

        Input -> end = high_resolution_clock::now();
        printf("\n%s", program::returnElapsedTime(Input));

        if(RUN2 -> kmc == true) program::writeLog(Input, BOX, RUN2, KMC);
        else program::writeLog(Input, BOX, RUN2);
    }

    else if(Input -> restart == true)
    {
        run_style *RES1 = new run_style(3, "Res1", 100, Input->dt, 10000, 10000, true, true, true, 1.0*BOX->boxLength_x);
        KMC_poisson *KMC;

        if(RES1 -> kmc == true)
        {
            KMC = new KMC_poisson(1e-3, 1.0, 0);
            RES1 -> integrateLangevin(ATOMS, BOX, INTERACTIONS, Input, KMC);
        }
        else
            RES1 -> integrateLangevin(ATOMS, BOX, INTERACTIONS, Input);

        sprintf(fname, "../%s/Data%d/frames/frame%ld.dat", Input->folder, Input->nr, long(RES1->maxSteps));
        program::writeFrame(ATOMS, Input, long(RES1->maxSteps), fname);

        Input -> end = high_resolution_clock::now();
        printf("\n%s", program::returnElapsedTime(Input));

        if(RES1 -> kmc == true) program::writeLog(Input, BOX, RES1, KMC);
        else program::writeLog(Input, BOX, RES1);        
    }

    delete[] ATOMS;
    delete BOX;
}