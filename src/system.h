// system.h

#ifndef SYSTEM_H
#define SYSTEM_H

#include "particles.h"
#include "interactions.h"

#define MAXCELL 900000

using namespace program;

namespace program {
    
    class SimBox {
    
    public:
    
    int dim;
    int nAtoms;
    float boxLength_x;
    float boxLength_y;
    float numFrac;
    float ke, pe, etot;
    float SPx, SPy;
    float temp;
    
    int Ncell_x, Ncell_y;
    float rcell_x, rcell_y;
    int MAPS[MAXCELL], LIST[MAXCELL], HEAD[MAXCELL];
    
    SimBox();
    ~SimBox();
    void initBox(atom_style* ATOMS, SimBox* BOX, pair_style* INTERACTION, sysInput *Input, char* fname = NULL);
    
    int cellindex(int ix, int iy);
    void buildCellMaps();
    void buildCellList(atom_style* ATOMS);
    
    void checkMinImage(float *dx = NULL, float *dy = NULL); 
    void checkPBC(atom_style *ATOMS);
    void setRandomConfig(atom_style *ATOMS);
    void assignMomenta(atom_style *ATOMS);
    void assignProperties(atom_style *ATOMS, sysInput *Input);
    void getBrownianForce(atom_style *ATOMS, bool zero, int step);
    };
}

#endif /* SYSTEM_H */