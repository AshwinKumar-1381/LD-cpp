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
    int nAtomTypes;
    float boxLength_x;
    float boxLength_y;
    float numFrac;
    float ke, pe, etot;
    float SVx, SVy;
    float temp;
    
    int Ncell_x, Ncell_y;
    float rcell_x, rcell_y;
    int MAPS[MAXCELL], LIST[MAXCELL], HEAD[MAXCELL];
    
    SimBox();
    ~SimBox();
    void initBox(atom_style* ATOMS, SimBox* BOX, interactions*** INTERACTIONS, sysInput *Input, char* fname = NULL);
    
    int cellindex(int ix, int iy);
    void buildCellMaps();
    void buildCellList(atom_style* ATOMS);
    
    void checkMinImage(float *dx = NULL, float *dy = NULL); 
    void checkPBC(atom_style *ATOMS);
    void setRandomConfig(atom_style *ATOMS);
    void setRandomRegionConfig(atom_style *ATOMS, sysInput *Input, float x_min, float x_max, float y_min, float y_max);
    void assignVelocities(atom_style *ATOMS);
    void assignProperties(atom_style *ATOMS, sysInput *Input, bool random = true);
    void getBrownianForce(atom_style *ATOMS, bool zero, int step);
    };
}

#endif /* SYSTEM_H */