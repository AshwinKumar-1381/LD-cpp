// system.h

#include "particles.h"

using namespace program;

namespace program {
    
    class SimBox {
    
    public:
    
    int nAtoms;
    float boxLength;
    float numFrac;
    float ke, pe;
    
    SimBox();
    void initBox(Atoms2D* ATOMS, float pfrac, float L, int N, char* fname);
    void buildCellList();
    void maps();
    void setRandomConfig();
    void readConfigFile(Atoms2D* ATOMS, char* fname);
    
    };
}
