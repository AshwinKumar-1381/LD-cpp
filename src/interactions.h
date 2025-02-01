// interactions.h

#ifndef INTERACTIONS_H
#define INTERACTIONS_H

namespace program{

    class WCA_2P
    {
        public:
    
        float epsilon;
        float sigma;
        float rcut;
    
        WCA_2P(float eps = 1.0, float sig = 1.0, float cut = 1.122462048);
        ~WCA_2P();
        float* get_forces(float r2ij);
    };

    class LJ_2P
    {
        public:

        float epsilon;
        float sigma;
        float rcut;

        LJ_2P(float eps = 1.0, float sig = 1.0, float cut = 1.122462048);
        ~LJ_2P();
        float* get_forces(float r2ij);
    };

    class interactions
    {
        public:
            
        const char *intType;
        any pair_style;

        interactions(const char *Type, vector<float> params);
        ~interactions();
        float getrc();
        float* getForce(float r2ij);
    };

    interactions*** createInteractions(int nAtomTypes);
    void mirrorInteractions(interactions ***Int, int nAtomTypes);
}

#endif /*INTERACTIONS_H*/