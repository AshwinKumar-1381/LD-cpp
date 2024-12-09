// interactions.h

#ifndef INTERACTIONS_H
#define INTERACTIONS_H

//using namespace program;

namespace program{

    typedef class WCA{
    
    public:
    
    float epsilon;
    float sigma;
    float rcut;
    
    WCA(float eps = 1.0, float sig = 1.0, float cut = 1.122462048);
    ~WCA();
    float energy(float r2ij);
    float force(float r2ij);
    
    }pair_style;
}

#endif /*INTERACTIONS_H*/