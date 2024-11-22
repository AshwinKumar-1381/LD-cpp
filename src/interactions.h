// interactions.h

#ifndef INTERACTIONS_H
#define INTERACTIONS_H

using namespace program;

namespace program{

    typedef class WCA{
    
    public:
    
    float epsilon;
    float sigma;
    
    WCA();
    ~WCA();
    float energy(float r2ij);
    float force(float r2ij);
    
    }pair_style;
}

#endif /*INTERACTIONS_H*/