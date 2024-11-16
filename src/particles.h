// particles.h

#ifndef PARTICLES_H
#define PARTICLES_H

using namespace program;

namespace program{
    
    typedef class Atoms2D{
    
    public:
    
    float rx, ry, rz;
    float r2x, r2y;
    float px, py;
    float p2x, p2y;
    float fx, fy;
    float bfx, bfy;
    int jumpx, jumpy;
    float m;
    float D;
    float Pe;
    char id;
    
    Atoms2D();
    ~Atoms2D();
    
    }atom_style;
}

#endif /*PARTICLES_H*/