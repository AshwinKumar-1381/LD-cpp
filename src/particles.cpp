// particles.cpp

#include "library.h"
#include "particles.h"

using namespace program;

program::Atoms2D::Atoms2D()
{
    rx=0.0, ry=0.0, rz=0.0;
    r2x=0.0, r2y=0.0;
    vx=0.0, vy=0.0;
    v2x=0.0, v2y=0.0;
    bfx=0.0, bfy=0.0;
    fx=0.0, fy=0.0;
    jumpx=0, jumpy=0;
    m=1.0; 
    D=1.0;
    Pe=0.0; 
    id='A';
    type=0;
}

program::Atoms2D::~Atoms2D(){}