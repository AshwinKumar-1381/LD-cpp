// interactions.cpp

#include "library.h"
#include "interactions.h"

using namespace program;

program::WCA::WCA()
{
	epsilon = 1.0;
	sigma = 1.0;
}

float program::WCA::energy(float r2ij)
{
	float rr2 = float(sigma*sigma/r2ij);
	float rr6 = float(rr2*rr2*rr2); 	
	return(4*epsilon*rr6*(rr6 - 1.0) + epsilon);
}

float program::WCA::force(float r2ij)
{
	float rr2 = float(sigma*sigma/r2ij);
	float rr6 = float(rr2*rr2*rr2);
	return(48.0*epsilon*rr2*rr6*(rr6 - 0.5));
}