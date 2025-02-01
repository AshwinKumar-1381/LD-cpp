// interactions.cpp

#include "library.h"
#include "interactions.h"

using namespace program;

interactions*** program::createInteractions(int nAtomTypes)
{
	nAtomTypes += 1;
	interactions ***Int = new interactions**[nAtomTypes];

	for(int i = 0; i < nAtomTypes; i++)
	{
		Int[i] = new interactions*[nAtomTypes];
		for(int j = 0; j < nAtomTypes; j++)
			Int[i][j] = NULL;
	}

	return(Int);
}

void program::mirrorInteractions(interactions ***Int, int nAtomTypes)
{
	for(int i = 1; i < nAtomTypes + 1; i++)
	{
		for(int j = 1; j < nAtomTypes + 1; j++)
		{
			if(Int[i][j] == NULL) Int[i][j] = Int[j][i];
		}
	}	
}

// ------------------ interactions class members ------------------
program::interactions::interactions(const char *Type, vector<float> params)
{
	intType = Type;
    if(strcmp(intType, "WCA_2P") == 0) 
        pair_style = WCA_2P(params[0], params[1], params[2]);

    else if(strcmp(intType, "LJ_2P") == 0) 
        pair_style = LJ_2P(params[0], params[1], params[2]);

    else 
    {
        printf("Pair style %s not found. Exiting ...", Type);
        exit(-1);
    }
}

program::interactions::~interactions(){}

float program::interactions::getrc()
{
	if(strcmp(intType, "WCA_2P") == 0) return(any_cast<WCA_2P>(pair_style).rcut);
	else if(strcmp(intType, "LJ_2P") == 0) return(any_cast<LJ_2P>(pair_style).rcut);
	else
	{
        printf("Pair style %s not found. Exiting ...", intType);
        exit(-1);
    }
}

float* program::interactions::getForce(float r2ij)
{
	float *retF;
	if(strcmp(intType, "WCA_2P") == 0) retF = any_cast<WCA_2P>(pair_style).get_forces(r2ij);
	else if(strcmp(intType, "LJ_2P") == 0) retF = any_cast<LJ_2P>(pair_style).get_forces(r2ij);
	else
	{
        printf("Pair style %s not found. Exiting ...", intType);
        exit(-1);
    }
    return(retF);
}

// ------------------ WCA_2P class members ------------------
program::WCA_2P::WCA_2P(float eps, float sig, float cut)
{
	epsilon = eps;
	sigma = sig;
	rcut = cut;
}

program::WCA_2P::~WCA_2P(){}

float* program::WCA_2P::get_forces(float r2ij)
{
	float rr2 = float(sigma*sigma/r2ij);
	float rr6 = float(rr2*rr2*rr2); 

	float *f = new float[2];
	f[0] = 4*epsilon*rr6*(rr6 - 1.0) + epsilon;		// energy
	f[1] = 48.0*epsilon*rr2*rr6*(rr6 - 0.5);  		// force
	return(f);
}

// ------------------ LJ_2P class members ------------------
program::LJ_2P::LJ_2P(float eps, float sig, float cut)
{
	epsilon = eps;
	sigma = sig;
	rcut = cut;
}

program::LJ_2P::~LJ_2P(){}

float* program::LJ_2P::get_forces(float r2ij){}