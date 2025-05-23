// system.cpp

#include "library.h" 
#include "fileIO.h"
#include "compute.h"
#include "random.h"
#include "system.h"

using namespace program;

program::SimBox::SimBox(){
	dim = 2;
    nAtoms = 1;
    nAtomTypes = 1;
    boxLength_x = 1.0;
    boxLength_y = boxLength_x;
    numFrac = 1.0;
    ke = 0.0;
    pe = 0.0;
    etot = 0.0;
    SVx = 0.0;
    SVy = 0.0;
    temp = 0.0;
    
    rcell_x = 1.122462048;
    rcell_y = rcell_x;
    Ncell_x = int(boxLength_x/rcell_x);
    Ncell_y = int(boxLength_y/rcell_y);

    for(int i = 0; i < MAXCELL; i++)
	{
		MAPS[i] = 0;
		LIST[i] = 0;
		HEAD[i] = 0;
	}
}

program::SimBox::~SimBox(){}

void program::SimBox::initBox(atom_style *ATOMS, SimBox *BOX, interactions ***INTERACTIONS, sysInput *Input, char *fname)
{   
    nAtoms = int(Input->N);
    nAtomTypes = Input->nAtomTypes;
    boxLength_x = Input->L;
    boxLength_y = Input->S*boxLength_x;
    numFrac = Input->pfrac;
    
    if(fname != NULL)
    {
    	program::readConfigFile(ATOMS, BOX, Input, fname);
    	assignProperties(ATOMS, Input, false);
    }
    else 
    {
        setRandomConfig(ATOMS);
        //setRandomRegionConfig(ATOMS, Input, 0.5*boxLength_x - 15.0, 0.5*boxLength_x + 15.0, 0.0, boxLength_y);
        assignProperties(ATOMS, Input, false);
        assignVelocities(ATOMS);
    }

    if(INTERACTIONS != NULL)
    {
    	rcell_x = program::getmaxrc(INTERACTIONS, nAtomTypes);
    	rcell_y = rcell_x;
    	buildCellMaps();
    
    	program::computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS);	
    }
    
    program::computeKineticEnergy(ATOMS, BOX);
    program::computeTemperature(BOX);
}

void program::SimBox::checkMinImage(float *dx, float *dy)
{
	if(dx != NULL)
	{
		if(*dx >= 0.5*boxLength_x) *dx -= boxLength_x;
		else if(*dx <= -0.5*boxLength_x) *dx += boxLength_x;
	}
	if(dy != NULL)
	{
		if(*dy >= 0.5*boxLength_y) *dy -= boxLength_y;
		else if(*dy <= -0.5*boxLength_y) *dy += boxLength_y; 
	}
}

void program::SimBox::checkPBC(atom_style *ATOMS)
{
	for(int i = 0; i < nAtoms; i++)
	{
		if(ATOMS[i].rx <= 0.0)
		{
			ATOMS[i].rx += boxLength_x;
			ATOMS[i].jumpx -= 1;
		}
		else if(ATOMS[i].rx >= boxLength_x)
		{
			ATOMS[i].rx -= boxLength_x;
			ATOMS[i].jumpx += 1;
		}
		if(ATOMS[i].ry <= 0.0)
		{
			ATOMS[i].ry += boxLength_y;
			ATOMS[i].jumpy -= 1; 
		}
		else if(ATOMS[i].ry >= boxLength_y)
		{
			ATOMS[i].ry -= boxLength_y;
			ATOMS[i].jumpy += 1;
		}
	}
}

void program::SimBox::setRandomConfig(atom_style *ATOMS)
{
	float zahl1;
	long idum;
	program::URN(&zahl1, &idum, 1);

	int i = 0;
	while(i < nAtoms)
	{
		int accept = 1;

		program::URN(&zahl1, &idum);
		ATOMS[i].rx = boxLength_x*zahl1;
		program::URN(&zahl1, &idum);
		ATOMS[i].ry = boxLength_y*zahl1;

		for(int j = 0; j < i; j++)
		{
			float dx = ATOMS[i].rx - ATOMS[j].rx;
			float dy = ATOMS[i].ry - ATOMS[j].ry;
			checkMinImage(&dx, &dy);

			if(sqrt(dx*dx + dy*dy) < 1.0)
			{
				accept = 0;
				break;
			}
		}

		if(accept == 1) i++;			
	}	
	printf("Generated Random configuration for %d atoms.\n", nAtoms);	
}

void program::SimBox::setRandomRegionConfig(atom_style *ATOMS, sysInput *Input, float xmin, float xmax, float ymin, float ymax)
{
	float zahl1;
	long idum;
	program::URN(&zahl1, &idum, 1);

	int i = 0;
	while(i < Input->PR*nAtoms)
	{
		int accept = 1;

		program::URN(&zahl1, &idum);
		ATOMS[i].rx = xmin + zahl1*(xmax - xmin);
		program::URN(&zahl1, &idum);
		ATOMS[i].ry = ymin + zahl1*(ymax - ymin);

		for(int j = 0; j < i; j++)
		{
			float dx = ATOMS[i].rx - ATOMS[j].rx;
			float dy = ATOMS[i].ry - ATOMS[j].ry;
			checkMinImage(&dx, &dy);

			if(sqrt(dx*dx + dy*dy) < 0.8)
			{
				accept = 0;
				break;
			}
		}

		if(accept == 1) i++;
	}

	while(i < nAtoms)
	{
		int accept = 1;

		program::URN(&zahl1, &idum);
		if(i < (1 + 0.5)*Input->PR*nAtoms)
			ATOMS[i].rx = 0.0 + zahl1*(xmin - 0.0);
		else
			ATOMS[i].rx = xmax + zahl1*(boxLength_x - xmax);

		program::URN(&zahl1, &idum);
		ATOMS[i].ry = ymin + zahl1*(ymax - ymin);

		for(int j = 0; j < i; j++)
		{
			float dx = ATOMS[i].rx - ATOMS[j].rx;
			float dy = ATOMS[i].ry - ATOMS[j].ry;
			checkMinImage(&dx, &dy);

			if(sqrt(dx*dx + dy*dy) < 0.8)
			{
				accept = 0;
				break;
			}
		}

		if(accept == 1) i++;
	}
	printf("Generated Random configuration for %d atoms.\n", nAtoms);
}

void program::SimBox::assignProperties(atom_style *ATOMS, sysInput *Input, bool random)
{
	int numA = 0, numB = 0;

	if(random == true)
	{
		float zahl1;
		long idum;
		program::URN(&zahl1, &idum, 2);
	
		int i = 0;
		while(i < nAtoms)
		{
			int accept = 0;
			program::URN(&zahl1, &idum);
		
			if((numA < nAtoms-int(Input->PR*nAtoms)) and (zahl1 < Input->PR or numB >= int(Input->PR*nAtoms)))
			{
				ATOMS[i].id = 'N';
				ATOMS[i].type = 1;
				ATOMS[i].Pe = Input->PeA;
				numA++;
				accept = 1;
			}

			else
			{
				ATOMS[i].id = 'O';
				ATOMS[i].type = 2;
				ATOMS[i].Pe = Input->PeB;
				numB++;
				accept = 1;
			}

			if(accept == 1) i++;
		}
	}

	else
	{
		for(int i = 0; i < Input->PR*nAtoms; i++)
		{
			ATOMS[i].id = 'O';
			ATOMS[i].type = 2;
			ATOMS[i].Pe = Input->PeB;
			numB++;
		}

		for(int i = Input->PR*nAtoms; i < nAtoms; i++)
		{
			ATOMS[i].id = 'N';
			ATOMS[i].type = 1;
			ATOMS[i].Pe = Input->PeA;
			numA++;
		}
	}

	for(int i = 0; i < nAtoms; i++)
	{
		ATOMS[i].D = Input->D_str;
		ATOMS[i].m = Input->m_str;
	}
	
	printf("Created %d atoms of type A and %d atoms of type B\n", numA, numB);
}

void program::SimBox::assignVelocities(atom_style *ATOMS)
{
	float gauss1, gauss2;
	long idum;
	program::GAUSS(&gauss1, &gauss2, &idum, 3);

	for(int i = 0; i < nAtoms; i++)
	{
		program::GAUSS(&gauss1, &gauss2, &idum);
		ATOMS[i].vx = sqrt(1.0/ATOMS[i].m)*gauss1;
		ATOMS[i].vy = sqrt(1.0/ATOMS[i].m)*gauss2;

		SVx += ATOMS[i].vx;
		SVy += ATOMS[i].vy;
	}

	for(int i = 0; i < nAtoms; i++)
	{
		ATOMS[i].vx -= (SVx/nAtoms);
		ATOMS[i].vy -= (SVy/nAtoms);
	}
	printf("Assigned velocities for %d atoms.\n", nAtoms);
}

int program::SimBox::cellindex(int ix, int iy)
{   
    if (ix >= Ncell_x) ix = ix - Ncell_x; 
	else if (ix <= -1) ix = ix + Ncell_x; 
	if (iy >= Ncell_y) iy = iy - Ncell_y;
	else if (iy <= -1) iy = iy + Ncell_y;

	if(Ncell_x >= Ncell_y)
		return (1 + ix + iy*Ncell_x);
	else
		return (1 + ix + iy*Ncell_y);
}

void program::SimBox::buildCellMaps()
{
    Ncell_x = int(boxLength_x/rcell_x);
    Ncell_y = int(boxLength_y/rcell_y);

    rcell_x = float(boxLength_x/Ncell_x);
    rcell_y = float(boxLength_y/Ncell_y);

	for (int ix = 0; ix < Ncell_x; ix++)
	{
		for (int iy = 0; iy < Ncell_y; iy++)
		{
			int imap = 4*(cellindex(ix,iy)-1);
			MAPS[imap+1] = cellindex(ix-1,iy);
			MAPS[imap+2] = cellindex(ix-1,iy+1);
			MAPS[imap+3] = cellindex(ix,iy+1);
			MAPS[imap+4] = cellindex(ix+1,iy+1);
		}
	}
    
    printf("Successfully constructed MAPS array with %d cells.\n", int(Ncell_x * Ncell_y));
}

void program::SimBox::buildCellList(atom_style *ATOMS)
{	
	for(int i = 0; i < MAXCELL; i++)
	{
		LIST[i] = 0;
		HEAD[i] = 0;
	}

	for (int i = 1; i <= nAtoms; i++)
	{
		int ii = i - 1;
		int ix = int(ATOMS[ii].rx/rcell_x);
		int iy = int(ATOMS[ii].ry/rcell_y);
		
		int icell = cellindex(ix, iy);	        
		LIST[i] = HEAD[icell];
		HEAD[icell] = i;		    
	}	
}

void program::SimBox::getBrownianForce(atom_style *ATOMS, bool zero, int step)
{
	float gauss1, gauss2;
	long idum;
	program::GAUSS(&gauss1, &gauss2, &idum, 4 + step);

	for(int i = 0; i < nAtoms; i++)
	{
		program::GAUSS(&gauss1, &gauss2, &idum);
		ATOMS[i].bfx = gauss1;
		ATOMS[i].bfy = gauss2;
	}

	if(zero == true)
	{
		float Sbfx = 0.0, Sbfy = 0.0;

		for(int i = 0; i < nAtoms; i++)
		{
			Sbfx += ATOMS[i].bfx;
			Sbfy += ATOMS[i].bfy;
		}

		for(int i = 0; i < nAtoms; i++)
		{
			ATOMS[i].bfx -= (Sbfx/nAtoms);
			ATOMS[i].bfy -= (Sbfy/nAtoms); 
		}
	}
}

void program::SimBox::addForce_x(atom_style *ATOMS, float field_loc_x)
{
	int sign_x;
	for(int i = 0; i < nAtoms; i++)
	{
		if(field_loc_x < ATOMS[i].rx) sign_x = -1;
		else sign_x = 1;

		ATOMS[i].fx += sign_x * ATOMS[i].Pe;
	}
}

void program::SimBox::slice(atom_style *ATOMS, float slice_x, float slice_y)
{
	if(slice_x == 0.0)
	{
		for(int i = 0; i < nAtoms; i++)
		{
			if(ATOMS[i].ry <= slice_y) ATOMS[i].id = 'N';
			else if(ATOMS[i].ry >= slice_y) ATOMS[i].id = 'O';	
		}
	}

	else if(slice_y == 0.0)
	{
		for(int i = 0; i < nAtoms; i++)
		{
			if(ATOMS[i].rx <= slice_x) ATOMS[i].id = 'N';
			else if(ATOMS[i].rx >= slice_x) ATOMS[i].id = 'O';	
		}
	}
}