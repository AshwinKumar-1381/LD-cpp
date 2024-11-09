// run.cpp

#include "library.h"
#include "run.h"
#include "compute.h"
#include "fileIO.h"

using namespace program;

program::runNVE::runNVE(int id, float t, float delta_t, int every1, int every2, bool val)
{
	runID = id;
	time = t;
	dt = delta_t;
	thermo_every = every1;
	traj_every = every2;
	norm = val;
	maxSteps = int(time/dt);
}

void program::runNVE::integrateNVE(atom_style *ATOMS, SimBox *BOX, pair_style *INTERACTION, sysInput *Input)
{
	maxSteps = ceil(time/dt);

	int fac;
	if(norm == true) fac = BOX->nAtoms;
	else fac = 1;

	for(int step = 0; step <= maxSteps; step++)
	{
		if(step == 0)
		{	
		 	printf("\n--- Run %d ---\n", runID);
		 	printf("\nstep pe ke etot\n");
		}

		if(step % thermo_every == 0)
		{
			printf("%d %f %f %f\n", step, BOX->pe/fac, BOX->ke/fac, BOX->etot/fac);
			program::writeThermo(BOX, Input, runID, fac, step);
		}

		float d2t = 0.5*dt;

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			float dtm = float(dt/ATOMS[i].m);

			ATOMS[i].p2x = ATOMS[i].px + d2t*ATOMS[i].fx;
			ATOMS[i].p2y = ATOMS[i].py + d2t*ATOMS[i].fy;

			ATOMS[i].rx = ATOMS[i].rx + dtm*ATOMS[i].p2x;
			ATOMS[i].ry = ATOMS[i].ry + dtm*ATOMS[i].p2y;
		}	

		BOX -> checkPBC(ATOMS);
		program::computeNonBondedInteractions(ATOMS, BOX, INTERACTION);

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			ATOMS[i].px = ATOMS[i].p2x + d2t*ATOMS[i].fx;
			ATOMS[i].py = ATOMS[i].p2y + d2t*ATOMS[i].fy;
		}

		program::computeKineticEnergy(ATOMS, BOX);	
	}
}

program