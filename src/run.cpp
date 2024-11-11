// run.cpp

#include "library.h"
#include "run.h"
#include "compute.h"
#include "fileIO.h"

using namespace program;

program::runNVE::runNVE(int id, float t, float delta_t, int thermo_val, int traj_val, bool norm_val)
{
	runID = id;
	time = t;
	dt = delta_t;
	thermo_every = thermo_val;
	traj_every = traj_val;
	norm = norm_val;
	maxSteps = ceil(time/dt);
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

program::runLangevin::runLangevin(int id, float t, float delta_t, int thermo_val, int traj_val, bool norm_val, bool zero_val)
{
	runID = id;
	time = t;
	dt = delta_t;
	thermo_every = thermo_val;
	traj_every = traj_val;
	norm = norm_val;
	zero = zero_val;
	maxSteps = ceil(time/dt);
}

void program::runLangevin::integrateLangevin(atom_style *ATOMS, SimBox *BOX, pair_style *INTERACTION, sysInput *Input)
{
	maxSteps = ceil(time/dt);

	int fac;
	if(norm == true) fac = BOX->nAtoms;
	else fac = 1;

	float d2t = 0.5*dt;
	float dtm = float(dt/ATOMS[0].m);
	float c1 = exp(-1.0*dtm);
	float c2 = sqrt(ATOMS[0].m*(1.0-c1*c1));

	for(int step = 0; step <= maxSteps; step++)
	{
		if(step == 0)
		{	
		 	printf("\n--- Run %d ---\n", runID);
		 	printf("\nstep pe ke etot temp\n");
		}

		if(step % thermo_every == 0)
		{
			printf("%d %f %f %f %f\n", step, BOX->pe/fac, BOX->ke/fac, BOX->etot/fac, BOX->temp);
			program::writeThermo(BOX, Input, runID, fac, step);
		}

		BOX -> getBrownianForce(ATOMS, zero, step);

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			ATOMS[i].p2x = ATOMS[i].px + d2t*ATOMS[i].fx;
			ATOMS[i].p2y = ATOMS[i].py + d2t*ATOMS[i].fy;

			ATOMS[i].r2x = ATOMS[i].rx + 0.5*dtm*ATOMS[i].p2x;
			ATOMS[i].r2y = ATOMS[i].ry + 0.5*dtm*ATOMS[i].p2y;

			ATOMS[i].p2x = ATOMS[i].p2x*c1 + ATOMS[i].bfx*c2;
			ATOMS[i].p2y = ATOMS[i].p2y*c1 + ATOMS[i].bfy*c2;

			ATOMS[i].rx = ATOMS[i].r2x + 0.5*dtm*ATOMS[i].p2x;
			ATOMS[i].ry = ATOMS[i].r2y + 0.5*dtm*ATOMS[i].p2y;
		}	

		BOX -> checkPBC(ATOMS);
		program::computeNonBondedInteractions(ATOMS, BOX, INTERACTION);

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			ATOMS[i].px = ATOMS[i].p2x + d2t*ATOMS[i].fx;
			ATOMS[i].py = ATOMS[i].p2y + d2t*ATOMS[i].fy;
		}

		program::computeKineticEnergy(ATOMS, BOX);	
		program::computeTemperature(BOX);
	}
}