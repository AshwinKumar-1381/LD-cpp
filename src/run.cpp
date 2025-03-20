// run.cpp

#include "library.h"
#include "run.h"
#include "compute.h"
#include "fileIO.h"
#include "random.h"

using namespace program;

program::runNVE::runNVE(int id, float t, float delta_t, int thermo_val, int traj_val, bool norm_val)
{
	runID = id;
	time = t;
	dt = delta_t;
	thermo_every = thermo_val;
	traj_every = traj_val;
	norm = norm_val;

	maxSteps = ceil((time + dt)/dt);
}

void program::runNVE::integrateNVE(atom_style *ATOMS, SimBox *BOX, interactions ***INTERACTIONS, sysInput *Input)
{
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
		program::computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS);

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			ATOMS[i].px = ATOMS[i].p2x + d2t*ATOMS[i].fx;
			ATOMS[i].py = ATOMS[i].p2y + d2t*ATOMS[i].fy;
		}

		program::computeKineticEnergy(ATOMS, BOX);	
	}
}

program::runLangevin::runLangevin(int id, float t, float delta_t, int thermo_val, int traj_val, bool norm_val, bool zero_val, bool kmc_val)
{
	runID = id;
	time = t;
	dt = delta_t;
	thermo_every = thermo_val;
	traj_every = traj_val;
	norm = norm_val;
	zero = zero_val;
	kmc = kmc_val;

	maxSteps = ceil((time + dt)/dt);
}

void program::runLangevin::integrateLangevin(atom_style *ATOMS, SimBox *BOX, interactions ***INTERACTIONS, sysInput *Input, KMC_poisson *KMC)
{
	int fac;
	if(norm == true) fac = BOX->nAtoms;
	else fac = 1;

	if(kmc == true) 
	{
		if(KMC -> dist == true)
			KMC -> initialize(Input, dt, ATOMS, BOX);
		else
			KMC -> initialize(Input, dt);
	}

	float d2t = 0.5*dt;
	float dtm = float(dt/ATOMS[0].m);
	float c1 = exp(-1.0*dtm);
	float c2 = sqrt(ATOMS[0].m*(1.0-c1*c1));
	float c3 = ATOMS[0].m*(1.0 - c1);

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

		if(step % traj_every == 0) 
			program::write2traj(ATOMS, Input, runID, step);
		
		if(kmc == true)
		{
			KMC -> Switch(ATOMS, BOX, Input, dt, step);
			if(step % thermo_every == 0) 
				program::writeKMC(KMC, Input, step);
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

			if(ATOMS[i].rx <= 0.5*BOX->boxLength_x) 
				ATOMS[i].p2x += 1.0*c3*ATOMS[i].Pe;
			else 
				ATOMS[i].p2x += -1.0*c3*ATOMS[i].Pe;

			ATOMS[i].rx = ATOMS[i].r2x + 0.5*dtm*ATOMS[i].p2x;
			ATOMS[i].ry = ATOMS[i].r2y + 0.5*dtm*ATOMS[i].p2y;
		}	

		BOX -> checkPBC(ATOMS);
		program::computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS);

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			ATOMS[i].px = ATOMS[i].p2x + d2t*ATOMS[i].fx;
			ATOMS[i].py = ATOMS[i].p2y + d2t*ATOMS[i].fy;
		}

		program::computeKineticEnergy(ATOMS, BOX);	
		program::computeTemperature(BOX);
	}

	if(kmc == true)
	{
		if(KMC -> dist == true)
		{
			KMC -> dist_msd_tau -> normalize(dt);

			sprintf(KMC->dist_msd_tau->fpathO, "../Data%d/kmcDist.dat", Input->nr);
			KMC -> dist_msd_tau -> write2file();
		}
	}
}

program::runBrownian::runBrownian(int id, float t, float delta_t, int thermo_val, int traj_val,
	bool norm_val, bool zero_val, bool kmc_val)
{
	runID = id;
	time = t;
	dt = delta_t;
	thermo_every = thermo_val;
	traj_every = traj_val;
	norm = norm_val;
	zero = zero_val;
	kmc = kmc_val;

	maxSteps = ceil((time + dt)/dt);
}

void program::runBrownian::integrateBrownian(atom_style *ATOMS, SimBox *BOX, interactions ***INTERACTIONS, sysInput *Input, KMC_poisson *KMC)
{
	int fac;
	if(norm == true) fac = BOX->nAtoms;
	else fac = 1;

	if(kmc == true) 
	{
		if(KMC -> dist == true)
			KMC -> initialize(Input, dt, ATOMS, BOX);
		else
			KMC -> initialize(Input, dt);
	}

	float c1 = sqrt(2.0*dt);

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

		if(step % traj_every == 0) 
			program::write2traj(ATOMS, Input, runID, step);
		
		if(kmc == true)
		{
			KMC -> Switch(ATOMS, BOX, Input, dt, step);
			if(step % thermo_every == 0) 
				program::writeKMC(KMC, Input, step);
		} 
		
		BOX -> getBrownianForce(ATOMS, zero, step);

		for(int i = 0; i < BOX->nAtoms; i++)
		{

			ATOMS[i].r2x = ATOMS[i].rx;
			ATOMS[i].r2y = ATOMS[i].ry;
		
			ATOMS[i].rx = ATOMS[i].r2x + ATOMS[i].fx*dt + c1*ATOMS[i].bfx;
			ATOMS[i].ry = ATOMS[i].r2y + ATOMS[i].fy*dt + c1*ATOMS[i].bfy;

			if(ATOMS[i].r2x <= 0.5*BOX->boxLength_x)
				ATOMS[i].rx += ATOMS[i].Pe*dt;
			else
				ATOMS[i].rx += -1.0*ATOMS[i].Pe*dt;

			// ATOMS[i].px = (ATOMS[i].rx - ATOMS[i].r2x)/dt;
			// ATOMS[i].py = (ATOMS[i].ry - ATOMS[i].r2y)/dt;
		}

		BOX -> checkPBC(ATOMS);
		program::computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS, true);

		// program::computeKineticEnergy(ATOMS, BOX);	
		// program::computeTemperature(BOX);
	}
}