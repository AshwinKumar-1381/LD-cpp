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

	float d2t = 0.5*dt;

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

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			float dt2m = float(0.5*dt/ATOMS[i].m);

			ATOMS[i].v2x = ATOMS[i].vx + dt2m*ATOMS[i].fx;
			ATOMS[i].v2y = ATOMS[i].vy + dt2m*ATOMS[i].fy;

			ATOMS[i].rx = ATOMS[i].rx + dt*ATOMS[i].v2x;
			ATOMS[i].ry = ATOMS[i].ry + dt*ATOMS[i].v2y;
		}	

		BOX -> checkPBC(ATOMS);
		program::computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS, false);

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			float dt2m = float(0.5*dt/ATOMS[i].m);

			ATOMS[i].vx = ATOMS[i].v2x + dt2m*ATOMS[i].fx;
			ATOMS[i].vy = ATOMS[i].v2y + dt2m*ATOMS[i].fy;
		}

		program::computeKineticEnergy(ATOMS, BOX);	
	}
}

program::runLangevin_D::runLangevin_D(int id, float t, float delta_t, int thermo_val, int traj_val, bool norm_val, bool zero_val, bool kmc_val, float field_x)
{
	runID = id;
	time = t;
	dt = delta_t;
	thermo_every = thermo_val;
	traj_every = traj_val;
	norm = norm_val;
	zero = zero_val;
	kmc = kmc_val;
	field_loc_x = field_x;

	maxSteps = ceil((time + dt)/dt);
}

void program::runLangevin_D::integrateLangevin(atom_style *ATOMS, SimBox *BOX, interactions ***INTERACTIONS, sysInput *Input, KMC_poisson *KMC)
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
			float dtm = float(dt/ATOMS[i].m);
			float c1 = exp(-1.0*dtm);
			float c2 = sqrt((1.0 - c1*c1)/ATOMS[i].m);
			float c3 = 1.0 - c1;
			
			long sign_x = long((field_loc_x - ATOMS[i].rx)/abs(field_loc_x - ATOMS[i].rx));

			ATOMS[i].v2x = ATOMS[i].vx + 0.5*dtm*ATOMS[i].fx;
			ATOMS[i].v2y = ATOMS[i].vy + 0.5*dtm*ATOMS[i].fy;

			ATOMS[i].r2x = ATOMS[i].rx + d2t*ATOMS[i].v2x;
			ATOMS[i].r2y = ATOMS[i].ry + d2t*ATOMS[i].v2y;

			ATOMS[i].v2x = ATOMS[i].v2x*c1 + ATOMS[i].bfx*c2 + sign_x*ATOMS[i].Pe*c3;
			ATOMS[i].v2y = ATOMS[i].v2y*c1 + ATOMS[i].bfy*c2;

			ATOMS[i].rx = ATOMS[i].r2x + d2t*ATOMS[i].v2x;
			ATOMS[i].ry = ATOMS[i].r2y + d2t*ATOMS[i].v2y;
		}	

		BOX -> checkPBC(ATOMS);
		program::computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS, false);

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			float dtm = float(dt/ATOMS[i].m);

			ATOMS[i].vx = ATOMS[i].v2x + 0.5*dtm*ATOMS[i].fx;
			ATOMS[i].vy = ATOMS[i].v2y + 0.5*dtm*ATOMS[i].fy;
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
	bool norm_val, bool zero_val, bool kmc_val, float field_x)
{
	runID = id;
	time = t;
	dt = delta_t;
	thermo_every = thermo_val;
	traj_every = traj_val;
	norm = norm_val;
	zero = zero_val;
	kmc = kmc_val;
	field_loc_x = field_x;

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
			long sign_x = long((field_loc_x - ATOMS[i].rx)/abs(field_loc_x - ATOMS[i].rx));
		
			ATOMS[i].rx = ATOMS[i].rx + (ATOMS[i].fx + sign_x*ATOMS[i].Pe)*dt + c1*ATOMS[i].bfx;
			ATOMS[i].ry = ATOMS[i].ry + ATOMS[i].fy*dt + c1*ATOMS[i].bfy;
		}

		BOX -> checkPBC(ATOMS);
		program::computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS, false);
	}
}