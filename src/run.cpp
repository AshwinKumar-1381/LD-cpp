// run.cpp

#include "library.h"
#include "run.h"
#include "compute.h"
#include "fileIO.h"
#include "random.h"

using namespace program;

program::runNVE::runNVE(int id, char run_name[50], float t, float delta_t, int thermo_val, int traj_val, bool norm_val)
{
	runID = id;
	name = {run_name};
	time = t;
	dt = delta_t;
	thermo_every = thermo_val;
	traj_every = traj_val;
	norm = norm_val;

	maxSteps = ceil(time/dt);
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
		 	printf("\n--- %s ---\n", name);
		 	printf("\nstep pe ke etot\n");
		}

		if(step % thermo_every == 0)
		{
			printf("%d %f %f %f\n", step, BOX->pe/fac, BOX->ke/fac, BOX->etot/fac);
			program::writeThermo(BOX, Input, runID, fac, step);
		}

		if(step % traj_every == 0) 
			program::write2traj(ATOMS, Input, runID, step);

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			float dt2m = float(0.5*dt/ATOMS[i].m);

			ATOMS[i].v2x = ATOMS[i].vx + dt2m*ATOMS[i].fx;
			ATOMS[i].v2y = ATOMS[i].vy + dt2m*ATOMS[i].fy;

			ATOMS[i].rx = ATOMS[i].rx + dt*ATOMS[i].v2x;
			ATOMS[i].ry = ATOMS[i].ry + dt*ATOMS[i].v2y;
		}	

		BOX -> checkPBC(ATOMS);

		if(INTERACTIONS != NULL)
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

program::runLangevin::runLangevin(int id, char run_name[50], float t, float delta_t, int thermo_val, int traj_val, bool norm_val, bool zero_val, bool kmc_val, float field_x)
{
	runID = id;
	name = {run_name};
	time = t;
	dt = delta_t;
	thermo_every = thermo_val;
	traj_every = traj_val;
	norm = norm_val;
	zero = zero_val;
	kmc = kmc_val;
	field_loc_x = field_x;

	maxSteps = ceil(time/dt);
}

void program::runLangevin::integrateLangevin(atom_style *ATOMS, SimBox *BOX, interactions ***INTERACTIONS, sysInput *Input, KMC_poisson *KMC)
{
	int fac;
	if(norm == true) fac = BOX->nAtoms;
	else fac = 1;

	if(Input->restart == true)
	{
		this->res_step = Input->res_step;
		maxSteps += res_step;
	}

	if(kmc == true) 
	{
		if(KMC -> dist == true)
			KMC -> initialize(Input, dt, ATOMS, BOX);
		else
			KMC -> initialize(Input, dt);
	}

	float d2t = 0.5*dt;

	for(int step = res_step; step <= maxSteps; step++)
	{
		if(step == res_step)
		{	
		 	printf("\n--- %s ---\n", name);
		 	printf("\nstep pe ke etot temp\n");

		 	if(INTERACTIONS != NULL)
		 		BOX -> addForce_x(ATOMS, field_loc_x);
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
			float c1 = exp(-1.0*dtm/ATOMS[i].D);
			float c2 = sqrt((1.0 - c1*c1)/ATOMS[i].m);

			ATOMS[i].v2x = ATOMS[i].vx + 0.5*dtm*ATOMS[i].fx;
			ATOMS[i].v2y = ATOMS[i].vy + 0.5*dtm*ATOMS[i].fy;

			ATOMS[i].r2x = ATOMS[i].rx + d2t*ATOMS[i].v2x;
			ATOMS[i].r2y = ATOMS[i].ry + d2t*ATOMS[i].v2y;

			ATOMS[i].v2x = ATOMS[i].v2x*c1 + ATOMS[i].bfx*c2;
			ATOMS[i].v2y = ATOMS[i].v2y*c1 + ATOMS[i].bfy*c2;

			ATOMS[i].rx = ATOMS[i].r2x + d2t*ATOMS[i].v2x;
			ATOMS[i].ry = ATOMS[i].r2y + d2t*ATOMS[i].v2y;
		}	

		BOX -> checkPBC(ATOMS);

		if(INTERACTIONS != NULL)
		{
			program::computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS, false);
			BOX -> addForce_x(ATOMS, field_loc_x);
		}

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

		 	if(INTERACTIONS != NULL)
		 		BOX -> addForce_x(ATOMS, field_loc_x);
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
			ATOMS[i].rx = ATOMS[i].rx + ATOMS[i].fx*dt + c1*ATOMS[i].bfx;
			ATOMS[i].ry = ATOMS[i].ry + ATOMS[i].fy*dt + c1*ATOMS[i].bfy;
		}

		BOX -> checkPBC(ATOMS);

		if(INTERACTIONS != NULL)
		{
			program::computeNonBondedInteractions(ATOMS, BOX, INTERACTIONS, false);
			BOX -> addForce_x(ATOMS, field_loc_x);
		}
	}
}