// run.cpp

#include "library.h"
#include "run.h"
#include "compute.h"

using namespace program;

program::runNVE::runNVE()
{
	runID = 0;
	time = 1.0;
	dt = 1.0;
	maxSteps = int(time/dt);
}

void program::runNVE::integrateNVE(atom_style *ATOMS, SimBox *BOX, pair_style *INTERACTION, sysInput *Input)
{
	time = Input->runtime1;
	dt = Input->dt1;
	maxSteps = int(time/dt);

	for(int step = 0; step <= maxSteps; step++)
	{
		
		float d2t = 0.5*dt;

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			float dtm = float(dt/ATOMS[i].m);

			ATOMS[i].p2x = ATOMS[i].px + d2t*ATOMS[i].fx;
			ATOMS[i].p2y = ATOMS[i].py + d2t*ATOMS[i].fy;

			ATOMS[i].rx = ATOMS[i].rx + dtm*ATOMS[i].p2x;
			ATOMS[i].ry = ATOMS[i].ry + dtm*ATOMS[i].p2y;
		}	

		program::computeNonBondedInteractions(ATOMS, BOX, INTERACTION);

		for(int i = 0; i < BOX->nAtoms; i++)
		{
			ATOMS[i].px = ATOMS[i].p2x + d2t*ATOMS[i].fx;
			ATOMS[i].py = ATOMS[i].p2y + d2t*ATOMS[i].fy;
		}

		program::computeKineticEnergy(ATOMS, BOX);

		if(step%100 == 0){
			printf("%d %f %f %f\n", step, BOX->ke/BOX->nAtoms, BOX->pe/BOX->nAtoms, BOX->etot/BOX->nAtoms);
			//printf("Step no. %d, ke = %f, pe = %f, etot = %f\n", step, BOX->ke/BOX->nAtoms, BOX->pe/BOX->nAtoms, BOX->etot/BOX->nAtoms);
		}
	}

}