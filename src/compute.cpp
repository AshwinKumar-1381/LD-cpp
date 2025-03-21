// compute.cpp

#include "library.h"
#include "compute.h"

void program::computeNonBondedInteractions(atom_style *ATOMS, SimBox *BOX, interactions ***INTERACTIONS, bool overlaps)
{
	BOX->buildCellList(ATOMS);

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		ATOMS[i].fx = 0.0;
		ATOMS[i].fy = 0.0;
	}

	BOX->pe = 0.0;
	int ncell = int(BOX->Ncell_x * BOX->Ncell_y);

	for(int icell = 1; icell <= ncell; icell++)
	{
		int i = BOX->HEAD[icell];
		while (i != 0)
		{
			int ii = i - 1;
			float rxi = ATOMS[ii].rx;
			float ryi = ATOMS[ii].ry;

			int j = BOX->LIST[i];
			while (j != 0)
			{
				int jj = j - 1;
				float dx = rxi - ATOMS[jj].rx;
				float dy = ryi - ATOMS[jj].ry;
				float r2ij = dx*dx + dy*dy;

				int typei = ATOMS[ii].type;
				int typej = ATOMS[jj].type;

				float rcut = INTERACTIONS[typei][typej]->getrc();

				if(overlaps == true)
				{
					if(r2ij <= 0.5*0.5*rcut*rcut)
						printf("\nOverlap detected between atom %d and atom %d at squared-distance of %f\n", ii, jj, r2ij);
				}

				if(r2ij <= rcut*rcut)
				{
					float *pairs = INTERACTIONS[typei][typej]->getForce(r2ij);  
					float EPOT = pairs[0];
					float F = pairs[1];

					dx *= F;
					dy *= F;

					ATOMS[ii].fx += dx;
					ATOMS[jj].fx -= dx;
					ATOMS[ii].fy += dy;
					ATOMS[jj].fy -= dy;
					BOX->pe += EPOT;

					delete[] pairs;
				}
				
				j = BOX->LIST[j];
			}
			i = BOX->LIST[i];
		}
	}

	int nNbors = 4;
	for(int icell = 1; icell <= ncell; icell++)
	{
		int i = BOX->HEAD[icell];
		while (i != 0)
		{
			int ii = i - 1;
			float rxi = ATOMS[ii].rx;
			float ryi = ATOMS[ii].ry;

			int icell_index = nNbors*(icell - 1);
			for(int nbor = 1; nbor <= nNbors; nbor++)
			{
				int jcell = BOX->MAPS[icell_index + nbor];

				int j = BOX->HEAD[jcell];
				while (j != 0)
				{
					int jj = j - 1;

					float dx = rxi - ATOMS[jj].rx;
					float dy = ryi - ATOMS[jj].ry;

					BOX->checkMinImage(&dx, &dy);
					float r2ij = dx*dx + dy*dy;

					int typei = ATOMS[ii].type;
					int typej = ATOMS[jj].type;

					float rcut = INTERACTIONS[typei][typej]->getrc();

					if(overlaps == true)
					{
						if(r2ij <= 0.5*0.5*rcut*rcut)
							printf("\nOverlap detected between atom %d and atom %d at squared-distance of %f\n", ii, jj, r2ij);
					}

					if(r2ij <= rcut*rcut)
					{
						float *pairs = INTERACTIONS[typei][typej]->getForce(r2ij);  
						float EPOT = pairs[0];
						float F = pairs[1];

						dx *= F;
						dy *= F;

						ATOMS[ii].fx += dx;
						ATOMS[jj].fx -= dx;
						ATOMS[ii].fy += dy;
						ATOMS[jj].fy -= dy;
						BOX->pe += EPOT;

						delete[] pairs;
					}
					j = BOX->LIST[j];
				}
			}
			i = BOX->LIST[i];
		}
	}
}

void program::computeKineticEnergy(atom_style *ATOMS, SimBox *BOX)
{
	BOX->ke = 0.0;
	BOX->etot = 0.0;

	for(int i = 0; i < BOX->nAtoms; i++)
	{
		float ekin = ATOMS[i].vx*ATOMS[i].vx + ATOMS[i].vy*ATOMS[i].vy;
		ekin *= ATOMS[i].m/2.0;
		BOX->ke += ekin; 
	}

	BOX->etot = BOX->pe + BOX->ke;
}

void program::computeTemperature(SimBox *BOX)
{
	BOX->temp = (2.0*BOX->ke)/(BOX->dim*BOX->nAtoms);
}