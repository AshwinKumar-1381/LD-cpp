// compute.h

#ifndef COMPUTE_H
#define COMPUTE_H

#include "system.h"

#define rcut 1.122462048 

using namespace program;

namespace program{

	void computeNonBondedInteractions(atom_style *ATOMS, SimBox *BOX, pair_style *INTERACTION);
	void computeKineticEnergy(atom_style *ATOMS, SimBox *BOX);
	void computeTemperature(SimBox *BOX);
}

#endif /* COMPUTE_H */