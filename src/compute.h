// compute.h

#ifndef COMPUTE_H
#define COMPUTE_H

#include "system.h"

using namespace program;

namespace program{

	void computeNonBondedInteractions(atom_style *ATOMS, SimBox *BOX, interactions ***INTERACTION);
	void computeKineticEnergy(atom_style *ATOMS, SimBox *BOX);
	void computeTemperature(SimBox *BOX);
}

#endif /* COMPUTE_H */