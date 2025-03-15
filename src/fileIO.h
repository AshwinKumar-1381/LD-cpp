// fileIO.h

#ifndef FILEIO_H
#define FILEIO_H

#include "system.h"
#include "kmc.h"
#include "run.h"

using namespace program;

namespace program{

	void readConfigFile(atom_style *ATOMS, SimBox *BOX, char *fname);
	void makeFolder(sysInput *Input);
	void writeLog(sysInput *Input, SimBox *BOX, run_style *RUN, KMC_poisson* KMC);
	void write2xyz(atom_style *ATOMS, sysInput *Input, float step, char *fname);
	void writeFrame(atom_style *ATOMS, sysInput *Input, char *fname);
	void writeThermo(SimBox *BOX, sysInput *Input, int runID, int fac, int step);
	void write2traj(atom_style *ATOMS, sysInput *Input, int runID, int step);
	void writeKMC(KMC_poisson *KMC, sysInput *Input, int step = -1);
}

#endif /* FILEIO_H */