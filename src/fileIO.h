// fileIO.h

#ifndef FILEIO_H
#define FILEIO_H

#include "system.h"

using namespace program;

namespace program{

	void readConfigFile(atom_style *ATOMS, SimBox *BOX, char *fname);

	void makeFolder(sysInput *Input);
	void writeLog(sysInput *Input);
	void write2xyz(atom_style *ATOMS, sysInput *Input, float step, char *fname);
	void writeFrame(atom_style *ATOMS, sysInput *Input, char *fname);
	void writeThermo(SimBox *BOX, sysInput *Input, int runID, int fac, int step);
	void write2traj(atom_style *ATOMS, sysInput *Input, int runID, int step);
}

#endif /* FILEIO_H */