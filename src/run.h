// run.h

#ifndef RUN_H
#define RUN_H

#include "system.h"

using namespace program;

namespace program{

	class runNVE{

		public:

		int runID;
		float time;
		float dt;
		int maxSteps;

		runNVE();
		~runNVE();
		void integrateNVE(atom_style *ATOMS, SimBox *BOX, pair_style *INTERACTION, sysInput *Input); 
	};
}

#endif /*RUN_H*/ 