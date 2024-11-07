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
		int thermo_every;
		int traj_every;
		bool norm;

		int maxSteps;

		runNVE(int id = 1, float t = 1.0, float delta_t = 1.0, int every1 = 1, int every2 = 1, bool val = true);
		~runNVE();

		void integrateNVE(atom_style *ATOMS, SimBox *BOX, pair_style *INTERACTION, sysInput *Input); 
	};

	class runLangevin{
	};
}

#endif /*RUN_H*/ 