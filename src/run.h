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

		runNVE(int id = 1, float t = 1.0, float delta_t = 1.0, int thermo_val = 1, int traj_val = 1, bool norm_val = true);
		~runNVE();

		void integrateNVE(atom_style *ATOMS, SimBox *BOX, pair_style *INTERACTION, sysInput *Input); 
	};

	class runLangevin{

		public:

		int runID;
		float time;
		float dt;
		int thermo_every;
		int traj_every;
		bool norm;
		bool zero;

		int maxSteps;

		runLangevin(int id = 1, float t = 1.0, float delta_t = 1.0, int thermo_val = 1, int traj_val = 1, bool norm_val = true, bool zero_val = true);
		~runLangevin();

		void integrateLangevin(atom_style *ATOMS, SimBox *BOX, pair_style *INTERACTION, sysInput *Input);
	};
}

#endif /*RUN_H*/ 