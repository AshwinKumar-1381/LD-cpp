// run.h

#ifndef RUN_H
#define RUN_H

#include "system.h"
#include "kmc.h"

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

		void integrateNVE(atom_style *ATOMS, SimBox *BOX, interactions ***INTERACTIONS, sysInput *Input); 
	};

	typedef class runLangevin_D{

		public:

		int runID;
		float time, dt;
		int thermo_every, traj_every, maxSteps;
		bool norm, zero, kmc;
		float field_loc_x; 

		runLangevin_D(int id = 1, float t = 1.0, float delta_t = 1.0, int thermo_val = 1, int traj_val = 1, bool norm_val = true, bool zero_val = true, bool kmc_val = false, float field_loc_x = 0.0);
		~runLangevin_D();

		void integrateLangevin(atom_style *ATOMS, SimBox *BOX, interactions ***INTERACTIONS, sysInput *Input, KMC_poisson *KMC = NULL);
	}run_style;

	typedef class runBrownian{

		public:

		int runID;
		float time, dt;
		int thermo_every, traj_every, maxSteps;
		bool norm, zero, kmc;
		float field_loc_x; 

		runBrownian(int id = 1, float t = 1.0, float delta_t = 1.0, int thermo_val = 1, int traj_val = 1, bool norm_val = true, bool zero_val = true, bool kmc_val = false, float field_loc_x = 0.0);
		~runBrownian();

		void integrateBrownian(atom_style *ATOMS, SimBox *BOX, interactions ***INTERACTIONS, sysInput *Input, KMC_poisson *KMC = NULL);
	};
}

#endif /*RUN_H*/ 