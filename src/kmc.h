// kmc.h

#ifndef KMC_H
#define KMC_H

#include "system.h"

using namespace program;

namespace program{

	struct nodeKMC{

		int pid;
		int tauStep;
		int numSwitches;
		nodeKMC *prev;
		nodeKMC *next;

		nodeKMC(int pid, int tauStep)
		{
			this->pid = pid;
			this->tauStep = tauStep;
			numSwitches = 0;
			prev = next = nullptr;
		}	
	};

	class KMC_poisson{

		public:

		float kmc_rate;
		float bias;
		int delay;
		bool verbose;
		long idum;
		float zahl1;
		nodeKMC *head; 
		int numA, numB;

		KMC_poisson(float rate_val = 1.0, float bias_val = 1.0, int delay_val = 0, bool v_val = false);
		~KMC_poisson();

		void initialize(sysInput *Input, float dt);
		int Sample(float dt);
		void insertNode(nodeKMC *newNode);
		void Switch(atom_style *ATOMS, sysInput *Input, float dt, int step); 
	};
}

#endif /*KMC_H*/