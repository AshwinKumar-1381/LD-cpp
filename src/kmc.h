// kmc.h

#ifndef KMC_H
#define KMC_H

#include "system.h"

using namespace program;

namespace program{

	struct nodeKMC{

		int pid;
		int tau;
		int tauStep;
		int numSwitches;
		nodeKMC *prev;
		nodeKMC *next;

		nodeKMC(int pid, int tau, int tauStep);
		~nodeKMC(){}
	};

	struct dispDist{

		float dr, dtau;				// bin widths
		float rmax, taumax;			// max bin vals
		int nTau, nR;				// num bins
		float **Dist;				// 2 - D matrix of size (nTau x nR)
		long nSamples;				// counts total no.of samples
		int delay;					// calculate distribution after this step

		char *fpathO;

		dispDist(float dr, float rmax, float dtau, float taumax, int delay);
		~dispDist(){}

		void normalize(float dt);
		void write2file();
	};

	class KMC_poisson{

		public:
		
		float kmc_rate;
		float bias;
		int delay;
		bool verbose;
		bool dist;

		long idum;
		float zahl1;
		int numA, numB;

		nodeKMC *head; 
		dispDist *dist_msd_tau;

		KMC_poisson(float rate_val = 1.0, float bias_val = 1.0, int delay_val = 0, bool v_val = false, bool dist_val = false);
		~KMC_poisson();

		void initialize(sysInput *Input, float dt, atom_style *ATOMS = NULL, SimBox *BOX = NULL);
		int Sample(float dt);
		void insertNode(nodeKMC *newNode);
		void logToDist(atom_style *ATOMS, SimBox *BOX, nodeKMC *node);
		void Switch(atom_style *ATOMS, SimBox *BOX, sysInput *Input, float dt, int step); 
	};
}

#endif /*KMC_H*/