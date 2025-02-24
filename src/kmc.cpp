// kmc.cpp

#include "library.h"
#include "kmc.h"
#include "random.h"

using namespace program;

program::nodeKMC::nodeKMC(int pid, int tau, int tauStep)
{
	this->pid = pid;
	this->tau = tau;
	this->tauStep = tauStep;
	numSwitches = 0;
	prev = next = nullptr;
}

program::dispDist::dispDist(float dr, float rmax, float dtau, float taumax, int delay)
{
	this->dr = dr;
	this->rmax = rmax;
	this->dtau = dtau;
	this->taumax = taumax;
	this->delay = delay;

	nR = ceil(rmax/dr);
	nTau = ceil(taumax/dtau);

	Dist = new float*[nTau];
	for(int i = 0; i < nTau; i++)
	{
		Dist[i] = new float[nR];
		for(int j = 0; j < nR; j++)
			Dist[i][j] = 0.0;
	}

	nSamples = 0;
	fpathO = new char[50];
}

void program::dispDist::normalize(float dt)
{
	for(int i = 0; i < nTau; i++)
	{
		for(int j = 0; j < nR; j++)
			Dist[i][j] /= (nSamples * dr * dtau * dt);
	}
}

void program::dispDist::write2file()
{
	FILE *fileO = fopen(fpathO, "w");

	if(fileO == NULL)
	{
		printf("Cannot create KMC distribution data file %s.\n", fpathO);
		exit(-1);
	}
	else
	{
		fprintf(fileO, "delR prob\n");
		for(int i = 0; i < nTau; i++)
		{
			fprintf(fileO, "delTau %f\n", (i+0.5)*dtau);
			for(int j = 0; j < nR; j++)
				fprintf(fileO, "%f %f\n", (j+0.5)*dr, Dist[i][j]);
		}
	}

	fclose(fileO);
}

program::KMC_poisson::KMC_poisson(float rate_val, float bias_val, int delay_val, bool v_val, bool dist_val)
{
	kmc_rate = rate_val;
	bias = bias_val;
	delay = delay_val;
	verbose = v_val;
	dist = dist_val;
}

int program::KMC_poisson::Sample(float dt)
{
	float tau;

	do{
		program::URN(&zahl1, &idum);
		tau = -1.0*log(zahl1);
		tau = ceil(tau/(kmc_rate * dt));

	} while(int(tau) < 1);
	
	return(int(tau));
}

void program::KMC_poisson::initialize(sysInput *Input, float dt, atom_style *ATOMS, SimBox *BOX)
{
	numB = Input->PR*Input->N;
	numA = Input->N - numB;

	program::URN(&zahl1, &idum, 5);

	for(int i = 0; i < Input->N; i++)
	{
		int tau = Sample(dt);

		if(i == 0)
		{	
			head = new nodeKMC(-1, -1, -1);
			nodeKMC *newNode = new nodeKMC(i, tau, delay + tau);

			head -> prev = nullptr;
			head -> next = newNode;
			newNode -> prev = head;
			newNode -> next = nullptr;
		}

		else
		{
			nodeKMC *newNode = new nodeKMC(i, tau, delay + tau);
			insertNode(newNode);
		}

		if(verbose == true)
			printf("\nHi from Node %d with tau %d\n", i, tau);
	}

	if(dist == true)
	{
		for(int i = 0; i < Input->N; i++)
		{
			float fac = 0.0;
			if(ATOMS[i].id == 'N') fac = 1.0;

			ATOMS[i].xsave = ATOMS[i].rx * fac;
			ATOMS[i].ysave = ATOMS[i].ry * fac;
		}

		float rmax = max(BOX->boxLength_x, BOX->boxLength_y);
		float tauAvg = ceil(1.0/(kmc_rate*dt)); 		// Avg. no.of steps required to switch
		dist_msd_tau = new dispDist(4.0, rmax, 0.1*tauAvg, 5.0*tauAvg, tauAvg);
	}
}

void program::KMC_poisson::logToDist(atom_style *ATOMS, SimBox *BOX, nodeKMC *node)
{
	int pid = node -> pid;
	float fac;

	if(ATOMS[pid].id == 'N')
	{
		if(node->tauStep >= dist_msd_tau->delay)
		{
			float dx = ATOMS[pid].rx + BOX->boxLength_x*ATOMS[pid].jumpx - ATOMS[pid].xsave;
			float dy = ATOMS[pid].ry + BOX->boxLength_y*ATOMS[pid].jumpy - ATOMS[pid].ysave;

			float delr = sqrt(dx*dx + dy*dy);
			float tau = node -> tau;

			if(delr <= dist_msd_tau->rmax and tau <= dist_msd_tau->taumax)
			{
				int i = int(tau/dist_msd_tau->dtau);
				int j = int(delr/dist_msd_tau->dr);

				dist_msd_tau -> Dist[i][j] += 1.0;
				dist_msd_tau -> nSamples += 1; 
			}	
		}
		fac = 0.0;
	}

	else if(ATOMS[pid].id == 'O') fac = 1.0;

	ATOMS[pid].xsave = ATOMS[pid].rx*fac;
	ATOMS[pid].ysave = ATOMS[pid].ry*fac;
}

void program::KMC_poisson::Switch(atom_style *ATOMS, SimBox *BOX, sysInput *Input, float dt, int step)
{
	if(step == (head->next)->tauStep)
	{
		nodeKMC *curr = head -> next;
		while(curr -> tauStep == step)
		{
			program::URN(&zahl1, &idum);

			if(ATOMS[curr->pid].id == 'N')
			{
				if(zahl1 <= bias/(1+bias))
				{
					if(dist == true) logToDist(ATOMS, BOX, curr);

					ATOMS[curr->pid].id = 'O';
					ATOMS[curr->pid].Pe = Input->PeB;
					numA--;
					numB++;	
				}
			}

			else if(ATOMS[curr->pid].id == 'O')
			{
				if(zahl1 <= 1.0/(1+bias))
				{
					if(dist == true) logToDist(ATOMS, BOX, curr);

					ATOMS[curr->pid].id = 'N';
					ATOMS[curr->pid].Pe = Input->PeA;
					numB--;
					numA++;	
				}
			}

			curr -> numSwitches++;
			curr -> tau = Sample(dt);
			curr -> tauStep = step + curr -> tau;

			if(verbose == true)
				printf("\nSwitching Node %d at step %d\nNext switch scheduled at step %d\n", curr->pid, step, curr->tauStep);

			curr = curr -> next;
			curr -> prev = head;
			nodeKMC *detachNode = head -> next;
			head -> next = curr;
			detachNode -> prev = nullptr;
			detachNode -> next = nullptr;

			insertNode(detachNode);
		}
	}
}

void program::KMC_poisson::insertNode(nodeKMC *newNode)
{
	nodeKMC *curr = head;
	int val = newNode -> tauStep;

	int accept = 0;
	while(curr -> next != nullptr)
	{
		curr = curr -> next;
		if(curr -> tauStep >= val)
		{
			if((curr->prev) -> tauStep == -1)
			{
				if(verbose==true) 
					printf("Inserting at first before %d\n", curr->tauStep);
				
				newNode -> prev = head;
				newNode -> next = curr;
				curr -> prev = newNode;
				head -> next = newNode;
			}

			else
			{
				if(verbose==true) 
					printf("Inserting in between %d and %d\n", (curr->prev)->tauStep, curr->tauStep);
				
				newNode -> next = curr;
				curr = curr -> prev;
				newNode -> prev = curr;
				curr -> next = newNode;
				curr = newNode -> next;
				curr -> prev = newNode;
			}

			accept = 1;
		}

		if(accept == 1) return;
	}

	if(accept == 0)
	{
		if(verbose==true) 
			printf("Inserting at last after %d\n", curr->tauStep);
		
		newNode -> prev = curr;
		newNode -> next = nullptr;
		curr -> next = newNode;
		return;
	}
}