// kmc.cpp

#include "library.h"
#include "kmc.h"
#include "random.h"

using namespace program;

program::KMC_poisson::KMC_poisson(float rate_val, float bias_val, int delay_val, bool v_val)
{
	kmc_rate = rate_val;
	bias = bias_val;
	delay = delay_val;
	verbose = v_val;
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

void program::KMC_poisson::initialize(sysInput *Input, float dt)
{
	numB = Input->PR*Input->N;
	numA = Input->N - numB;

	int *switchTimes = new int[Input->N];

	program::URN(&zahl1, &idum, 5);

	for(int i = 0; i < Input->N; i++) 
		switchTimes[i] = delay + Sample(dt); 

	for(int i = 0; i < Input->N; i++)
	{
		if(verbose==true)
			printf("\nHi from Node %d with tau %d\n", i, switchTimes[i]);

		if(i == 0)
		{	
			head = new nodeKMC(-1, -1);
			nodeKMC *newNode = new nodeKMC(i, switchTimes[i]);

			head -> prev = nullptr;
			head -> next = newNode;
			newNode -> prev = head;
			newNode -> next = nullptr;
		}

		else
		{
			nodeKMC *newNode = new nodeKMC(i, switchTimes[i]);
			insertNode(newNode);
		}
	}

	delete switchTimes;
}

void program::KMC_poisson::Switch(atom_style *ATOMS, sysInput *Input, float dt, int step)
{
	if(step == (head->next)->tauStep)
	{
		nodeKMC *curr = head -> next;
		while(curr -> tauStep == step)
		{
			program::URN(&zahl1, &idum);

			if(zahl1 <= bias/(1+bias))
			{
				if(ATOMS[curr->pid].id == 'N')
				{
					ATOMS[curr->pid].id = 'O';
					ATOMS[curr->pid].Pe = Input->PeB;
					numA--;
					numB++;
				}
			}
				
			else
			{
				if(ATOMS[curr->pid].id == 'O')
				{
					ATOMS[curr->pid].id = 'N';
					ATOMS[curr->pid].Pe = Input->PeA;
					numB--;
					numA++;
				}
			}

			curr -> numSwitches++;
			curr -> tauStep = step + Sample(dt);

			if(verbose == true)
				printf("\nSwitching Node %d at step %d\nNext switch at step %d\n", curr->pid, step, curr->tauStep);

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