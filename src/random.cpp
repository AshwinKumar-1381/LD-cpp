// random.cpp

#include "library.h"
#include "random.h"

void program::URN(float *zahl1, long *idum, int num)
{
	if(num > 0)
	{
		time_t t = time(&t);
		*idum = -(t%1000000) + 22067*num;
	}

	*idum ^= MASK;
	long k = (*idum)/IQ;
	*idum = IA*(*idum-k*IQ)-IR*k;
	if(*idum < 0) *idum += IM;

	*zahl1 = AM*(*idum);
	*idum ^= MASK;
}

void::program::GAUSS(float *gauss1, float *gauss2, long *idum, int num)
{
	if(num > 0)
	{
		time_t t = time(&t);
		*idum = -(t%1000000) + 22607*num;
	}

	float v1, v2, rsq;

	do {
		*idum ^= MASK;
		long k = (*idum)/IQ;
		*idum = IA*(*idum-k*IQ)-IR*k;
		if(*idum < 0) *idum += IM;
		float zahl1 = AM*(*idum);
		*idum ^= MASK;

		*idum ^= MASK;
		k = (*idum)/IQ;
		*idum = IA*(*idum-k*IQ)-IR*k;
		if(*idum < 0) *idum += IM;
		float zahl2 = AM*(*idum);
		*idum ^= MASK;

		v1 = 2.0*zahl1 - 1.0;
		v2 = 2.0*zahl2 - 1.0;
		rsq = v1*v1 + v2*v2;
	} while (rsq >= 1.0 || rsq == 0.0);

	float fac = sqrt(-2.0*log(rsq)/rsq);
	*gauss2 = v1*fac;
	*gauss1 = v2*fac;
}