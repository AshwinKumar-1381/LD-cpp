// random.h

#define pi 3.141592654

// Random-Number-Generators
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

using namespace program;

namespace program {
    
    void gauss();
    void URN();
}
