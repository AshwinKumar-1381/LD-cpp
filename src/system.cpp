// system.cpp

#include "library.h" 
#include "system.h"

using namespace program;

program::SimBox::SimBox(){
    nAtoms = 1;
    boxLength = 1.0;
    numFrac = 1.0;
    ke = 0.0;
    pe = 0.0;
}

void program::SimBox::initBox(Atoms2D *ATOMS, float pfrac, float L, int N, char *fname=NULL)
{    
    nAtoms = N;
    boxLength = L;
    numFrac = pfrac;
    
    if(fname != NULL) readConfigFile(ATOMS, fname);
}

void program::SimBox::readConfigFile(Atoms2D *ATOMS, char *fname){

    printf("%s\n",fname);
    
    char *pipeString = (char*) malloc(500*sizeof(char));
    FILE *initFile = fopen(fname, "r");
    
    if(initFile == NULL) 
    {
        printf("Can not open %s . File does not exist.\n", fname);
        exit(-1);
    }
    while(!feof(initFile))
    {
        fgets(pipeString, 500, initFile);
        fgets(pipeString, 500, initFile);
        
        for(int i = 0; i < nAtoms; i++)
        {
            fgets(pipeString, 500, initFile);
            //printf("Atom %d - %c %f %f \n", i+1, ATOMS[i].id, ATOMS[i].rx, ATOMS[i].ry);
            sscanf(pipeString, "%c %f %f %*f", &(ATOMS[i].id), &(ATOMS[i].rx), &(ATOMS[i].ry));
            //printf("Atom %d - %c %f %f \n", i+1, ATOMS[i].id, ATOMS[i].rx, ATOMS[i].ry);
        }
        
        //for(int i = 0; i < nAtoms; i++) printf("Atom %d - %c %f %f \n", i+1, ATOMS[i].id, ATOMS[i].rx, ATOMS[i].ry);
    
    //exit(-1);    
    }
    
    fclose(initFile);
    
    //for (int i = 0; i < nAtoms; i++) printf("Atom %d - %c %f %f \n", i+1, ATOMS[i].id, ATOMS[i].rx, ATOMS[i].ry);
    
}

