// fileIO.cpp

#include "library.h"
#include "fileIO.h"

using namespace program;

void program::readConfigFile(atom_style *ATOMS, SimBox *BOX, char *fname)
{
    char *pipeString = new char[500];
    FILE *initFile = fopen(fname, "r");
    
    if(initFile == NULL) 
    {
        printf("Can not open %s. File does not exist.\n", fname);
        exit(-1);
    }
    else
    {
        fgets(pipeString, 500, initFile);
        
        for(int i = 0; i < BOX->nAtoms; i++)
        {
            fgets(pipeString, 500, initFile);
            sscanf(pipeString, "%c %f %f %f %f", &ATOMS[i].id, &ATOMS[i].rx, &ATOMS[i].ry, &ATOMS[i].px, &ATOMS[i].py);
        }   
    }
    
    fclose(initFile);
    delete(pipeString);

    printf("Read configurations for %d atoms from %s.\n", BOX->nAtoms, fname);
}

void program::makeFolder(sysInput *Input)
{
    char *cmd = new char[50];
    char *fname = new char[50];
    sprintf(fname, "../Data%d", Input->nr);

    sprintf(cmd, "rm -r %s", fname);
    std::system(cmd);
    sprintf(cmd, "mkdir %s", fname);
    std::system(cmd);
    delete(cmd);
    delete(fname);
}

void program::writeLog(sysInput *Input)
{
    char *fname = new char[50];
    sprintf(fname, "../Data%d/log.dat", Input->nr);
    
    FILE *log = fopen(fname, "w");
    if(log == NULL)
    {
        printf("Could not create Log file. Exiting...\n");
        exit(-1);
    }
    else
    {
        fprintf(log, "nr %d\n", Input->nr);
        fprintf(log, "lj_ep %f\n", Input->nr);
        fprintf(log, "L %f\n", Input->L);
        fprintf(log, "pfrac %f\n", Input->pfrac);
        fprintf(log, "m_str %f\n", Input->m_str);
        fprintf(log, "PR %f\n", Input->PR);
        fprintf(log, "Pe_A %f\n", Input->PeA);
        fprintf(log, "Pe_B %f\n", Input->PeB);
        fprintf(log, "N %d\n", Input->N);
    }
    fclose(log);
    delete(fname);
}

void program::write2xyz(atom_style *ATOMS, sysInput *Input, float step, char *fname)
{
    remove(fname);
    FILE *traj = fopen(fname, "a+");
    if(traj == NULL)
    {
        printf("Cannot create initial trajectory file. Exiting...\n");
        exit(-1);
    }

    fprintf(traj, "%d\n", Input->N);
    fprintf(traj, "Timestep %f\n", step);

    for(int i = 0; i < Input->N; i++) 
        fprintf(traj, "%c %f %f %f\n", ATOMS[i].id, ATOMS[i].rx, ATOMS[i].ry);

    fclose(traj);
}

void program::writeFrame(atom_style *ATOMS, sysInput *Input, char *fname)
{
    remove(fname);
    FILE *frame = fopen(fname, "a+");
    if(frame == NULL)
    {
        printf("Cannot create initial trajectory file. Exiting...\n");
        exit(-1);
    }

    fprintf(frame, "ID rx ry px py\n");

    for(int i = 0; i < Input->N; i++) 
        fprintf(frame, "%c %f %f %f %f\n", ATOMS[i].id, ATOMS[i].rx, ATOMS[i].ry, ATOMS[i].px, ATOMS[i].py);

    fclose(frame);
}

void program::writeThermo(SimBox *BOX, sysInput *Input, int runID, int fac, int step)
{
    char *fname = new char[500];
    FILE *thermo;
    sprintf(fname, "../Data%d/thermo%d.dat", Input->nr, runID);

    if(step == 0)
    {
        remove(fname);
        thermo = fopen(fname, "w");
        if(thermo == NULL)
        {
            printf("Cannot create %s. Exiting...\n", fname);
            exit(-1);
        }
        else fprintf(thermo, "step pe ke etot temp\n");
    }

    else 
    {
        thermo = fopen(fname, "a+");
        if(thermo == NULL)
        {
            printf("Cannot open %s. File does not exist.\n", fname);
            exit(-1);
        }
    }

    fprintf(thermo, "%d %f %f %f %f\n", step, BOX->pe/fac, BOX->ke/fac, BOX->etot/fac, BOX->temp);
    fclose(thermo);
}   