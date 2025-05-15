// fileIO.cpp

#include "library.h"
#include "fileIO.h"

using namespace program;

void program::readConfigFile(atom_style *ATOMS, SimBox *BOX, sysInput *Input, char *fname)
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
        int nAtoms;

        fgets(pipeString, 500, initFile);
        sscanf(pipeString, "%d %*s", &nAtoms);

        if(nAtoms == BOX->nAtoms)
        {
            fgets(pipeString, 500, initFile);
            sscanf(pipeString, "%*s %d", &Input->res_step);
        
            for(int i = 0; i < BOX->nAtoms; i++)
            {
                fgets(pipeString, 500, initFile);
                sscanf(pipeString, "%c %f %f %f %f %f", &ATOMS[i].id, &ATOMS[i].rx, &ATOMS[i].ry, &ATOMS[i].rz, &ATOMS[i].vx, &ATOMS[i].vy);
            }  

            fclose(initFile);
            delete(pipeString);

            printf("Read configurations for %d atoms from %s.\n", BOX->nAtoms, fname);  
        }
        else
        {
            printf("Error: Atom count mismatch between initial config. file and input script. Exiting...\n");
            exit(-1);
        }
    }
}

void program::makeFolder(sysInput *Input)
{
    if(Input->restart == false)
    {
        char *cmd = new char[50];
        char *fname = new char[50];
        sprintf(fname, "../%s/Data%d", Input->folder, Input->nr);

        sprintf(cmd, "rm -r %s", fname);
        std::system(cmd);
        sprintf(cmd, "mkdir -p %s/frames", fname);
        std::system(cmd);
        delete(cmd);
        delete(fname);
    }
}

void program::writeLog(sysInput *Input, SimBox *BOX, run_style *RUN, KMC_poisson* KMC)
{
    char *fname = new char[50];
    sprintf(fname, "../%s/Data%d/log%s.dat", Input->folder, Input->nr, RUN->name);
    remove(fname);

    FILE *log = fopen(fname, "w");
    if(log == NULL)
    {
        printf("Could not create Log file. Exiting...\n");
        exit(-1);
    }
    else
    {   
        time_t startTime = high_resolution_clock::to_time_t(Input->begin);
        time_t endTime = high_resolution_clock::to_time_t(Input->end);

        fprintf(log, "RUN PARAMS\n\n");
        fprintf(log, "Start Time : %s", ctime(&startTime));
        fprintf(log, "End Time : %s", ctime(&endTime));

        fprintf(log, "%s\n", program::returnElapsedTime(Input));
        fprintf(log, "nr = %d\n", Input->nr);
        fprintf(log, "Name : %s\n", RUN->name);
        fprintf(log, "Total time = %f\n", RUN->time);
        fprintf(log, "dt = %f\n", RUN->dt);
        fprintf(log, "Nsteps = %d\n\n", RUN->maxSteps - RUN->res_step);

        fprintf(log, "SIM PARAMS\n\n");
        fprintf(log, "Lx = %f\n", BOX->boxLength_x);
        fprintf(log, "Ly = %f\n", BOX->boxLength_y);
        fprintf(log, "pfrac = %f\n", Input->pfrac);
        fprintf(log, "N = %d\n", Input->N);
        fprintf(log, "PR = %f\n", Input->PR);
        fprintf(log, "m_str = %f\n", Input->m_str);
        fprintf(log, "D_str = %f\n", Input->D_str);
        fprintf(log, "Pe_A = %f\n", Input->PeA);
        fprintf(log, "Pe_B = %f\n\n", Input->PeB);

        if(KMC != NULL)
        {
            fprintf(log, "KMC PARAMS\n\n");
            fprintf(log, "lambda = %f\n", KMC->kmc_rate);
            fprintf(log, "bias param = %f\n", KMC->bias);
        }
    }

    fclose(log);
    delete(fname);
}

void program::write2xyz(atom_style *ATOMS, sysInput *Input, long step, char *fname)
{
    remove(fname);
    FILE *traj = fopen(fname, "w");
    if(traj == NULL)
    {
        printf("Cannot create initial trajectory file. Exiting...\n");
        exit(-1);
    }

    fprintf(traj, "%d\n", Input->N);
    fprintf(traj, " Atoms. Timestep: %ld\n", step);

    for(int i = 0; i < Input->N; i++) 
        fprintf(traj, "%c %f %f %f\n", ATOMS[i].id, ATOMS[i].rx, ATOMS[i].ry, ATOMS[i].rz);

    fclose(traj);
}

void program::writeFrame(atom_style *ATOMS, sysInput *Input, long step, char *fname)
{
    remove(fname);
    FILE *frame = fopen(fname, "a+");
    if(frame == NULL)
    {
        printf("Cannot create initial trajectory file. Exiting...\n");
        exit(-1);
    }

    fprintf(frame, "%d Atoms\n", Input->N);
    fprintf(frame, "Timestep: %ld\n", step);
    fprintf(frame, "ID rx ry rz vx vy\n");

    for(int i = 0; i < Input->N; i++) 
        fprintf(frame, "%c %f %f %f %f %f\n", ATOMS[i].id, ATOMS[i].rx, ATOMS[i].ry, ATOMS[i].rz, ATOMS[i].vx, ATOMS[i].vy);

    fclose(frame);
}

void program::writeThermo(SimBox *BOX, sysInput *Input, int runID, int fac, int step)
{
    char *fname = new char[500];
    FILE *thermo;
    sprintf(fname, "../%s/Data%d/thermo%d.dat", Input->folder, Input->nr, runID);

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

void program::write2traj(atom_style *ATOMS, sysInput *Input, int runID, int step)
{
    char *fname = new char[500];
    FILE *traj;
    sprintf(fname, "../%s/Data%d/traj%d.xyz", Input->folder, Input->nr, runID);

    if(step == 0)
    {
        remove(fname);
        traj = fopen(fname, "w");
        if(traj == NULL)
        {
            printf("Cannot create %s. Exiting...\n", fname);
            exit(-1);
        }
    }

    else
    {
        traj = fopen(fname, "a+");
        if(traj == NULL)
        {
            printf("Cannot open %s. File does not exist.\n", fname);
            exit(-1);
        }
    }

    fprintf(traj, "%d\n", Input->N);
    fprintf(traj, " Atoms. Timestep: %d\n", step);

    for(int i = 0; i < Input -> N; i++) 
        fprintf(traj, "%c %f %f %f %f %f\n", ATOMS[i].id, ATOMS[i].rx, ATOMS[i].ry, ATOMS[i].rz, ATOMS[i].vx, ATOMS[i].vy);

    fclose(traj);
}

void program::writeKMC(KMC_poisson *KMC, sysInput *Input, int step)
{
    char *fname = new char[500];
    FILE *file;
    sprintf(fname, "../%s/Data%d/kmc.dat", Input->folder, Input->nr);

    if(step == 0)
    {
        char *fname2 = new char[500];
        sprintf(fname2, "../%s/Data%d/kmc_switch_times.dat", Input->folder, Input->nr);
        remove(fname2);

        file = fopen(fname2, "w");
        if(file == NULL)
        {
            printf("Cannot create %s file. Exiting...\n", fname2);
            exit(-1);
        }
        else
        {
            fprintf(file, "pid tau\n");

            nodeKMC *curr = KMC->head;
            while(curr->next != NULL)
            {
                curr = curr->next;
                fprintf(file, "%d %d\n", curr->pid, curr->tauStep);
            }
        }

        fclose(file);

        remove(fname);
    
        file = fopen(fname, "w");
        if(file == NULL)
        {
            printf("Cannot create %s file. Exiting...\n", fname);
            exit(-1);
        }
        else
            fprintf(file, "Step numA numB\n");
    }

    else
    {
        file = fopen(fname, "a+");
        if(file == NULL)
        {
            printf("Cannot open %s file. Exiting...\n", fname);
            exit(-1);
        }
    }

    fprintf(file, "%d %d %d\n", step, KMC->numA, KMC->numB);
    fclose(file);
}

char* program::returnElapsedTime(sysInput *Input)
{
    auto dur = duration_cast<seconds>(Input->end - Input->begin);

    int t_time = dur.count();
    int time[3] = {0, 0, 0}; 
    int fac1 = 3600;
    for (int i = 0; i<3; i++)
    {
        time[i] = int(t_time/fac1);
        t_time -= time[i]*fac1;
        fac1 /= 60;    
    }
    
    char *strtime = new char[50];
    sprintf(strtime, "Total simulation time = %d hrs %d mins %d secs\n", time[0], time[1], time[2]);

    return(strtime);
}