#include <mpi.h>

#include "initLB.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *iproc, int *jproc,
        int *kproc, int *timesteps, int *timestepsPerPlotting, int argc, char **argv)
{
    /* Checking if there are exactly 2 parameters in the input. Else give Error message. */
    if (argc != 2) {
        printf("You must provide a filename!\n");
        return -1;
    }
    const char *filename = argv[1];
    /* Reading all the input parameters */
    READ_INT(filename, *xlength);
    READ_DOUBLE(filename, *tau);
    read_double(filename, "vwallx", &velocityWall[0]);
    read_double(filename, "vwally", &velocityWall[1]);
    read_double(filename, "vwallz", &velocityWall[2]);
    READ_INT(filename, *iproc);
    READ_INT(filename, *jproc);
    READ_INT(filename, *kproc);
    READ_INT(filename, *timesteps);
    READ_INT(filename, *timestepsPerPlotting);

    return 0;
}

void broadcastInitialValues(int *xlength, double *tau, double *velocityWall, int *iproc, int *jproc,
        int *kproc, int *timesteps)
{
    MPI_Bcast(xlength, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /* Note: 3 doubles for velocityWall */
    MPI_Bcast(velocityWall, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(iproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(jproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(kproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(timesteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

#if 0 /* Maybe remove this later on */
typedef enum
{
    NONE = 0x0, LEFT = 0x1, RIGHT = 0x2, TOP = 0x4, BOTTOM = 0x8, FRONT = 0x10, BACK = 0x20
}BOUNDARY;

int checkBoundary(int rank, int xlength)
{
    int result = NONE;
    if (rank / (xlength * xlength) == 0)
    result |= BOTTOM;
    if (rank / (xlength * xlength) == xlength - 1)
    result |= TOP;
    if ((rank / xlength) % xlength == 0)
    result |= FRONT;
    if ((rank / xlength) % xlength == xlength - 1)
    result |= BACK;
    if (rank % xlength == 0)
    result |= LEFT;
    if (rank % xlength == xlength - 1)
    result |= RIGHT;
    return result;
}
#endif

void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength,
        const int * const sublength, int rank, int number_of_ranks)
{
    /* TODO: Change the memory allocation according to the lattice sizes for each MPI process */
    const int xl = sublength[0];
    const int yl = sublength[1];
    const int zl = sublength[2];

    /* Set all values of flagField to 0, i.e. FLUID state. We will apply boundary conditions
     in the following */
    memset(flagField, FLUID, (xl + 2) * (yl + 2) * (zl + 2) * sizeof(*flagField));

    /* TODO: Check if this is working correctly and improve performance */
    /* The values for Boundary on Z = 0 are set to No_Slip, whereas Z=Zmax is set to MOVING_WALL */
    for (int x = 0; x < xl + 2; ++x) {
        for (int y = 0; y < yl + 2; ++y) {
            /* Subdomains in bottom part of domain */
            if (rank / (xlength * xlength) == 0) {
                flagField[fidx(sublength, x, y, 0)] = NO_SLIP;
                flagField[fidx(sublength, x, y, zl + 1)] = PARALLEL_BOUNDARY;
            }
            else if (rank / (xlength * xlength) == xlength - 1) {
                /* Subdomains in upper part of domain - Top is the moving wall, bottom is a no-slip boundary */
                flagField[fidx(sublength, x, y, 0)] = PARALLEL_BOUNDARY;
                flagField[fidx(sublength, x, y, zl + 1)] = MOVING_WALL;
            }
            else { /* Inner cells have parallel boundaries */
                flagField[fidx(sublength, x, y, 0)] = PARALLEL_BOUNDARY;
                flagField[fidx(sublength, x, y, zl + 1)] = PARALLEL_BOUNDARY;
            }
        }
    }

    /* The values for Boundary on Y = 0 and Ymax plane set to No_Slip */
    for (int x = 0; x < xl + 2; ++x) {
        for (int z = 0; z < zl + 2; ++z) {
            if ((rank / xlength) % xlength == 0) {
                flagField[fidx(sublength, x, 0, z)] = NO_SLIP;
                flagField[fidx(sublength, x, yl + 1, z)] = PARALLEL_BOUNDARY;
            }
            else if ((rank / xlength) % xlength == xlength - 1) {
                flagField[fidx(sublength, x, 0, z)] = PARALLEL_BOUNDARY;
                flagField[fidx(sublength, x, yl + 1, z)] = NO_SLIP;
            }
            else {
                flagField[fidx(sublength, x, 0, z)] = PARALLEL_BOUNDARY;
                flagField[fidx(sublength, x, yl + 1, z)] = PARALLEL_BOUNDARY;
            }
        }
    }
    /* The values for Boundary on Y = 0 and Ymax plane set to No_Slip */
    for (int y = 0; y < yl + 2; ++y) {
        for (int z = 0; z < zl + 2; ++z) {
            if (rank % xlength == 0) {
                flagField[fidx(sublength, 0, y, z)] = NO_SLIP;
                flagField[fidx(sublength, xl + 1, y, z)] = PARALLEL_BOUNDARY;
            }
            else if (rank % xlength == xlength - 1) {
                flagField[fidx(sublength, 0, y, z)] = PARALLEL_BOUNDARY;
                flagField[fidx(sublength, xl + 1, y, z)] = NO_SLIP;
            }
            else {
                flagField[fidx(sublength, 0, y, z)] = PARALLEL_BOUNDARY;
                flagField[fidx(sublength, xl + 1, y, z)] = PARALLEL_BOUNDARY;
            }
        }
    }

    /* Stream and Collide Fields are initialized to the respective lattice weights of the Cell */
    for (int x = 0; x < xl + 2; ++x) {
        for (int y = 0; y < yl + 2; ++y) {
            for (int z = 0; z < zl + 2; ++z) {
                for (int i = 0; i < Q; ++i) {
                    collideField[idx(sublength, x, y, z, i)] = LATTICEWEIGHTS[i];
                    streamField[idx(sublength, x, y, z, i)] = LATTICEWEIGHTS[i];
                }
            }
        }
    }
}

