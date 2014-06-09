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
        int *kproc, int *timesteps, int *timestepsPerPlotting)
{
    MPI_Bcast(xlength, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /* Note: 3 doubles for velocityWall */
    MPI_Bcast(velocityWall, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(iproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(jproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(kproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(timesteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(timestepsPerPlotting, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void initialiseFields(double *collideField, double *streamField, int *flagField,
        const int * const sublength, int rank, int iproc, int jproc, int kproc)
{
    /* TODO: Change the memory allocation according to the lattice sizes for each MPI process */
    const int xl = sublength[0];
    const int yl = sublength[1];
    const int zl = sublength[2];
    /* Minor value and major value to set for a plane (i.e. FRONT and BACK or BOTTOM and TOP) */
    STATE minor, major;

    /* Set all values of flagField to 0, i.e. FLUID state. We will apply boundary conditions
     in the following */
    memset(flagField, FLUID, (xl + 2) * (yl + 2) * (zl + 2) * sizeof(*flagField));

    /* TODO: Check if this is working correctly and improve performance */

    /* XY plane */
    if (checkBoundary(rank, iproc, jproc, kproc) & BOTTOM) {
        printf("Rank #%d - BOTTOM\n", rank);
        /* Subdomains in bottom part of domain */
        minor = NO_SLIP;
        major = PARALLEL_BOUNDARY;
    }
    if (checkBoundary(rank, iproc, jproc, kproc) & TOP) {
        printf("Rank #%d - TOP\n", rank);
        /* Subdomains in upper part of domain - Top is the moving wall, bottom is a no-slip boundary */
        minor = PARALLEL_BOUNDARY;
        major = MOVING_WALL;
    }
    else {
        /* Inner cells have parallel boundaries */
        minor = PARALLEL_BOUNDARY;
        major = PARALLEL_BOUNDARY;
    }
    for (int x = 0; x < xl + 2; ++x) {
        for (int y = 0; y < yl + 2; ++y) {
            flagField[fidx(sublength, x, y, 0)] = minor;
            flagField[fidx(sublength, x, y, zl + 1)] = major;
        }
    }

    /* XZ plane */
    if (checkBoundary(rank, iproc, jproc, kproc) & FRONT) {
        printf("Rank #%d - FRONT\n", rank);
        minor = NO_SLIP;
        major = PARALLEL_BOUNDARY;
    }
    else if (checkBoundary(rank, iproc, jproc, kproc) & BACK) {
        printf("Rank #%d - BACK\n", rank);
        minor = PARALLEL_BOUNDARY;
        major = NO_SLIP;
    }
    else {
        minor = PARALLEL_BOUNDARY;
        major = PARALLEL_BOUNDARY;
    }
    for (int x = 0; x < xl + 2; ++x) {
        for (int z = 0; z < zl + 2; ++z) {
            flagField[fidx(sublength, x, 0, z)] = minor;
            flagField[fidx(sublength, x, yl + 1, z)] = major;
        }
    }

    /* YZ plane */
    if (checkBoundary(rank, iproc, jproc, kproc) & LEFT) {
        printf("Rank #%d - LEFT\n", rank);
        minor = NO_SLIP;
        major = PARALLEL_BOUNDARY;
    }
    else if (checkBoundary(rank, iproc, jproc, kproc) & RIGHT) {
        printf("Rank #%d - RIGHT\n", rank);
        minor = PARALLEL_BOUNDARY;
        major = NO_SLIP;
    }
    else {
        minor = PARALLEL_BOUNDARY;
        major = PARALLEL_BOUNDARY;
    }
    /* The values for Boundary on Y = 0 and Ymax plane set to No_Slip */
    for (int y = 0; y < yl + 2; ++y) {
        for (int z = 0; z < zl + 2; ++z) {
            flagField[fidx(sublength, 0, y, z)] = minor;
            flagField[fidx(sublength, xl + 1, y, z)] = major;
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

