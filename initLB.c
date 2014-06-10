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

void allocateBufferSpace(double **sendBuffer, double **readBuffer, const int * const sublength)
{
    const int xl = sublength[0];
    const int yl = sublength[1];
    const int zl = sublength[2];
    /* We allocate space for 5 PDFs in every plane direction times the size of the plane */
    /* Send buffers */
    sendBuffer[0] = malloc(5 * sizeof(double) * (yl + 2) * (zl + 2));
    sendBuffer[1] = malloc(5 * sizeof(double) * (yl + 2) * (zl + 2));
    sendBuffer[2] = malloc(5 * sizeof(double) * (xl + 2) * (yl + 2));
    sendBuffer[3] = malloc(5 * sizeof(double) * (xl + 2) * (yl + 2));
    sendBuffer[4] = malloc(5 * sizeof(double) * (xl + 2) * (zl + 2));
    sendBuffer[5] = malloc(5 * sizeof(double) * (xl + 2) * (zl + 2));
    /* Read buffers */
    readBuffer[0] = malloc(5 * sizeof(double) * (yl + 2) * (zl + 2));
    readBuffer[1] = malloc(5 * sizeof(double) * (yl + 2) * (zl + 2));
    readBuffer[2] = malloc(5 * sizeof(double) * (xl + 2) * (yl + 2));
    readBuffer[3] = malloc(5 * sizeof(double) * (xl + 2) * (yl + 2));
    readBuffer[4] = malloc(5 * sizeof(double) * (xl + 2) * (zl + 2));
    readBuffer[5] = malloc(5 * sizeof(double) * (xl + 2) * (zl + 2));
}

void cleanup(double *collideField, double *streamField, int *flagField, double **sendBuffer,
        double **readBuffer)
{
    free(collideField);
    free(streamField);
    free(flagField);
    for (int i = 0; i < 6; ++i) {
        free(sendBuffer[i]);
        free(readBuffer[i]);
    }
}

void initialiseFields(double *collideField, double *streamField, int *flagField,
        const int * const sublength, int rank, int iproc, int jproc, int kproc)
{
    const int xl = sublength[0];
    const int yl = sublength[1];
    const int zl = sublength[2];

    /* Minor value and major value to set for a plane (i.e. FRONT and BACK or BOTTOM and TOP) */
    STATE minor, major;

    /* Set all values of flagField to 0, i.e. FLUID state. We will apply boundary conditions
     in the following */
    memset(flagField, FLUID, (xl + 2) * (yl + 2) * (zl + 2) * sizeof(*flagField));

    /* XY plane */
    if (checkBoundary(rank, iproc, jproc, kproc) & BOTTOM) {
        /* Subdomains in bottom part of domain */
        minor = NO_SLIP;
        major = PARALLEL_BOUNDARY;
    }
    if (checkBoundary(rank, iproc, jproc, kproc) & TOP) {
        minor = PARALLEL_BOUNDARY;
        major = MOVING_WALL;
    }
    else {
        /* Inner cells have parallel boundaries */
        minor = PARALLEL_BOUNDARY;
        major = PARALLEL_BOUNDARY;
    }
    for (int y = 0; y < yl + 2; ++y) {
        for (int x = 0; x < xl + 2; ++x) {
            flagField[fidx(sublength, x, y, 0)] = minor;
            flagField[fidx(sublength, x, y, zl + 1)] = major;
        }
    }

    /* XZ plane */
    if (checkBoundary(rank, iproc, jproc, kproc) & FRONT) {
        minor = NO_SLIP;
        major = PARALLEL_BOUNDARY;
    }
    else if (checkBoundary(rank, iproc, jproc, kproc) & BACK) {
        minor = PARALLEL_BOUNDARY;
        major = NO_SLIP;
    }
    else {
        minor = PARALLEL_BOUNDARY;
        major = PARALLEL_BOUNDARY;
    }
    for (int z = 0; z < zl + 2; ++z) {
        for (int x = 0; x < xl + 2; ++x) {
            flagField[fidx(sublength, x, 0, z)] = minor;
            flagField[fidx(sublength, x, yl + 1, z)] = major;
        }
    }

    /* YZ plane */
    if (checkBoundary(rank, iproc, jproc, kproc) & LEFT) {
        minor = NO_SLIP;
        major = PARALLEL_BOUNDARY;
    }
    else if (checkBoundary(rank, iproc, jproc, kproc) & RIGHT) {
        minor = PARALLEL_BOUNDARY;
        major = NO_SLIP;
    }
    else {
        minor = PARALLEL_BOUNDARY;
        major = PARALLEL_BOUNDARY;
    }
    for (int z = 0; z < zl + 2; ++z) {
        for (int y = 0; y < yl + 2; ++y) {
            flagField[fidx(sublength, 0, y, z)] = minor;
            flagField[fidx(sublength, xl + 1, y, z)] = major;
        }
    }

    /* Stream and Collide Fields are initialized to the respective lattice weights of the Cell */
    int index = 0;
    for (int z = 0; z < zl + 2; ++z) {
        for (int y = 0; y < yl + 2; ++y) {
            for (int x = 0; x < xl + 2; ++x) {
                for (int i = 0; i < Q; ++i) {
                    collideField[index] = LATTICEWEIGHTS[i];
                    streamField[index] = LATTICEWEIGHTS[i];
                    ++index;
                }
            }
        }
    }
}

