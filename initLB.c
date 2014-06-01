#include <mpi.h>

#include "initLB.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *tau, double *velocityWall, int *iproc, int *jproc, int *kproc,
		int *timesteps, int *timestepsPerPlotting, int argc, char **argv){
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

void broadcastInitialValues(int *xlength, double *tau, double *velocityWall,
		int *iproc, int *jproc, int *kproc, int *timesteps) {
	MPI_Bcast(xlength, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	/* Note: 3 doubles for velocityWall */
	MPI_Bcast(velocityWall, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(iproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(jproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(kproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(timesteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
}


void initialiseFields(double *collideField, double *streamField, int *flagField, int xl, int yl, int zl, int rank, int number_of_ranks){
    /* Set all values of flagField to 0, i.e. FLUID state. We will apply boundary conditions
        in the following */
    memset(flagField, FLUID, (xl+2)*(yl+2)*(zl+2) * sizeof(*flagField));

    /* TODO: If the rank of the process indicates that a boundary is being treated, set the respective conditions */
#if 0
    for (int x = 0; x < xl + 2; ++x) {
        for (int y = 0; y < yl + 2; ++y) {
            flagField[fidx(xlength, x, y, 0)] = NO_SLIP;
            flagField[fidx(xlength, x, y, xlength+1)] = MOVING_WALL;
        }
    }

/* The values for Boundary on Y = 0 and Ymax plane set to No_Slip */
    for (int x = 0; x < xl + 2; ++x) {
        for (int z = 0; z < zl + 2; ++z) {
            flagField[fidx(xlength, x, 0, z)] = NO_SLIP;
            flagField[fidx(xlength, x, xlength+1, z)] = NO_SLIP;
        }
    }
/* The values for Boundary on Y = 0 and Ymax plane set to No_Slip */
    for (int y = 0; y < xl + 2; ++y) {
        for (int z = 0; z < zl + 2; ++z) {
            flagField[fidx(xlength, 0, y, z)] = NO_SLIP;
            flagField[fidx(xlength, xlength+1, y, z)] = NO_SLIP;
        }
    }
#endif
    /* Stream and Collide Fields are initialized to the respective lattice weights of the Cell */
    const int lengths[3] = {xl, yl, zl};
    for (int x = 0; x < xl + 2; ++x) {
        for (int y = 0; y < yl + 2; ++y) {
            for (int z = 0; z < zl + 2; ++z) {
                for (int i = 0; i < Q; ++i) {
                    collideField[idx(lengths, x, y, z, i)] = LATTICEWEIGHTS[i];
                    streamField[idx(lengths, x, y, z, i)] = LATTICEWEIGHTS[i];
                }
            }
        }
    }
}

