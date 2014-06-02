#include "initLB.h"
#include "LBDefinitions.h"

int readParameters(int *xlength, double *Tau, double *velocityWall,
		int *Timesteps, int *Time_per_plot, int argc, char *argv[]) {
	/* Checking if there are exactly 2 parameters in the input. Else give Error message. */
	if (argc != 2) {
		printf("You must provide a filename!\n");
		return -1;
	}
	const char *filename = argv[1];
	/* Reading all the input parameters */
	read_int(filename, "Length_in_X", &xlength[0]);
	read_int(filename, "Length_in_Y", &xlength[1]);
	read_int(filename, "Length_in_Z", &xlength[2]);
	read_double(filename, "Wall_Vel_X", &velocityWall[0]);
	read_double(filename, "Wall_Vel_Y", &velocityWall[1]);
	read_double(filename, "Wall_Vel_Z", &velocityWall[2]);
	READ_DOUBLE(filename, *Tau);
	READ_INT(filename, *Timesteps);
	READ_INT(filename, *Time_per_plot);

	return 0;
}

void initialiseFields(double *collideField, double *streamField, int *flagField,
		int *xlength) {

	/* Set all values of flagField to 0, i.e. FLUID state. We will apply boundary conditions
	 in the following */
	memset(flagField, FLUID,
			(xlength[2] + 2) * (xlength[1] + 2) * (xlength[0] + 2)
					* sizeof(*flagField));

	/* The values for Boundary on Z = 0 set to No_Slip and Zmax plane set to Moving_Wall */
	for (int x = 0; x < xlength[0] + 2; ++x) {
		for (int y = 0; y < xlength[1] + 2; ++y) {
			flagField[fidx(xlength, x, y, 0)] = NO_SLIP;
			flagField[fidx(xlength, x, y, xlength[2] + 1)] = MOVING_WALL;
		}
	}

	/* The values for Boundary on Y = 0 and Ymax plane set to No_Slip */
	for (int x = 0; x < xlength[0] + 2; ++x) {
		for (int z = 0; z < xlength[2] + 2; ++z) {
			flagField[fidx(xlength, x, 0, z)] = NO_SLIP;
			flagField[fidx(xlength, x, xlength[1] + 1, z)] = NO_SLIP;
		}
	}
	/* The values for Boundary on Y = 0 and Ymax plane set to No_Slip */
	for (int y = 0; y < xlength[1] + 2; ++y) {
		for (int z = 0; z < xlength[2] + 2; ++z) {
			flagField[fidx(xlength, 0, y, z)] = NO_SLIP;
			flagField[fidx(xlength, xlength[0] + 1, y, z)] = NO_SLIP;
		}
	}
	/* Stream and Collide Fields are initialized to the respective Lattice-Weights of the Cell */
	for (int z = 0; z < xlength[2] + 2; ++z) {
		for (int y = 0; y < xlength[1] + 2; ++y) {
			for (int x = 0; x < xlength[0] + 2; ++x) {
				for (int i = 0; i < Q; ++i) {
					collideField[idx(xlength, x, y, z, i)] = LATTICEWEIGHTS[i];
					streamField[idx(xlength, x, y, z, i)] = LATTICEWEIGHTS[i];
				}
			}
		}
	}
}

