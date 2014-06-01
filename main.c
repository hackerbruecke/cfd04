#ifndef _MAIN_C_
#define _MAIN_C_

#include <mpi.h>
#include <stdio.h>

#include "LBDefinitions.h"
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "visualLB.h"
#include "boundary.h"
#include "initLB.h"

int main(int argc, char *argv[]) {
	double 	*collideField = NULL;
	double	*streamField = NULL;
	int 	*flagField = NULL;
	int 	xlength;				/* Domain length (constant => cubic domain) */
	double 	tau;					/* Relaxation parameter for BGK approximation */
	double 	velocityWall[D];		/* 3-dimensional wall velocity for lid */
	int 	timesteps;				/* Timesteps to perform */
	int 	timestepsPerPlotting;	/* Plots are generated every timestepsPerPlotting time steps */
	/* MPI */
	int 	rank;					/* Rank of this process */
	int 	number_of_ranks;		/* Number of MPI ranks */
	int 	iproc, jproc, kproc;	/* Subdivision of each domain dimension */
	double 	*sendBuffer[6];			/* Send and receive buffer for exchanging PDFs among processes */
	double 	*readBuffer[6];

	/* Initialize MPI */
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &number_of_ranks);

	if (rank == 0) {
		printf("LBM simulation by Krivokapic, Mody, Malcher - CFD Lab SS2014\n");
		printf("============================================================\n");
		printf("Reading in parameters...\n");
		if (readParameters(&xlength, &tau, velocityWall, &iproc, &jproc, &kproc, &timesteps,
				&timestepsPerPlotting, argc, argv))
		{
			printf("Reading in parameters failed. Aborting program!\n");
			exit(EXIT_FAILURE);
		}
		printf("...done\n");
		printf("Starting LBM with MPI using %d ranks...\n", number_of_ranks);
	}
	/* TODO: Broadcast values that were read from file to other processes*/
	broadcastInitialValues(&xlength, &tau, velocityWall, &iproc, &jproc, &kproc, &timesteps);

    /* TODO: Change the memory allocation according to the lattice sizes for each MPI process */
	/* Compute cuboid lengths for subdomains which are handed over to the processes */
	int xl = xlength / iproc + 2;
	int yl = xlength / jproc + 2;
	int zl = xlength / kproc + 2;

	/* Subdomain volume including ghost layer */
	const int volume_ghost = (xl+2)*(yl+2)*(zl+2);
	collideField = malloc(sizeof(*collideField) * Q * volume_ghost);
	streamField  = malloc(sizeof(*streamField) * Q * volume_ghost);
	flagField    = malloc(sizeof(*flagField) * volume_ghost);
	/* TODO: MPI version */
    /* Initialize pointers */
	initialiseFields(collideField, streamField, flagField, xl, yl, zl, rank, number_of_ranks);
    
    for (int t = 0; t < timesteps; ++t) {
		double *swap = NULL;
    	/* TODO: Extraction, Swap, Injection
    	 * 1. Extraction step: Extract relevant pdfs to be sent. One buffer for all sides
    	 * => Buffer size = 5(c_i's) * 6(lattice planes) = 30 pdfs for each send/recv
    	 * 2. Swap step: Send and receive pdfs. Send from i to j must be followed by recv in i from j.
    	 * => 1. Send i->j
    	 *    2. Recv j->i
    	 * 3. Injection step: Write received pdfs into the corresponding boundary layers
    	 */
		memset(sendBuffer, 0, 6 * sizeof(*sendBuffer));
		memset(readBuffer, 0, 6 * sizeof(*readBuffer));


		doStreaming(collideField, streamField, flagField, xlength);
		swap = collideField;
		collideField = streamField;
		streamField = swap;

		doCollision(collideField, flagField, &tau, xlength);
		treatBoundary(collideField, flagField, velocityWall, xlength);
        
        /* Write output to vtk file for postprocessing */
/*
		if (rank == 0 && t % timestepsPerPlotting == 0) {
			printf("%d %%\r", (int)((double)t/timesteps*100));
			fflush(stdout);
			writeVtkOutput(collideField, flagField, "lbm_out", t, xlength);
		}
*/
		/* TODO: Wait for all processes to finish before starting with new loop */
		MPI_Barrier(MPI_COMM_WORLD);
	}

    if (rank == 0) {
		printf("Done!\n============================================================\n");
		printf("LBM simulation completed for timesteps\n", timesteps);
		printf("Freeing allocated memory...\n");
    }
    free(collideField);
    free(streamField);
    free(flagField);
    MPI_Finalize();
	return 0;
}

#endif
