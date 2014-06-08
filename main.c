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

int main(int argc, char *argv[])
{
    double *collideField = NULL;
    double *streamField = NULL;
    int *flagField = NULL;
    /* Domain length (constant => cubic domain) */
    int xlength;
    /* Relaxation parameter for BGK approximation */
    double tau;
    /* 3-dimensional wall velocity for lid */
    double velocityWall[D];
    /* Timesteps to perform */
    int timesteps;
    /* Plots are generated every timestepsPerPlotting time steps */
    int timestepsPerPlotting;
    /* MPI */
    /* Rank of this process */
    int rank;
    /* Number of MPI ranks */
    int number_of_ranks;
    /* Subdivision of each domain dimension */
    int iproc, jproc, kproc;
    /* Send and receive buffer for exchanging PDFs among processes */
    double *sendBuffer[6];
    double *readBuffer[6];

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_ranks);

    if (rank == 0) {
        printf("LBM simulation by Krivokapic, Mody, Malcher - CFD Lab SS2014\n");
        printf("============================================================\n");
        printf("Reading in parameters...\n");
        if (readParameters(&xlength, &tau, velocityWall, &iproc, &jproc, &kproc, &timesteps,
                &timestepsPerPlotting, argc, argv)) {
            printf("Reading in parameters failed. Aborting program!\n");
            fflush(stdout);
            MPI_Finalize();
            exit(EXIT_FAILURE);
        }
        printf("...done\n");
        printf("Starting LBM with MPI using %d ranks...\n", number_of_ranks);
    }
    broadcastInitialValues(&xlength, &tau, velocityWall, &iproc, &jproc, &kproc, &timesteps,
            &timestepsPerPlotting);
    /* TODO: Add checks for ijk==ranks and xlength%ijk != 0 */
    if (iproc * jproc * kproc != number_of_ranks) {
        if (rank == 0) {
            printf("iproc*jproc*kproc != number of ranks! Aborting program!\n");
            fflush(stdout);
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
    /* TODO: Broadcast values that were read from file to other processes*/

    /* Compute cuboid lengths for subdomains which are handed over to the processes */
#if 0
    printf("Values of proc #%d: xlength=%d, ijkproc=(%d, %d, %d), timesteps=%d\n",
            rank, xlength, iproc, jproc, kproc, timesteps);
#endif
    /* Length of subdomains (NOT including ghost cells) */
    const int sublength[3] = { xlength / iproc, xlength / jproc, xlength / kproc };
    /*const int sublength[3] = { xlength / iproc + 2, xlength / jproc + 2, xlength / kproc + 2 };*/
    /* Subdomain volume including ghost layer */
    const int volume_ghost = (sublength[0] + 2) * (sublength[1] + 2) * (sublength[2] + 2);
    collideField = malloc(sizeof(*collideField) * Q * volume_ghost);
    streamField = malloc(sizeof(*streamField) * Q * volume_ghost);
    flagField = malloc(sizeof(*flagField) * volume_ghost);
#if 0
    printf("Sublength = (%d,%d,%d), process %d\n", sublength[0], sublength[1], sublength[2], rank);
#endif
    /* Initialize pointers */
    initialiseFields(collideField, streamField, flagField, sublength, rank, iproc, jproc, kproc);

    const int xl = sublength[0];
    const int yl = sublength[1];
    const int zl = sublength[2];
    /* Send buffers */
    sendBuffer[0] = malloc(5 * sizeof(double) * (yl+2)*(zl+2));
    sendBuffer[1] = malloc(5 * sizeof(double) * (yl+2)*(zl+2));
    sendBuffer[2] = malloc(5 * sizeof(double) * (xl+2)*(yl+2));
    sendBuffer[3] = malloc(5 * sizeof(double) * (xl+2)*(yl+2));
    sendBuffer[4] = malloc(5 * sizeof(double) * (xl+2)*(zl+2));
    sendBuffer[5] = malloc(5 * sizeof(double) * (xl+2)*(zl+2));
    /* Read buffers */
    readBuffer[0] = malloc(5 * sizeof(double) * (yl+2)*(zl+2));
    readBuffer[1] = malloc(5 * sizeof(double) * (yl+2)*(zl+2));
    readBuffer[2] = malloc(5 * sizeof(double) * (xl+2)*(yl+2));
    readBuffer[3] = malloc(5 * sizeof(double) * (xl+2)*(yl+2));
    readBuffer[4] = malloc(5 * sizeof(double) * (xl+2)*(zl+2));
    readBuffer[5] = malloc(5 * sizeof(double) * (xl+2)*(zl+2));

    /*MPI_Barrier(MPI_COMM_WORLD);*/
    for (int t = 0; t < timesteps; ++t) {
        /* TODO: Extraction, Swap, Injection: Check if the right PDF field (collideField) was chosen
         * 1. Extraction step: Extract relevant pdfs to be sent. One buffer for all sides
         * => Buffer size = 5(c_i's) * 6(lattice planes) = 30 pdfs for each send/recv
         * 2. Swap step: Send and receive pdfs. Send from i to j must be followed by recv in i from j.
         * => 1. Send i->j
         *    2. Recv j->i
         * 3. Injection step: Write received pdfs into the corresponding boundary layers
         */
#if 0
        printf("Rank #%d of %d before exchanging\n", rank, number_of_ranks);
        fflush(stdout);
#endif
        exchangePdfs(sendBuffer, readBuffer, collideField, sublength, rank, iproc, jproc, kproc);

        doStreaming(collideField, streamField, flagField, sublength);

        double *swap = NULL;
        swap = collideField;
        collideField = streamField;
        streamField = swap;

        doCollision(collideField, flagField, &tau, sublength);
        treatBoundary(collideField, flagField, velocityWall, sublength);
#if 1
        /* Write output to vtk file for postprocessing */
        if (t % timestepsPerPlotting == 0) {
            if (rank == 0) {
                printf("%d %%\r", (int) ((double) t / timesteps * 100));
                fflush(stdout);
            }
            writeVtkOutput(collideField, flagField, "lbm_out", rank, iproc, jproc, kproc, t,
                    sublength);
        }
#endif
        /* TODO: Wait for all processes to finish before starting with new loop */
        /*MPI_Barrier(MPI_COMM_WORLD);*/
    }

    if (rank == 0) {
        printf("Done!\n============================================================\n");
        printf("LBM simulation completed for %d timesteps\n", timesteps);
        printf("Freeing allocated memory...\n");
    }
    free(collideField);
    free(streamField);
    free(flagField);
    for (int i = 0; i < 6; ++i) {
        free(sendBuffer[i]);
        free(readBuffer[i]);
    }
    MPI_Finalize();
    return 0;
}

#endif
