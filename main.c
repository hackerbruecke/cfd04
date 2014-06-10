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
    /* Time steps to perform */
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
        /* Error indicator */
        int error = 0;
        printf("LBM simulation by Krivokapic, Mody, Malcher - CFD Lab SS2014\n");
        printf("============================================================\n");
        printf("Reading in parameters...\n");
        if (readParameters(&xlength, &tau, velocityWall, &iproc, &jproc, &kproc, &timesteps,
                &timestepsPerPlotting, argc, argv)) {
            printf("Reading in parameters failed. Aborting program!\n");
            fflush(stdout);
            error = 1;
        }
        else if (iproc * jproc * kproc != number_of_ranks) {
            printf("iproc*jproc*kproc != number of ranks! Aborting program!\n");
            fflush(stdout);
            error = 1;
        }
        else if (xlength % iproc || xlength % jproc || xlength % kproc) {
            printf("iproc/jproc/kproc must divide xlength without remainder!\n");
            fflush(stdout);
            error = 1;
        }
        if (error) {
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }
        /* No errors, continue */
        printf("...done\n");
        printf("Starting LBM with MPI using %d ranks...\n", number_of_ranks);
    }
    /* Broadcast values read in by rank 0 to others */
    broadcastInitialValues(&xlength, &tau, velocityWall, &iproc, &jproc, &kproc, &timesteps,
            &timestepsPerPlotting);
    /* Length of subdomains (NOT including ghost cells) */
    const int sublength[3] = { xlength / iproc, xlength / jproc, xlength / kproc };
    /* Subdomain volume including ghost layer */
    const int volume_ghost = (sublength[0] + 2) * (sublength[1] + 2) * (sublength[2] + 2);
    collideField = malloc(sizeof(*collideField) * Q * volume_ghost);
    streamField  = malloc(sizeof(*streamField)  * Q * volume_ghost);
    flagField    = malloc(sizeof(*flagField)    * volume_ghost);

    /* Initialize pointers for PDF fields */
    initialiseFields(collideField, streamField, flagField, sublength, rank, iproc, jproc, kproc);
    /* Allocate space for send/receive buffer */
    allocateBufferSpace(sendBuffer, readBuffer, sublength);

    for (int t = 0; t < timesteps; ++t) {
        /*
         * 1. Extraction step: Extract relevant pdfs to be sent. One buffer for all sides
         * => Buffer size = 5(c_i's) * 6(lattice planes)*elements = 30*elements pdfs for each send/recv
         * 2. Swap step: Send and receive pdfs. Send from i to j must be followed by recv in i from j.
         * => 1. Send i->j
         *    2. Recv j->i
         * 3. Injection step: Write received pdfs into the corresponding boundary layers
         */
        exchangePdfs(sendBuffer, readBuffer, collideField, sublength, rank, iproc, jproc, kproc);
        /* Do the actual streaming step after PDFs from parallel boundaries were received. */
        doStreaming(collideField, streamField, flagField, sublength);

        double *swap = NULL;
        swap = collideField;
        collideField = streamField;
        streamField = swap;

        doCollision(collideField, flagField, &tau, sublength);
        treatBoundary(collideField, flagField, velocityWall, sublength);

        /* Write output to vtk file for postprocessing */
        if (t % timestepsPerPlotting == 0) {
            writeVtkOutput(collideField, flagField, "lbm_out", rank, iproc, jproc, kproc, t,
                    sublength);
            if (rank == 0) {
                printf("%d %%\r", (int) ((double) t / timesteps * 100));
                fflush(stdout);
                writeParallelVtkFile("lbm_out", iproc, jproc, kproc, sublength, t);
            }
        }
    }

    if (rank == 0) {
        printf("Done!\n============================================================\n");
        printf("LBM simulation completed for %d timesteps\n", timesteps);
        printf("Freeing allocated memory...\n");
    }
    cleanup(collideField, streamField, flagField, sendBuffer, readBuffer);
    MPI_Finalize();
    return 0;
}

#endif
