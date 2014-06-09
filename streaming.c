#include <mpi.h>
#include <stdio.h>

#include "streaming.h"
#include "LBDefinitions.h"

void exchangePdfs(double **sendBuffer, double **readBuffer, double *collideField,
        const int * const sublength, int rank, int iproc, int jproc, int kproc)
{
    MPI_Status status;
    const int xl = sublength[0];
    const int yl = sublength[1];
    const int zl = sublength[2];
    int i;
    /* ====================================
     * ========= X direction ==============
     * ==================================== */
    /* If the current rank does NOT process a right boundary subdomain, we must exchange PDFs
     * with the next (right) neighbor. */
#if 0
    printf("Checks for rank #%d\n", rank);
#endif
    if ((checkBoundary(rank, iproc, jproc, kproc) & RIGHT) == 0) {
        i = 0;
        /* Collect PDFs to send */
        for (int y = 0; y < yl + 2; ++y) {
            for (int z = 0; z < zl + 2; ++z) {
                /* Copy all right-pointing PDFs from rightmost lattices into the buffer */
                sendBuffer[1][i++] = collideField[idx(sublength, xl + 1, y, z, 3)];
                sendBuffer[1][i++] = collideField[idx(sublength, xl + 1, y, z, 7)];
                sendBuffer[1][i++] = collideField[idx(sublength, xl + 1, y, z, 10)];
                sendBuffer[1][i++] = collideField[idx(sublength, xl + 1, y, z, 13)];
                sendBuffer[1][i++] = collideField[idx(sublength, xl + 1, y, z, 17)];
            }
        }
        /* Send PDFs to right neighbor. Since we are ignoring rightmost
         * elements, it is safe to send the PDFs to rank + 1. */
        MPI_Send(&sendBuffer[1][0], 5 * (yl + 2) * (zl + 2), MPI_DOUBLE, rank + 1, 0,
                MPI_COMM_WORLD);

        /* Receive the elements of the right neighbor and set the left-facing
         * PDFs of the rightmost lattices of this rank to the buffer values.
         */
        MPI_Recv(&readBuffer[1][0], 5 * (yl + 2) * (zl + 2), MPI_DOUBLE, rank + 1, 0,
                MPI_COMM_WORLD, &status);
        i = 0;
        /* Assign received PDFs to lattices */
        for (int y = 0; y < yl + 2; ++y) {
            for (int z = 0; z < zl + 2; ++z) {
                collideField[idx(sublength, xl + 1, y, z, 1)] = readBuffer[1][i++];
                collideField[idx(sublength, xl + 1, y, z, 5)] = readBuffer[1][i++];
                collideField[idx(sublength, xl + 1, y, z, 8)] = readBuffer[1][i++];
                collideField[idx(sublength, xl + 1, y, z, 11)] = readBuffer[1][i++];
                collideField[idx(sublength, xl + 1, y, z, 15)] = readBuffer[1][i++];
            }
        }
    }
    /* If the current rank does NOT process a left boundary subdomain, we must exchange PDFs
     * with the previous (left) neighbor. */
    if ((checkBoundary(rank, iproc, jproc, kproc) & LEFT) == 0) {
        i = 0;
        for (int y = 0; y < yl + 2; ++y) {
            for (int z = 0; z < zl + 2; ++z) {
                /* Copy all right-pointing PDFs into the buffer */
                sendBuffer[0][i++] = collideField[idx(sublength, 0, y, z, 1)];
                sendBuffer[0][i++] = collideField[idx(sublength, 0, y, z, 5)];
                sendBuffer[0][i++] = collideField[idx(sublength, 0, y, z, 8)];
                sendBuffer[0][i++] = collideField[idx(sublength, 0, y, z, 11)];
                sendBuffer[0][i++] = collideField[idx(sublength, 0, y, z, 15)];
            }
        }
        /* Send PDFs to left neighbor. Since we are ignoring leftmost
         * elements, it is safe to send the PDFs to rank - 1. */
        MPI_Send(&sendBuffer[0][0], 5 * (yl + 2) * (zl + 2), MPI_DOUBLE, rank - 1, 0,
                MPI_COMM_WORLD);

        MPI_Recv(&readBuffer[0][0], 5 * (yl + 2) * (zl + 2), MPI_DOUBLE, rank - 1, 0,
                MPI_COMM_WORLD, &status);
        i = 0;
        /* Assign received PDFs to lattices */
        for (int y = 0; y < yl + 2; ++y) {
            for (int z = 0; z < zl + 2; ++z) {
                collideField[idx(sublength, 0, y, z, 3)] = readBuffer[0][i++];
                collideField[idx(sublength, 0, y, z, 7)] = readBuffer[0][i++];
                collideField[idx(sublength, 0, y, z, 10)] = readBuffer[0][i++];
                collideField[idx(sublength, 0, y, z, 13)] = readBuffer[0][i++];
                collideField[idx(sublength, 0, y, z, 17)] = readBuffer[0][i++];
            }
        }
    }
    /* ====================================
     * ========= Y direction ==============
     * ==================================== */
    /* If the current rank does NOT process a back boundary subdomain, we must exchange PDFs
     * with the previous neighbor. */
    if ((checkBoundary(rank, iproc, jproc, kproc) & BACK) == 0) {
        i = 0;
        for (int x = 0; x < xl + 2; ++x) {
            for (int z = 0; z < zl + 2; ++z) {
                /* Copy all back-pointing PDFs from the backmost lattices into the buffer */
                sendBuffer[5][i++] = collideField[idx(sublength, 0, yl + 1, z, 4)];
                sendBuffer[5][i++] = collideField[idx(sublength, 0, yl + 1, z, 11)];
                sendBuffer[5][i++] = collideField[idx(sublength, 0, yl + 1, z, 12)];
                sendBuffer[5][i++] = collideField[idx(sublength, 0, yl + 1, z, 13)];
                sendBuffer[5][i++] = collideField[idx(sublength, 0, yl + 1, z, 18)];
            }
        }
        /* Send PDFs to back neighbor. Since we are ignoring backmost
         * elements, it is safe to send the PDFs to rank + xlength. */
        MPI_Send(&sendBuffer[5][0], 5 * (xl + 2) * (zl + 2), MPI_DOUBLE, rank + iproc, 0,
                MPI_COMM_WORLD);
        /* Receive the elements of the back neighbor and set the front-facing
         * PDFs of the backmost lattices of this rank to the buffer values.
         */
        MPI_Recv(&readBuffer[5][0], 5 * (xl + 2) * (zl + 2), MPI_DOUBLE, rank + iproc, 0,
                MPI_COMM_WORLD, &status);
        i = 0;
        for (int x = 0; x < xl + 2; ++x) {
            for (int z = 0; z < zl + 2; ++z) {
                collideField[idx(sublength, x, yl + 1, z, 0)] = readBuffer[5][i++];
                collideField[idx(sublength, x, yl + 1, z, 5)] = readBuffer[5][i++];
                collideField[idx(sublength, x, yl + 1, z, 6)] = readBuffer[5][i++];
                collideField[idx(sublength, x, yl + 1, z, 7)] = readBuffer[5][i++];
                collideField[idx(sublength, x, yl + 1, z, 14)] = readBuffer[5][i++];
            }
        }
    }
    /* If the current rank does NOT process a front boundary subdomain, we must exchange PDFs
     * with the neighbor behind. */
    if ((checkBoundary(rank, iproc, jproc, kproc) & FRONT) == 0) {
        i = 0;
        for (int x = 0; x < xl + 2; ++x) {
            for (int z = 0; z < zl + 2; ++z) {
                /* Copy all front-pointing PDFs into the buffer */
                sendBuffer[4][i++] = collideField[idx(sublength, x, 0, z, 0)];
                sendBuffer[4][i++] = collideField[idx(sublength, x, 0, z, 5)];
                sendBuffer[4][i++] = collideField[idx(sublength, x, 0, z, 6)];
                sendBuffer[4][i++] = collideField[idx(sublength, x, 0, z, 7)];
                sendBuffer[4][i++] = collideField[idx(sublength, x, 0, z, 14)];
            }
        }
        /* Send PDFs to front neighbor. Since we are ignoring frontmost
         * elements, it is safe to send the PDFs to rank - xlength. */
        MPI_Send(&sendBuffer[4][0], 5 * (xl + 2) * (zl + 2), MPI_DOUBLE, rank - iproc, 0,
                MPI_COMM_WORLD);

        /* Receive the elements of the previous neighbor and set the back-facing
         * PDFs of the frontmost lattices of this rank to the buffer values.
         */
        MPI_Recv(&readBuffer[4][0], 5 * (xl + 2) * (zl + 2), MPI_DOUBLE, rank - iproc, 0,
                MPI_COMM_WORLD, &status);
        i = 0;
        for (int x = 0; x < xl + 2; ++x) {
            for (int z = 0; z < zl + 2; ++z) {
                collideField[idx(sublength, x, 0, z, 4)] = readBuffer[4][i++];
                collideField[idx(sublength, x, 0, z, 11)] = readBuffer[4][i++];
                collideField[idx(sublength, x, 0, z, 12)] = readBuffer[4][i++];
                collideField[idx(sublength, x, 0, z, 13)] = readBuffer[4][i++];
                collideField[idx(sublength, x, 0, z, 18)] = readBuffer[4][i++];
            }
        }
    }
    /* ====================================
     * ========= Z direction ==============
     * ==================================== */
    /* If the current rank does NOT process a bottom boundary subdomain, we must exchange
     * PDFs with the -xlength*xlength (lower) neighbor. */
    if ((checkBoundary(rank, iproc, jproc, kproc) & BOTTOM) == 0) {
        i = 0;
        for (int x = 0; x < xl + 2; ++x) {
            for (int y = 0; y < yl + 2; ++y) {
                /* Copy all down-pointing PDFs from bottom lattices into the buffer */
                sendBuffer[3][i++] = collideField[idx(sublength, x, y, 0, 0)];
                sendBuffer[3][i++] = collideField[idx(sublength, x, y, 0, 1)];
                sendBuffer[3][i++] = collideField[idx(sublength, x, y, 0, 2)];
                sendBuffer[3][i++] = collideField[idx(sublength, x, y, 0, 3)];
                sendBuffer[3][i++] = collideField[idx(sublength, x, y, 0, 4)];
            }
        }
        /* Send PDFs to lower neighbor. Since we are ignoring bottom
         * elements, it is safe to send the PDFs to rank - xlength*xlength. */
        MPI_Send(&sendBuffer[3][0], 5 * (xl + 2) * (yl + 2), MPI_DOUBLE, rank - iproc * jproc, 0,
                MPI_COMM_WORLD);

        /* Receive the elements of the lower neighbor and set the top-facing
         * PDFs of the lowest lattices of this rank to the buffer values.
         */
        MPI_Recv(&readBuffer[3][0], 5 * (xl + 2) * (yl + 2), MPI_DOUBLE, rank - iproc * jproc, 0,
                MPI_COMM_WORLD, &status);
        i = 0;
        for (int x = 0; x < xl + 2; ++x) {
            for (int y = 0; y < yl + 2; ++y) {
                collideField[idx(sublength, x, y, 0, 14)] = readBuffer[3][i++];
                collideField[idx(sublength, x, y, 0, 15)] = readBuffer[3][i++];
                collideField[idx(sublength, x, y, 0, 16)] = readBuffer[3][i++];
                collideField[idx(sublength, x, y, 0, 17)] = readBuffer[3][i++];
                collideField[idx(sublength, x, y, 0, 18)] = readBuffer[3][i++];
            }
        }
    }
    /* If the current rank does NOT process a top boundary subdomain, we must exchange
     * PDFs with the +xlength*xlength (upper) neighbor. */
    if ((checkBoundary(rank, iproc, jproc, kproc) & TOP) == 0) {
        i = 0;
        for (int x = 0; x < xl + 2; ++x) {
            for (int y = 0; y < yl + 2; ++y) {
                /* Copy all up-pointing PDFs from the top lattices into the buffer */
                sendBuffer[2][i++] = collideField[idx(sublength, x, y, zl + 1, 14)];
                sendBuffer[2][i++] = collideField[idx(sublength, x, y, zl + 1, 15)];
                sendBuffer[2][i++] = collideField[idx(sublength, x, y, zl + 1, 16)];
                sendBuffer[2][i++] = collideField[idx(sublength, x, y, zl + 1, 17)];
                sendBuffer[2][i++] = collideField[idx(sublength, x, y, zl + 1, 18)];

            }
        }
        /* Send PDFs to upper neighbor. Since we are ignoring top
         * elements, it is safe to send the PDFs to rank + xlength*xlength. */
        MPI_Send(&sendBuffer[2][0], 5 * (xl + 2) * (yl + 2), MPI_DOUBLE, rank + iproc * jproc, 0,
                MPI_COMM_WORLD);

        /* Receive the elements of the upper neighbor and set the down-facing
         * PDFs of the highest lattices of this rank to the buffer values.
         */
        MPI_Recv(&readBuffer[2][0], 5 * (xl + 2) * (yl + 2), MPI_DOUBLE, rank + iproc * jproc, 0,
                MPI_COMM_WORLD, &status);
        i = 0;
        for (int x = 0; x < xl + 2; ++x) {
            for (int y = 0; y < yl + 2; ++y) {
                collideField[idx(sublength, x, y, zl + 1, 0)] = readBuffer[2][i++];
                collideField[idx(sublength, x, y, zl + 1, 1)] = readBuffer[2][i++];
                collideField[idx(sublength, x, y, zl + 1, 2)] = readBuffer[2][i++];
                collideField[idx(sublength, x, y, zl + 1, 3)] = readBuffer[2][i++];
                collideField[idx(sublength, x, y, zl + 1, 4)] = readBuffer[2][i++];
            }
        }
    }
}

void doStreaming(double *collideField, double *streamField, int *flagField,
        const int * const sublength)
{
    int dx, dy, dz;
    double fi;
    /*Setting distribution function for each moving direction/lattice velocity of every particle*/
    for (int z = 1; z < sublength[0] + 1; ++z) {
        for (int y = 1; y < sublength[1] + 1; ++y) {
            for (int x = 1; x < sublength[2] + 1; ++x) {
                for (int i = 0; i < Q; ++i) {

                    /*dx = c_i_x*dt, dt = 1*/
                    dx = LATTICEVELOCITIES[i][0];
                    dy = LATTICEVELOCITIES[i][1];
                    dz = LATTICEVELOCITIES[i][2];

                    /*New value for our distribution function (DF) of the index 'i'

                     (We set it to DF(i) of the next particle, whose i-th lattice velocity
                     points towards considered particle (x,y,z))

                     Position of that next particle is given by (x-dx, y-dy, z-dz)*/

                    fi = collideField[idx(sublength, x - dx, y - dy, z - dz, i)];
                    /*fi = collideField[Q * ((z-dz)*(sublength+2)*(sublength+2) + (y-dy)*(sublength+2) + x-dx) + i];*/
                    streamField[idx(sublength, x, y, z, i)] = fi;
                    /*streamField[Q * (z*(sublength+2)*(sublength+2) + y*(sublength+2) + x) + i] = fi;*/
                }
            }
        }
    }

}

