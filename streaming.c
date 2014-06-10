#include <mpi.h>

#include "streaming.h"
#include "LBDefinitions.h"

/* Helper enum to identify neighbor direction. Note: Do not replace with RIGHT, LEFT etc.
 * as the latter enums contain powers of two instead of subsequent numbers
 */
typedef enum
{
    LE, RI, TO, BO, FR, BA
} DIRECTION;

void doRightSendRecv(double **sendBuffer, double **readBuffer, double *collideField,
        const int * const sublength, int rank)
{
    MPI_Status status;
    const int xl = sublength[0];
    const int yl = sublength[1];
    const int zl = sublength[2];
    int i = 0;
    /* Collect PDFs to send */
    for (int y = 0; y < yl + 2; ++y) {
        for (int z = 0; z < zl + 2; ++z) {
            /* Copy all right-pointing PDFs from rightmost lattices into the buffer */
            sendBuffer[RI][i++] = collideField[idx(sublength, xl, y, z, 3)];
            sendBuffer[RI][i++] = collideField[idx(sublength, xl, y, z, 7)];
            sendBuffer[RI][i++] = collideField[idx(sublength, xl, y, z, 10)];
            sendBuffer[RI][i++] = collideField[idx(sublength, xl, y, z, 13)];
            sendBuffer[RI][i++] = collideField[idx(sublength, xl, y, z, 17)];
        }
    }
    /* Send PDFs to right neighbor. Since we are ignoring rightmost
     * elements, it is safe to send the PDFs to rank + 1. */
    MPI_Send(&sendBuffer[RI][0], 5 * (yl + 2) * (zl + 2), MPI_DOUBLE, rank + 1, 0,
    MPI_COMM_WORLD);

    /* Receive the elements of the right neighbor and set the left-facing
     * PDFs of the rightmost lattices of this rank to the buffer values.
     */
    MPI_Recv(&readBuffer[RI][0], 5 * (yl + 2) * (zl + 2), MPI_DOUBLE, rank + 1, 0,
    MPI_COMM_WORLD, &status);
    i = 0;
    /* Assign received PDFs to lattices */
    for (int y = 0; y < yl + 2; ++y) {
        for (int z = 0; z < zl + 2; ++z) {
            collideField[idx(sublength, xl + 1, y, z, 1)] = readBuffer[RI][i++];
            collideField[idx(sublength, xl + 1, y, z, 5)] = readBuffer[RI][i++];
            collideField[idx(sublength, xl + 1, y, z, 8)] = readBuffer[RI][i++];
            collideField[idx(sublength, xl + 1, y, z, 11)] = readBuffer[RI][i++];
            collideField[idx(sublength, xl + 1, y, z, 15)] = readBuffer[RI][i++];
        }
    }
}

void doLeftSendRecv(double **sendBuffer, double **readBuffer, double *collideField,
        const int * const sublength, int rank)
{
    MPI_Status status;
    const int yl = sublength[1];
    const int zl = sublength[2];
    int i = 0;
    for (int y = 0; y < yl + 2; ++y) {
        for (int z = 0; z < zl + 2; ++z) {
            /* Copy all left-pointing PDFs into the buffer */
            sendBuffer[LE][i++] = collideField[idx(sublength, 1, y, z, 1)];
            sendBuffer[LE][i++] = collideField[idx(sublength, 1, y, z, 5)];
            sendBuffer[LE][i++] = collideField[idx(sublength, 1, y, z, 8)];
            sendBuffer[LE][i++] = collideField[idx(sublength, 1, y, z, 11)];
            sendBuffer[LE][i++] = collideField[idx(sublength, 1, y, z, 15)];
        }
    }
    /* Send PDFs to left neighbor. Since we are ignoring leftmost
     * elements, it is safe to send the PDFs to rank - 1. */
    MPI_Send(&sendBuffer[LE][0], 5 * (yl + 2) * (zl + 2), MPI_DOUBLE, rank - 1, 0,
    MPI_COMM_WORLD);

    MPI_Recv(&readBuffer[LE][0], 5 * (yl + 2) * (zl + 2), MPI_DOUBLE, rank - 1, 0,
    MPI_COMM_WORLD, &status);
    i = 0;
    /* Assign received PDFs to lattices */
    for (int y = 0; y < yl + 2; ++y) {
        for (int z = 0; z < zl + 2; ++z) {
            collideField[idx(sublength, 0, y, z, 3)] = readBuffer[LE][i++];
            collideField[idx(sublength, 0, y, z, 7)] = readBuffer[LE][i++];
            collideField[idx(sublength, 0, y, z, 10)] = readBuffer[LE][i++];
            collideField[idx(sublength, 0, y, z, 13)] = readBuffer[LE][i++];
            collideField[idx(sublength, 0, y, z, 17)] = readBuffer[LE][i++];
        }
    }
}

void doBackSendRecv(double **sendBuffer, double **readBuffer, double *collideField,
        const int * const sublength, int rank, int iproc)
{
    MPI_Status status;
    const int xl = sublength[0];
    const int yl = sublength[1];
    const int zl = sublength[2];
    int i = 0;
    for (int x = 0; x < xl + 2; ++x) {
        for (int z = 0; z < zl + 2; ++z) {
            /* Copy all back-pointing PDFs from the backmost lattices into the buffer */
            sendBuffer[BA][i++] = collideField[idx(sublength, x, yl, z, 4)];
            sendBuffer[BA][i++] = collideField[idx(sublength, x, yl, z, 11)];
            sendBuffer[BA][i++] = collideField[idx(sublength, x, yl, z, 12)];
            sendBuffer[BA][i++] = collideField[idx(sublength, x, yl, z, 13)];
            sendBuffer[BA][i++] = collideField[idx(sublength, x, yl, z, 18)];
        }
    }
    /* Send PDFs to back neighbor. Since we are ignoring backmost
     * elements, it is safe to send the PDFs to rank + xlength. */
    MPI_Send(&sendBuffer[BA][0], 5 * (xl + 2) * (zl + 2), MPI_DOUBLE, rank + iproc, 0,
    MPI_COMM_WORLD);
    /* Receive the elements of the back neighbor and set the front-facing
     * PDFs of the backmost lattices of this rank to the buffer values.
     */
    MPI_Recv(&readBuffer[BA][0], 5 * (xl + 2) * (zl + 2), MPI_DOUBLE, rank + iproc, 0,
    MPI_COMM_WORLD, &status);
    i = 0;
    for (int x = 0; x < xl + 2; ++x) {
        for (int z = 0; z < zl + 2; ++z) {
            collideField[idx(sublength, x, yl+1, z, 0)] = readBuffer[BA][i++];
            collideField[idx(sublength, x, yl+1, z, 5)] = readBuffer[BA][i++];
            collideField[idx(sublength, x, yl+1, z, 6)] = readBuffer[BA][i++];
            collideField[idx(sublength, x, yl+1, z, 7)] = readBuffer[BA][i++];
            collideField[idx(sublength, x, yl+1, z, 14)] = readBuffer[BA][i++];
        }
    }
}

void doFrontSendRecv(double **sendBuffer, double **readBuffer, double *collideField,
        const int * const sublength, int rank, int iproc)
{
    MPI_Status status;
    const int xl = sublength[0];
    const int zl = sublength[2];
    int i = 0;
    for (int x = 0; x < xl + 2; ++x) {
        for (int z = 0; z < zl + 2; ++z) {
            /* Copy all front-pointing PDFs into the buffer */
            sendBuffer[FR][i++] = collideField[idx(sublength, x, 1, z, 0)];
            sendBuffer[FR][i++] = collideField[idx(sublength, x, 1, z, 5)];
            sendBuffer[FR][i++] = collideField[idx(sublength, x, 1, z, 6)];
            sendBuffer[FR][i++] = collideField[idx(sublength, x, 1, z, 7)];
            sendBuffer[FR][i++] = collideField[idx(sublength, x, 1, z, 14)];
        }
    }
    /* Send PDFs to front neighbor. Since we are ignoring frontmost
     * elements, it is safe to send the PDFs to rank - xlength. */
    MPI_Send(&sendBuffer[FR][0], 5 * (xl + 2) * (zl + 2), MPI_DOUBLE, rank - iproc, 0,
    MPI_COMM_WORLD);

    /* Receive the elements of the previous neighbor and set the back-facing
     * PDFs of the frontmost lattices of this rank to the buffer values.
     */
    MPI_Recv(&readBuffer[FR][0], 5 * (xl + 2) * (zl + 2), MPI_DOUBLE, rank - iproc, 0,
    MPI_COMM_WORLD, &status);
    i = 0;
    for (int x = 0; x < xl + 2; ++x) {
        for (int z = 0; z < zl + 2; ++z) {
            collideField[idx(sublength, x, 0, z, 4)] = readBuffer[FR][i++];
            collideField[idx(sublength, x, 0, z, 11)] = readBuffer[FR][i++];
            collideField[idx(sublength, x, 0, z, 12)] = readBuffer[FR][i++];
            collideField[idx(sublength, x, 0, z, 13)] = readBuffer[FR][i++];
            collideField[idx(sublength, x, 0, z, 18)] = readBuffer[FR][i++];
        }
    }
}

void doBottomSendRecv(double **sendBuffer, double **readBuffer, double *collideField,
        const int * const sublength, int rank, int iproc, int jproc)
{
    MPI_Status status;
    const int xl = sublength[0];
    const int yl = sublength[1];

    int i = 0;
    for (int x = 0; x < xl + 2; ++x) {
        for (int y = 0; y < yl + 2; ++y) {
            /* Copy all down-pointing PDFs from bottom lattices into the buffer */
            sendBuffer[BO][i++] = collideField[idx(sublength, x, y, 1, 0)];
            sendBuffer[BO][i++] = collideField[idx(sublength, x, y, 1, 1)];
            sendBuffer[BO][i++] = collideField[idx(sublength, x, y, 1, 2)];
            sendBuffer[BO][i++] = collideField[idx(sublength, x, y, 1, 3)];
            sendBuffer[BO][i++] = collideField[idx(sublength, x, y, 1, 4)];
        }
    }
    /* Send PDFs to lower neighbor. Since we are ignoring bottom
     * elements, it is safe to send the PDFs to rank - xlength*xlength. */
    MPI_Send(&sendBuffer[BO][0], 5 * (xl + 2) * (yl + 2), MPI_DOUBLE, rank - iproc * jproc, 0,
    MPI_COMM_WORLD);

    /* Receive the elements of the lower neighbor and set the top-facing
     * PDFs of the lowest lattices of this rank to the buffer values.
     */
    MPI_Recv(&readBuffer[BO][0], 5 * (xl + 2) * (yl + 2), MPI_DOUBLE, rank - iproc * jproc, 0,
    MPI_COMM_WORLD, &status);
    i = 0;
    for (int x = 0; x < xl + 2; ++x) {
        for (int y = 0; y < yl + 2; ++y) {
            collideField[idx(sublength, x, y, 0, 14)] = readBuffer[BO][i++];
            collideField[idx(sublength, x, y, 0, 15)] = readBuffer[BO][i++];
            collideField[idx(sublength, x, y, 0, 16)] = readBuffer[BO][i++];
            collideField[idx(sublength, x, y, 0, 17)] = readBuffer[BO][i++];
            collideField[idx(sublength, x, y, 0, 18)] = readBuffer[BO][i++];
        }
    }
}

void doTopSendRecv(double **sendBuffer, double **readBuffer, double *collideField,
        const int * const sublength, int rank, int iproc, int jproc)
{
    MPI_Status status;
    const int xl = sublength[0];
    const int yl = sublength[1];
    const int zl = sublength[2];
    int i = 0;

    for (int x = 0; x < xl + 2; ++x) {
        for (int y = 0; y < yl + 2; ++y) {
            /* Copy all up-pointing PDFs from the top lattices into the buffer */
            sendBuffer[TO][i++] = collideField[idx(sublength, x, y, zl, 14)];
            sendBuffer[TO][i++] = collideField[idx(sublength, x, y, zl, 15)];
            sendBuffer[TO][i++] = collideField[idx(sublength, x, y, zl, 16)];
            sendBuffer[TO][i++] = collideField[idx(sublength, x, y, zl, 17)];
            sendBuffer[TO][i++] = collideField[idx(sublength, x, y, zl, 18)];

        }
    }
    /* Send PDFs to upper neighbor. Since we are ignoring top
     * elements, it is safe to send the PDFs to rank + xlength*xlength. */
    MPI_Send(&sendBuffer[TO][0], 5 * (xl + 2) * (yl + 2), MPI_DOUBLE, rank + iproc * jproc, 0,
    MPI_COMM_WORLD);

    /* Receive the elements of the upper neighbor and set the down-facing
     * PDFs of the highest lattices of this rank to the buffer values.
     */
    MPI_Recv(&readBuffer[TO][0], 5 * (xl + 2) * (yl + 2), MPI_DOUBLE, rank + iproc * jproc, 0,
    MPI_COMM_WORLD, &status);
    i = 0;
    for (int x = 0; x < xl + 2; ++x) {
        for (int y = 0; y < yl + 2; ++y) {
            collideField[idx(sublength, x, y, zl + 1, 0)] = readBuffer[TO][i++];
            collideField[idx(sublength, x, y, zl + 1, 1)] = readBuffer[TO][i++];
            collideField[idx(sublength, x, y, zl + 1, 2)] = readBuffer[TO][i++];
            collideField[idx(sublength, x, y, zl + 1, 3)] = readBuffer[TO][i++];
            collideField[idx(sublength, x, y, zl + 1, 4)] = readBuffer[TO][i++];
        }
    }
}

void exchangePdfs(double **sendBuffer, double **readBuffer, double *collideField,
        const int * const sublength, int rank, int iproc, int jproc, int kproc)
{
    /* ====================================
     * ========= X direction ==============
     * ==================================== */
    /* If the current rank does NOT process a right boundary subdomain, we must exchange PDFs
     * with the next (right) neighbor. */
    if ((checkBoundary(rank, iproc, jproc, kproc) & RIGHT) == 0) {
        doRightSendRecv(sendBuffer, readBuffer, collideField, sublength, rank);
    }
    /* If the current rank does NOT process a left boundary subdomain, we must exchange PDFs
     * with the previous (left) neighbor. */
    if ((checkBoundary(rank, iproc, jproc, kproc) & LEFT) == 0) {
        doLeftSendRecv(sendBuffer, readBuffer, collideField, sublength, rank);
    }
    /* ====================================
     * ========= Y direction ==============
     * ==================================== */
    /* If the current rank does NOT process a back boundary subdomain, we must exchange PDFs
     * with the previous neighbor. */
    if ((checkBoundary(rank, iproc, jproc, kproc) & BACK) == 0) {
        doBackSendRecv(sendBuffer, readBuffer, collideField, sublength, rank, iproc);
    }
    /* If the current rank does NOT process a front boundary subdomain, we must exchange PDFs
     * with the neighbor behind. */
    if ((checkBoundary(rank, iproc, jproc, kproc) & FRONT) == 0) {
        doFrontSendRecv(sendBuffer, readBuffer, collideField, sublength, rank, iproc);
    }
    /* ====================================
     * ========= Z direction ==============
     * ==================================== */
    /* If the current rank does NOT process a top boundary subdomain, we must exchange
     * PDFs with the +xlength*xlength (upper) neighbor. */
    if ((checkBoundary(rank, iproc, jproc, kproc) & TOP) == 0) {
        doTopSendRecv(sendBuffer, readBuffer, collideField, sublength, rank, iproc, jproc);
    }
    /* If the current rank does NOT process a bottom boundary subdomain, we must exchange
     * PDFs with the -xlength*xlength (lower) neighbor. */
    if ((checkBoundary(rank, iproc, jproc, kproc) & BOTTOM) == 0) {
        doBottomSendRecv(sendBuffer, readBuffer, collideField, sublength, rank, iproc, jproc);
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

