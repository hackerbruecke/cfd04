#include <stdio.h>
#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

/* Check whether (x,y,z) is inside the domain bounds, including the ghost layer */
static inline int inBounds(int x, int y, int z, const int * const sublength)
{
    return x > 0 && x < sublength[0] + 1 && y > 0 && y < sublength[1] + 1 && z > 0
            && z < sublength[2] + 1;
}

/* Computes the no-slip condition for cell (x,y,z)
 * if the respective neighbor determined by i is within bounds and a fluid cell */
static inline void computeNoSlipCondition(const int * const sublength, int x, int y, int z, int i,
        int dx, int dy, int dz, double *collideField, const int * const flagField)
{
    /* Check if the i-th directions is facing the fluid particles and update accordingly */
    if (inBounds(x + dx, y + dy, z + dz, sublength)
            && flagField[fidx(sublength, x + dx, z + dy, z + dz)] == FLUID) {
        collideField[idx(sublength, x, y, z, i)] = collideField[idx(sublength, x + dx, y + dy,
                z + dz, inv(i))];
    }
}

/* Computed the moving-wall condition for cell (x,y,z)
 * if the respective neighbor determined by i is within bounds and a fluid cell */
static inline void computeMovingWallCondition(const int * const sublength, int x, int y, int z,
        int i, int dx, int dy, int dz, double *collideField, const int * const flagField,
        const double * const wallVelocity)
{
    if (inBounds(x + dx, y + dy, z + dz, sublength)
            && flagField[fidx(sublength, x + dx, y + dy, z + dz)] == FLUID) {
        double *currentCell = &collideField[idx(sublength, x + dx, y + dy, z + dz, 0)];

        double density;
        computeDensity(currentCell, &density);
        double cu = 0; /*dot product of c_i and velocity of wall*/
        for (int d = 0; d < D; ++d) {
            cu += LATTICEVELOCITIES[i][d] * wallVelocity[d];
        }
        double finv = collideField[idx(sublength, x + dx, y + dy, z + dz, inv(i))];
        collideField[idx(sublength, x, y, z, i)] = finv
                + 2 * LATTICEWEIGHTS[i] * density * cu / (C_S * C_S);
    }
}

void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity,
        const int * const sublength)
{
    const int xmax = sublength[0]+1;
    const int ymax = sublength[1]+1;
    const int zmax = sublength[2]+1;
    /*Going through the whole boundary domain*/

    /*For each boundary cell/particle we update all the distribution functions,
     that point towards fluid, according to boundary conditions*/
    for (int i = 0; i < Q; ++i) {
        /*dx = c_i_x*dt, dt = 1*/
        int dx = LATTICEVELOCITIES[i][0];
        int dy = LATTICEVELOCITIES[i][1];
        int dz = LATTICEVELOCITIES[i][2];
        for (int x = 0; x < xmax+1; ++x) {
            for (int y = 0; y < ymax+1; ++y) {
                if (flagField[fidx(sublength, x, y, 0)] == NO_SLIP)
                    computeNoSlipCondition(sublength, x, y, 0, i, dx, dy, dz, collideField, flagField);
                if (flagField[fidx(sublength, x, y, zmax)] == MOVING_WALL)
                    computeMovingWallCondition(sublength, x, y, zmax, i, dx, dy, dz,
                            collideField, flagField, wallVelocity);
            }
        }

        for (int x = 0; x < xmax+1; ++x) {
            for (int z = 0; z < zmax+1; ++z) {
                if (flagField[fidx(sublength, x, 0, z)] == NO_SLIP) 
                    computeNoSlipCondition(sublength, x, 0, z, i, dx, dy, dz, collideField, flagField);
                if (flagField[fidx(sublength, x, ymax, z)] == NO_SLIP)
                    computeNoSlipCondition(sublength, x, ymax, z, i, dx, dy, dz,
                            collideField, flagField);
            }
        }

        for (int y = 0; y < ymax+1; ++y) {
            for (int z = 0; z < zmax+1; ++z) {
                if (flagField[fidx(sublength, 0, y, z)] == NO_SLIP)
                    computeNoSlipCondition(sublength, 0, y, z, i, dx, dy, dz, collideField, flagField);
                if (flagField[fidx(sublength, xmax, y, z)] == NO_SLIP)
                    computeNoSlipCondition(sublength, xmax, y, z, i, dx, dy, dz,
                            collideField, flagField);
            }
        }
    }
}
