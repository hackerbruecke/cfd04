#include "collision.h"
#include "LBDefinitions.h"
#include <stdio.h>
#include <math.h>

void computePostCollisionDistributions(double *currentCell, const double * const tau,
        const double * const feq)
{
    /* Updating all(19) post collision distribution fields using BGK approx rule in Current Cell */
    for (int i = 0; i < Q; ++i) {
        currentCell[i] -= (currentCell[i] - feq[i]) / *tau;
    }
}

void doCollision(double *collideField, int *flagField, const double * const tau,
        const int * const sublength)
{
    double density;
    double velocity[D];
    double feq[Q];
    double *currentCell;

    for (int z = 1; z < sublength[0] + 1; ++z) {
        for (int y = 1; y < sublength[1] + 1; ++y) {
            for (int x = 1; x < sublength[2] + 1; ++x) {
                currentCell = &collideField[idx(sublength, x, y, z, 0)];
                /*Updating values for velocity, density and Feq for Current cell*/
                computeDensity(currentCell, &density);
                computeVelocity(currentCell, &density, velocity);
                computeFeq(&density, velocity, feq);
                computePostCollisionDistributions(currentCell, tau, feq);
            }
        }
    }
}
