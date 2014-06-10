#include "collision.h"
#include "LBDefinitions.h"

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
    for (int z = 1; z < sublength[0] + 1; ++z) {
        for (int y = 1; y < sublength[1] + 1; ++y) {
            for (int x = 1; x < sublength[2] + 1; ++x) {
                double *currentCell = &collideField[idx(sublength, x, y, z, 0)];
                /*Updating values for velocity, density and Feq for Current cell*/
                double density;
                computeDensity(currentCell, &density);
                double velocity[D];
                computeVelocity(currentCell, &density, velocity);
                double feq[Q];
                computeFeq(&density, velocity, feq);
                computePostCollisionDistributions(currentCell, tau, feq);
            }
        }
    }
}
