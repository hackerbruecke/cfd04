#include <string.h>

#include "LBDefinitions.h"
#include "computeCellValues.h"

void computeDensity(const double * const currentCell, double *density)
{
    /* Sum of all fi in the Cell to compute density */
    *density = 0;
    for (int i = 0; i < Q; ++i) {
        *density += currentCell[i];
    }
}

void computeVelocity(const double * const currentCell, const double * const density,
        double *velocity)
{
    /* Initializing velocity vector to 0 */
    memset(velocity, 0, sizeof *velocity * D);

    /* Compute velocity by momentum equation */
    for (int d = 0; d < D; ++d) {
        for (int q = 0; q < Q; ++q) {
            velocity[d] += currentCell[q] * LATTICEVELOCITIES[q][d];
        }
    }
    for (int d = 0; d < D; ++d) {
        velocity[d] /= *density;
    }
}

void computeFeq(const double * const density, const double * const velocity, double *feq)
{
    double c_dot_u;
    double u_dot_u;

    for (int q = 0; q < Q; ++q) {
        c_dot_u = 0.0;
        u_dot_u = 0.0;
        for (int d = 0; d < D; ++d) {
            /* Compute dot products needed for the computation of feq */
            u_dot_u += velocity[d] * velocity[d];
            c_dot_u += LATTICEVELOCITIES[q][d] * velocity[d];
        }
        feq[q] = LATTICEWEIGHTS[q] * (*density)
                * (1 + (c_dot_u) / (C_S * C_S) + (c_dot_u) * (c_dot_u) / (2 * C_S * C_S * C_S * C_S)
                        - u_dot_u / (2 * C_S * C_S));
    }
}

