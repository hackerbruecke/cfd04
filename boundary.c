#include <stdio.h>
#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"

/* Check whether (x,y,z) is inside the domain bounds, including the ghost layer */
static inline int inBounds(int x, int y, int z, int xlength) {
	return  x > 0 && x < xlength+1 && y > 0 && y < xlength+1 && z > 0 && z < xlength+1;
}

/* Computes the no-slip condition for cell (x,y,z)
 * if the respective neighbor determined by i is within bounds and a fluid cell */
static inline void computeNoSlipCondition(int xlength, int x, int y, int z, int i, int dx, int dy, int dz,
		double *collideField, const int * const flagField) {
	/* Check if the i-th directions is facing the fluid particles and update accordingly */
	if (inBounds(x+dx, y+dy, z+dz, xlength) && flagField[fidx(xlength, x+dx, z+dy, z+dz)] == FLUID) {
		collideField[idx(xlength, x, y, z, i)] = collideField[idx(xlength, x + dx, y + dy, z + dz, inv(i))];
	}
}

/* Computed the moving-wall condition for cell (x,y,z)
 * if the respective neighbor determined by i is within bounds and a fluid cell */
static inline void computeMovingWallCondition(int xlength, int x, int y, int z, int i, int dx, int dy, int dz,
		double *collideField, const int * const flagField, const double * const wallVelocity) {
	if (inBounds(x+dx, y+dy, z+dz, xlength) && flagField[fidx(xlength, x+dx, y+dy, z+dz)] == FLUID) {
		double *currentCell = &collideField[idx(xlength, x+dx, y+dy, z+dz, 0)];

		double density;
		computeDensity(currentCell, &density);
		double cu = 0; /*dot product of c_i and velocity of wall*/
		for (int d = 0; d < D; ++d) {
			cu += LATTICEVELOCITIES[i][d] * wallVelocity[d];
		}
		double finv = collideField[idx(xlength, x + dx, y + dy, z + dz, inv(i))];
		collideField[idx(xlength, x, y, z, i)] = finv + 2*LATTICEWEIGHTS[i]*density*cu/(C_S*C_S);
	}
}

void treatBoundary(double *collideField, int* flagField,
		const double * const wallVelocity, int xlength) {
	/*Going through the whole boundary domain*/

	/*For each boundary cell/particle we update all the distribution functions,
	 that point towards fluid, according to boundary conditions*/
	for (int i = 0; i < Q; ++i) {
		/*dx = c_i_x*dt, dt = 1*/
		int dx = LATTICEVELOCITIES[i][0];
		int dy = LATTICEVELOCITIES[i][1];
		int dz = LATTICEVELOCITIES[i][2];
		for (int x = 0; x < xlength + 2; ++x) {
			for (int y = 0; y < xlength + 2; ++y) {
				computeNoSlipCondition(xlength, x, y, 0, i, dx, dy, dz, collideField, flagField);
				computeMovingWallCondition(xlength, x, y, xlength+1, i, dx, dy, dz, collideField, flagField, wallVelocity);
			}
		}

		for (int x = 0; x < xlength + 2; ++x) {
			for (int z = 0; z < xlength + 2; ++z) {
				computeNoSlipCondition(xlength, x, 0, z, i, dx, dy, dz, collideField, flagField);
				computeNoSlipCondition(xlength, x, xlength+1, z, i, dx, dy, dz, collideField, flagField);
			}
		}

		for (int y = 0; y < xlength + 2; ++y) {
			for (int z = 0; z < xlength + 2; ++z) {
				computeNoSlipCondition(xlength, 0, y, z, i, dx, dy, dz, collideField, flagField);
				computeNoSlipCondition(xlength, xlength+1, y, z, i, dx, dy, dz, collideField, flagField);
			}
		}
	}
}
