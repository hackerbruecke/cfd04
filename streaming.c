#include "streaming.h"
#include "LBDefinitions.h"

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

