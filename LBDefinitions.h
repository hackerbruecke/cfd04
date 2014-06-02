#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

/* Dimensions of D3Q19 model */
#define D 3
#define Q 19

/* Three-dimensional lattice velocities (c_i) sorted by index i */
static const int LATTICEVELOCITIES[19][3] =
{
    {  0, -1, -1 }, { -1,  0, -1 }, {  0,  0, -1 }, {  1,  0, -1 }, {  0,  1, -1 },
    { -1, -1,  0 }, {  0, -1,  0 }, {  1, -1,  0 }, { -1,  0,  0 }, {  0,  0,  0 },
    {  1,  0,  0 }, { -1,  1,  0 }, {  0,  1,  0 }, {  1,  1,  0 }, {  0, -1,  1 },
    { -1,  0,  1 }, {  0,  0,  1 }, {  1,  0,  1 }, {  0,  1,  1 }
};

/* Lattice weights sorted by directions c_i */
static const double LATTICEWEIGHTS[19] =
{
    1.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36,
    1.0/36, 2.0/36, 1.0/36, 2.0/36, 12.0/36,
    2.0/36, 1.0/36, 2.0/36, 1.0/36, 1.0/36,
    1.0/36, 2.0/36, 1.0/36, 1.0/36
};

/* Speed of sound */
static const double C_S = 0.57735026919l; 

/* Type/State of lattice element */
typedef enum {
	FLUID = 0,
	NO_SLIP = 1,
	MOVING_WALL = 2
} STATE;

/*
 * Cell index
 */
static inline int idx(int* xlength, int x, int y, int z, int i) {
    return Q * (z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x) + i;
}
/*
 * Flag index
 */
static inline int fidx(int* xlength, int x, int y, int z) {
    return z * (xlength[0] + 2) * (xlength[1] + 2) + y * (xlength[0] + 2) + x;
}

/*
 * Inverse index
 */
static inline int inv(int idx) {
    return 18 - idx;
}
#endif

