#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_

/* Dimensions of D3Q19 model */
#define D 3
#define Q 19

/* Three-dimensional lattice velocities (c_i) sorted by index i */
static const int LATTICEVELOCITIES[19][3] = { { 0, -1, -1 }, { -1, 0, -1 }, { 0, 0, -1 },
        { 1, 0, -1 }, { 0, 1, -1 }, { -1, -1, 0 }, { 0, -1, 0 }, { 1, -1, 0 }, { -1, 0, 0 }, { 0, 0,
                0 }, { 1, 0, 0 }, { -1, 1, 0 }, { 0, 1, 0 }, { 1, 1, 0 }, { 0, -1, 1 },
        { -1, 0, 1 }, { 0, 0, 1 }, { 1, 0, 1 }, { 0, 1, 1 } };

/* Lattice weights sorted by directions c_i */
static const double LATTICEWEIGHTS[19] = { 1.0 / 36, 1.0 / 36, 2.0 / 36, 1.0 / 36, 1.0 / 36, 1.0
        / 36, 2.0 / 36, 1.0 / 36, 2.0 / 36, 12.0 / 36, 2.0 / 36, 1.0 / 36, 2.0 / 36, 1.0 / 36, 1.0
        / 36, 1.0 / 36, 2.0 / 36, 1.0 / 36, 1.0 / 36 };

/* Speed of sound */
static const double C_S = 0.57735026919l;

/* Type/State of lattice element */
typedef enum
{
    FLUID = 0, NO_SLIP = 1, MOVING_WALL = 2, PARALLEL_BOUNDARY
} STATE;

#if 1 /* Maybe remove this later on */
typedef enum
{
    NONE    = 0,
    LEFT    = 1 << 0,
    RIGHT   = 1 << 1,
    TOP     = 1 << 2,
    BOTTOM  = 1 << 3,
    FRONT   = 1 << 4,
    BACK    = 1 << 5
} BOUNDARY;

static inline int checkBoundary(int rank, int iproc, int jproc, int kproc)
{
    int result = NONE;
    if (rank / (iproc * jproc) == 0)
        result |= BOTTOM;
    if (rank / (iproc * jproc) == kproc-1)
        result |= TOP;
    if ((rank / iproc) % jproc == 0)
        result |= FRONT;
    if ((rank / iproc) % jproc == jproc-1)
        result |= BACK;
    if (rank % iproc == 0)
        result |= LEFT;
    if (rank % iproc == iproc-1)
        result |= RIGHT;
    return result;
}
#endif

/*
 * Cell index
 */
static inline int idx(const int * const sublength, int x, int y, int z, int i)
{
    return Q * (z * (sublength[1] + 2) * (sublength[0] + 2) + y * (sublength[0] + 2) + x) + i;
}
/*
 * Flag index
 */
static inline int fidx(const int * const sublength, int x, int y, int z)
{
    return z * (sublength[1] + 2) * (sublength[0] + 2) + y * (sublength[0] + 2) + x;
}

/*
 * Subdomain index. Returns the rank corresponding to the x,y,z coordinates of a subdomain element
 */
static inline int sidx(const int * const xlength, int x, int y, int z)
{
    return x + xlength[0] * y + xlength[0] * xlength[1] * z;
}

/*
 * Inverse index
 */
static inline int inv(int idx)
{
    return 18 - idx;
}
#endif

