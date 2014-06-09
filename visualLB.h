#ifndef _VISUALLB_H_
#define _VISUALLB_H_

void writeParallelVtkFile(const char * const filename, int iproc, int jproc, int kproc,
        const int * const sublength, unsigned int t);

/** writes the density and velocity field (derived from the distributions in collideField)
 *  to a file determined by 'filename' and timestep 't'. You can re-use parts of the code
 *  from visual.c (VTK output for Navier-Stokes solver) and modify it for 3D datasets.
 */
void writeVtkOutput(const double * const collideField, const int * const flagField,
        const char * filename, int rank, int iproc, int jproc, int kproc, unsigned int t,
        const int * const xlength);

#endif

