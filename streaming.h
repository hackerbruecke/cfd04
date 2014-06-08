#ifndef _STREAMING_H_
#define _STREAMING_H_

void exchangePdfs(double **sendBuffer, double **readBuffer, double *collideField,
        const int * const sublength, int rank, int iproc, int jproc, int kproc);

/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
void doStreaming(double *collideField, double *streamField, int *flagField,
        const int * const sublength);

#endif

