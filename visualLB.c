#include "visualLB.h"
#include <stdio.h>

#include "LBDefinitions.h"
#include "helper.h"
#include "computeCellValues.h"

void writeParallelVtkFile(const char * const filename, int xlength, int iproc, int jproc, int kproc,
        const int * const sublength, unsigned int t)
{
    char fn[80];
    char sourceFile[256];
    FILE* fp = NULL;
    sprintf(fn, "%s_par.%d.pvts", filename, t);
    fp = fopen(fn, "w");
    if (fp == NULL) {
        ERROR("Failed to open file!");
        return;
    }
    /* write header */
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "   <PStructuredGrid WholeExtent=\"0 %d 0 %d 0 %d \" GhostLevel=\"#\">\n",
            xlength-1, xlength-1 , xlength-1);

    fprintf(fp, "       <PPointData Scalars=\"density\" Vectors=\"velocity\">\n");
    fprintf(fp, "           <DataArray type=\"Float32\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\" />\n");
    fprintf(fp, "           <DataArray type=\"Float32\" Name=\"density\" />\n");
    fprintf(fp, "       </PPointData>\n");

    fprintf(fp, "\t\t<PCellData></PCellData>\n");

    fprintf(fp, "\t\t<PPoints>\n");
    fprintf(fp, "\t\t\t<DataArray  type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"ascii\" />\n");
    fprintf(fp, "\t\t</PPoints>\n");

    for (int k = 0; k < kproc; ++k) {
        for (int j = 0; j < jproc; ++j) {
            for (int i = 0; i < iproc; ++i) {

                int xa = i * sublength[0];
                xa > 1 ? --xa : xa;
                int xb = (i + 1) * sublength[0]-1;

                int ya = j * sublength[1];
                ya > 1 ? --ya : ya;
                int yb = (j + 1) * sublength[1]-1;

                int za = k * sublength[2];
                za > 1 ? --za : za;
                int zb = (k + 1) * sublength[2]-1;

                sprintf(sourceFile, "%s_%d.%d.vts", filename, i + iproc * j + iproc * jproc * k, t);
                fprintf(fp, "       <Piece Extent=\"%d %d %d %d %d %d\" Source=\"%s\">\n",
                        xa, xb, ya, yb, za, zb, sourceFile);
                fprintf(fp, "\t\t</Piece>\n");
            }
        }
    }

    fprintf(fp, "    </PStructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");
}

void writeVtkOutput(const double * const collideField, const int * const flagField,
        const char * filename, int rank, int iproc, int jproc, int kproc, unsigned int t,
        const int * const sublength)
{
    char fn[256];
    /* Save filename as a combination of passed filename and timestep */
    sprintf(fn, "%s_%d.%d.vts", filename, rank, t);

    FILE *fp = fopen(fn, "w");
    if (fp == NULL) {
        ERROR("Failed to open file!");
        return;
    }

    const int xl = sublength[0];
    const int yl = sublength[1];
    const int zl = sublength[2];
    fprintf(fp, "<?xml version=\"1.0\"?>\n");
    fprintf(fp, "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(fp, "   <StructuredGrid WholeExtent=\"0 %d 0 %d 0 %d\" >\n", xl-1, yl-1, zl-1);
    fprintf(fp, "   <Piece Extent=\"0 %d 0 %d 0 %d\" >\n", xl - 1, yl - 1, zl - 1);
    fprintf(fp, "       <PointData Scalars=\"density\" Vectors=\"velocity\">\n");
    /* ====================== Velocity ===================== */
    fprintf(fp,
            "           <DataArray Name=\"velocity\" type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n");

    double density;
    double vel[D];
    const double *currentCell;
    /* Compute (macroscopic) velocities for all cells */
    for (int z = 1; z < zl + 1; ++z) {
        for (int y = 1; y < yl + 1; ++y) {
            for (int x = 1; x < xl + 1; ++x) {
                currentCell = &collideField[idx(sublength, x, y, z, 0)];
                computeDensity(currentCell, &density);
                computeVelocity(currentCell, &density, vel);
                fprintf(fp, "%f %f %f\n", vel[0], vel[1], vel[2]);
            }
        }
    }
    fprintf(fp, "           </DataArray>\n");
    /* ====================== Density ===================== */
    fprintf(fp, "           <DataArray Name=\"density\" type=\"Float32\">\n");

    /* Compute density for each cell */
    for (int z = 1; z < zl + 1; ++z) {
        for (int y = 1; y < yl + 1; ++y) {
            for (int x = 1; x < xl + 1; ++x) {
                currentCell = &collideField[idx(sublength, x, y, z, 0)];
                computeDensity(currentCell, &density);
                fprintf(fp, "%f\n", density);
            }
        }
    }

    fprintf(fp, "           </DataArray>\n");
    fprintf(fp, "       </PointData>\n");
    fprintf(fp, "       <CellData></CellData>\n");
    fprintf(fp, "       <Points>\n");
    /* ====================== Lattice points ===================== */
    fprintf(fp,
            "           <DataArray type=\"Float32\" Name=\"coordinate\" NumberOfComponents=\"3\" format=\"ascii\">\n");
    /* Print lattice points */
    for (int z = 1; z < zl + 1; ++z) {
        for (int y = 1; y < yl + 1; ++y) {
            for (int x = 1; x < xl + 1; ++x) {
                fprintf(fp, "%d %d %d\n", x, y, z);
            }
        }
    }
    fprintf(fp, "           </DataArray>\n");
    fprintf(fp, "       </Points>\n");
    fprintf(fp, "    </Piece>\n");
    fprintf(fp, "    </StructuredGrid>\n");
    fprintf(fp, "</VTKFile>\n");

    if (fclose(fp)) {
        ERROR("Failed to close file!");
        return;
    }
}

