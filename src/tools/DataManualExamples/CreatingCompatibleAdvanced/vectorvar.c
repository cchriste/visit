#include <silo.h>
#include <stdio.h>

void
write_curv3d(DBfile *dbfile)
{
    /* Write a curvilinear mesh. */
    float x[2][3][4] = {
       {{0.,1.,2.,3.},{0.,1.,2.,3.}, {0.,1.,2.,3.}},
       {{0.,1.,2.,3.},{0.,1.,2.,3.}, {0.,1.,2.,3.}}
    };
    float y[2][3][4] = {
       {{0.5,0.,0.,0.5},{1.,1.,1.,1.}, {1.5,2.,2.,1.5}},
       {{0.5,0.,0.,0.5},{1.,1.,1.,1.}, {1.5,2.,2.,1.5}}
    };
    float z[2][3][4] = {
       {{0.,0.,0.,0.},{0.,0.,0.,0.},{0.,0.,0.,0.}},
       {{1.,1.,1.,1.},{1.,1.,1.,1.},{1.,1.,1.,1.}}
    };
    int dims[] = {4, 3, 2};
    int ndims = 3;
    float *coords[] = {(float*)x, (float*)y, (float*)z};
    DBPutQuadmesh(dbfile, "quadmesh", NULL, coords, dims, ndims,
                  DB_FLOAT, DB_NONCOLLINEAR, NULL);
}

#define NX 4
#define NY 3
#define NZ 2

void
write_zonecent_quadvar(DBfile *dbfile)
{
    int i, dims[3], ndims = 3;
    int ncells = (NX-1)*(NY-1)*(NZ-1);
    float *comp[3];
    char *varnames[] = {"zonal_comp0","zonal_comp1","zonal_comp2"};
    comp[0] = (float *)malloc(sizeof(float)*ncells);
    comp[1] = (float *)malloc(sizeof(float)*ncells);
    comp[2] = (float *)malloc(sizeof(float)*ncells);
    for(i = 0; i < ncells; ++i)
    {
        comp[0][i] = (float)i; /*vector component 0*/
        comp[1][i] = (float)i; /*vector component 1*/
        comp[2][i] = (float)i; /*vector component 2*/
    }
    dims[0] = NX-1; dims[1] = NY-1; dims[2] = NZ-1;
    DBPutQuadvar(dbfile, "zonal", "quadmesh",
                 3, varnames, comp, dims,
                 ndims, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);
    free(comp[0]);
    free(comp[1]);
    free(comp[2]);
}

void
write_nodecent_quadvar(DBfile *dbfile)
{
    int i, dims[3], ndims = 3;
    int nnodes = NX*NY*NZ;
    float *comp[3];
    char *varnames[] = {"nodal_comp0","nodal_comp1","nodal_comp2"};
    comp[0] = (float *)malloc(sizeof(float)*nnodes);
    comp[1] = (float *)malloc(sizeof(float)*nnodes);
    comp[2] = (float *)malloc(sizeof(float)*nnodes);
    for(i = 0; i < nnodes; ++i)
    {
        comp[0][i] = (float)i; /*vector component 0*/
        comp[1][i] = (float)i; /*vector component 1*/
        comp[2][i] = (float)i; /*vector component 2*/
    }
    dims[0] = NX; dims[1] = NY; dims[2] = NZ;
    DBPutQuadvar(dbfile, "nodal", "quadmesh",
                 3, varnames, comp, dims,
                 ndims, NULL, 0, DB_FLOAT, DB_NODECENT, NULL);
    free(comp[0]);
    free(comp[1]);
    free(comp[2]);
}

int
main(int argc, char *argv[])
{
    DBfile *dbfile = NULL;
    /* Open the Silo file */
    dbfile = DBCreate("vectorvar.silo", DB_CLOBBER, DB_LOCAL,
        "Vector quadvars in 3D", DB_HDF5);
    if(dbfile == NULL)
    {
        fprintf(stderr, "Could not create Silo file!\n");
        return -1;
    }

    /* Add other Silo calls here. */
    write_curv3d(dbfile);
    write_zonecent_quadvar(dbfile);
    write_nodecent_quadvar(dbfile);

    /* Close the Silo file. */
    DBClose(dbfile);
    return 0;
}
