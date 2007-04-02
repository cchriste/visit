/* SIMPLE PARALLEL SIMULATION SKELETON */
#include <VisItControlInterface_V1.h>
#include <stdio.h>
#include <string.h>

void simulate_one_timestep(void);
void read_input_deck(void) { }
int  simulation_done(void)   { return 0; }

/* Is the simulation in run mode (not waiting for VisIt input) */
static int    runFlag = 1;
static int    simcycle = 0;
static double simtime = 0.;

/* Callback function for control commands. */
void ControlCommandCallback(const char *cmd,
    int int_data, float float_data,
    const char *string_data)
{
    if(strcmp(cmd, "halt") == 0)
        runFlag = 0;
    else if(strcmp(cmd, "step") == 0)
        simulate_one_timestep();
    else if(strcmp(cmd, "run") == 0)
        runFlag = 1;
}

void mainloop(void)
{
    int blocking, visitstate, err = 0;

    do
    {
        blocking = runFlag ? 0 : 1;
        /* Get input from VisIt or timeout so the simulation can run. */
        visitstate = VisItDetectInput(blocking, -1);

        /* Do different things depending on the output from VisItDetectInput. */
        if(visitstate >= -5 && visitstate <= -1)
        {
            fprintf(stderr, "Can't recover from error!\n");
            err = 1;
        }
        else if(visitstate == 0)
        {
            /* There was no input from VisIt, return control to sim. */
            simulate_one_timestep();
        }
        else if(visitstate == 1)
        {
            /* VisIt is trying to connect to sim. */
            if(VisItAttemptToCompleteConnection())
            {
                fprintf(stderr, "VisIt connected\n");
                VisItSetCommandCallback(ControlCommandCallback);
            }
            else
                fprintf(stderr, "VisIt did not connect\n");
        }
        else if(visitstate == 2)
        {
            /* VisIt wants to tell the engine something. */
            runFlag = 0;
            if(!VisItProcessEngineCommand())
            {
                /* Disconnect on an error or closed connection. */
                VisItDisconnect();
                /* Start running again if VisIt closes. */
                runFlag = 1;
            }
        }
    } while(!simulation_done() && err == 0);
}

int main(int argc, char **argv)
{
    /* Initialize environment variables. */
    VisItSetupEnvironment();
    /* Write out .sim file that VisIt uses to connect. */
    VisItInitializeSocketAndDumpSimFile("scalar",
        "Demonstrates scalar data access function",
        "/path/to/where/sim/was/started",
        NULL, NULL);

    /* Read input problem setup, geometry, data. */
    read_input_deck();

    /* Call the main loop. */
    mainloop();

    return 0;
}

/* SIMULATE ONE TIME STEP */
#include <unistd.h>
void simulate_one_timestep(void)
{
    ++simcycle;
    simtime += 0.0134;
    printf("Simulating time step: cycle=%d, time=%lg\n", simcycle, simtime);
    sleep(1);
}

/* DATA ACCESS FUNCTIONS */
#include <VisItDataInterface_V1.h>

VisIt_SimulationMetaData *VisItGetMetaData(void)
{
    /* Create a metadata object with no variables. */
    size_t sz = sizeof(VisIt_SimulationMetaData);
    VisIt_SimulationMetaData *md = 
        (VisIt_SimulationMetaData *)malloc(sz);
    memset(md, 0, sz);

    /* Set the simulation state. */
    md->currentMode = runFlag ? VISIT_SIMMODE_RUNNING : VISIT_SIMMODE_STOPPED;
    md->currentCycle = simcycle;
    md->currentTime = simtime;

#define NDOMAINS 1
    /* Allocate enough room for 2 meshes in the metadata. */
    md->numMeshes = 2;
    sz = sizeof(VisIt_MeshMetaData) * md->numMeshes;
    md->meshes = (VisIt_MeshMetaData *)malloc(sz);
    memset(md->meshes, 0, sz);

    /* Set the first mesh's properties.*/
    md->meshes[0].name = strdup("mesh2d");
    md->meshes[0].meshType = VISIT_MESHTYPE_RECTILINEAR;
    md->meshes[0].topologicalDimension = 2;
    md->meshes[0].spatialDimension = 2;
    md->meshes[0].numBlocks = NDOMAINS;
    md->meshes[0].blockTitle = strdup("Domains");
    md->meshes[0].blockPieceName = strdup("domain");
    md->meshes[0].numGroups = 0;
    md->meshes[0].units = strdup("cm");
    md->meshes[0].xLabel = strdup("Width");
    md->meshes[0].yLabel = strdup("Height");
    md->meshes[0].zLabel = strdup("Depth");

    /* Set the second mesh's properties.*/
    md->meshes[1].name = strdup("mesh3d");
    md->meshes[1].meshType = VISIT_MESHTYPE_CURVILINEAR;
    md->meshes[1].topologicalDimension = 3;
    md->meshes[1].spatialDimension = 3;
    md->meshes[1].numBlocks = NDOMAINS;
    md->meshes[1].blockTitle = strdup("Domains");
    md->meshes[1].blockPieceName = strdup("domain");
    md->meshes[1].numGroups = 0;
    md->meshes[1].units = strdup("Miles");
    md->meshes[1].xLabel = strdup("Width");
    md->meshes[1].yLabel = strdup("Height");
    md->meshes[1].zLabel = strdup("Depth");

    /* Add some scalar variables. */
    md->numScalars = 2;
    sz = sizeof(VisIt_ScalarMetaData) * md->numScalars;
    md->scalars = (VisIt_ScalarMetaData *)malloc(sz);
    memset(md->scalars, 0, sz);

    /* Add a zonal variable on mesh2d. */
    md->scalars[0].name = strdup("zonal");
    md->scalars[0].meshName = strdup("mesh2d");
    md->scalars[0].centering = VISIT_VARCENTERING_ZONE;

    /* Add a nodal variable on mesh3d. */
    md->scalars[1].name = strdup("nodal");
    md->scalars[1].meshName = strdup("mesh3d");
    md->scalars[1].centering = VISIT_VARCENTERING_NODE;

    /* Add a curve variable. */
    md->numCurves = 1;
    sz = sizeof(VisIt_CurveMetaData) * md->numCurves;
    md->curves = (VisIt_CurveMetaData *)malloc(sz);
    memset(md->curves, 0, sz);

    md->curves[0].name = strdup("sine");
    md->curves[0].xUnits = strdup("radians");
    md->curves[0].xLabel = strdup("angle");
    md->curves[0].yLabel = strdup("amplitude");

    /* Add a material variable. */
    md->numMaterials = 1;
    sz = sizeof(VisIt_MaterialMetaData) * md->numMaterials;
    md->materials = (VisIt_MaterialMetaData *)malloc(sz);
    memset(md->materials, 0, sz);

    md->materials[0].name = strdup("mat");
    md->materials[0].meshName = strdup("mesh2d");
    md->materials[0].numMaterials = 3;
    /* Allocate memory to store the list of material names. */
    md->materials[0].materialNames = (const char **)malloc(sizeof(char *) *
        md->materials[0].numMaterials);
    md->materials[0].materialNames[0] = strdup("Iron");
    md->materials[0].materialNames[1] = strdup("Copper");
    md->materials[0].materialNames[2] = strdup("Nickel");

    /* Add some expressions. */
    md->numExpressions = 2;
    sz = sizeof(VisIt_ExpressionMetaData) * md->numExpressions;
    md->expressions = (VisIt_ExpressionMetaData *)malloc(sz);
    memset(md->expressions, 0, sz);

    md->expressions[0].name = strdup("zvec");
    md->expressions[0].definition = strdup("{zonal, zonal, zonal}");
    md->expressions[0].vartype = VISIT_VARTYPE_VECTOR;

    md->expressions[1].name = strdup("nid");
    md->expressions[1].definition = strdup("nodeid(mesh3d)");
    md->expressions[1].vartype = VISIT_VARTYPE_SCALAR;

    /* Add some custom commands. */
    md->numGenericCommands = 3;
    sz = sizeof(VisIt_SimulationControlCommand) * md->numGenericCommands;
    md->genericCommands = (VisIt_SimulationControlCommand *)malloc(sz);
    memset(md->genericCommands, 0, sz);

    md->genericCommands[0].name = strdup("halt");
    md->genericCommands[0].argType = VISIT_CMDARG_NONE;
    md->genericCommands[0].enabled = 1;

    md->genericCommands[1].name = strdup("step");
    md->genericCommands[1].argType = VISIT_CMDARG_NONE;
    md->genericCommands[1].enabled = 1;

    md->genericCommands[2].name = strdup("run");
    md->genericCommands[2].argType = VISIT_CMDARG_NONE;
    md->genericCommands[2].enabled = 1;

    return md;
}

/* Rectilinear mesh */
float rmesh_x[] = {0., 1., 2.5, 5.};
float rmesh_y[] = {0., 2., 2.25, 2.55,  5.};
int   rmesh_dims[] = {4, 5, 1};
int   rmesh_ndims = 2;
float zonal[] = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.};

/* Curvilinear mesh */
float cmesh_x[2][3][4] = {
   {{0.,1.,2.,3.},{0.,1.,2.,3.}, {0.,1.,2.,3.}},
   {{0.,1.,2.,3.},{0.,1.,2.,3.}, {0.,1.,2.,3.}}
};
float cmesh_y[2][3][4] = {
   {{0.5,0.,0.,0.5},{1.,1.,1.,1.}, {1.5,2.,2.,1.5}},
   {{0.5,0.,0.,0.5},{1.,1.,1.,1.}, {1.5,2.,2.,1.5}}
};
float cmesh_z[2][3][4] = {
   {{0.,0.,0.,0.},{0.,0.,0.,0.},{0.,0.,0.,0.}},
   {{1.,1.,1.,1.},{1.,1.,1.,1.},{1.,1.,1.,1.}}
};
int cmesh_dims[] = {4, 3, 2};
int cmesh_ndims = 3;
double nodal[2][3][4] = {
   {{1.,2.,3.,4.},{5.,6.,7.,8.},{9.,10.,11.,12}},
   {{13.,14.,15.,16.},{17.,18.,19.,20.},{21.,22.,23.,24.}}
};

VisIt_MeshData *VisItGetMesh(int domain, const char *name)
{
    VisIt_MeshData *mesh = NULL;
    size_t sz = sizeof(VisIt_MeshData);

    if(strcmp(name, "mesh2d") == 0)
    {
        /* Allocate VisIt_MeshData. */
        mesh = (VisIt_MeshData *)malloc(sz);
        memset(mesh, 0, sz);
        /* Make VisIt_MeshData contain a VisIt_RectilinearMesh. */
        sz = sizeof(VisIt_RectilinearMesh);
        mesh->rmesh = (VisIt_RectilinearMesh *)malloc(sz);
        memset(mesh->rmesh, 0, sz);

        /* Tell VisIt which mesh object to use. */
        mesh->meshType = VISIT_MESHTYPE_RECTILINEAR;

        /* Set the mesh's number of dimensions. */
        mesh->rmesh->ndims = rmesh_ndims;

        /* Set the mesh dimensions. */
        mesh->rmesh->dims[0] = rmesh_dims[0];
        mesh->rmesh->dims[1] = rmesh_dims[1];
        mesh->rmesh->dims[2] = rmesh_dims[2];

        mesh->rmesh->baseIndex[0] = 0;
        mesh->rmesh->baseIndex[1] = 0;
        mesh->rmesh->baseIndex[2] = 0;

        mesh->rmesh->minRealIndex[0] = 0;
        mesh->rmesh->minRealIndex[1] = 0;
        mesh->rmesh->minRealIndex[2] = 0;
        mesh->rmesh->maxRealIndex[0] = rmesh_dims[0]-1;
        mesh->rmesh->maxRealIndex[1] = rmesh_dims[1]-1;
        mesh->rmesh->maxRealIndex[2] = rmesh_dims[2]-1;

        /* Let VisIt use the simulation's copy of the mesh coordinates. */
        mesh->rmesh->xcoords = VisIt_CreateDataArrayFromFloat(
           VISIT_OWNER_SIM, rmesh_x);
        mesh->rmesh->ycoords = VisIt_CreateDataArrayFromFloat(
           VISIT_OWNER_SIM, rmesh_y);
    }
    else if(strcmp(name, "mesh3d") == 0)
    {
        /* Allocate VisIt_MeshData. */
        mesh = (VisIt_MeshData *)malloc(sz);
        memset(mesh, 0, sz);
        /* Make VisIt_MeshData contain a VisIt_CurvilinearMesh. */
        sz = sizeof(VisIt_CurvilinearMesh);
        mesh->cmesh = (VisIt_CurvilinearMesh *)malloc(sz);
        memset(mesh->cmesh, 0, sz);

        /* Tell VisIt which mesh object to use. */
        mesh->meshType = VISIT_MESHTYPE_CURVILINEAR;

        /* Set the mesh's number of dimensions. */
        mesh->cmesh->ndims = cmesh_ndims;

        /* Set the mesh dimensions. */
        mesh->cmesh->dims[0] = cmesh_dims[0];
        mesh->cmesh->dims[1] = cmesh_dims[1];
        mesh->cmesh->dims[2] = cmesh_dims[2];

        mesh->cmesh->baseIndex[0] = 0;
        mesh->cmesh->baseIndex[1] = 0;
        mesh->cmesh->baseIndex[2] = 0;

        mesh->cmesh->minRealIndex[0] = 0;
        mesh->cmesh->minRealIndex[1] = 0;
        mesh->cmesh->minRealIndex[2] = 0;
        mesh->cmesh->maxRealIndex[0] = cmesh_dims[0]-1;
        mesh->cmesh->maxRealIndex[1] = cmesh_dims[1]-1;
        mesh->cmesh->maxRealIndex[2] = cmesh_dims[2]-1;

        /* Let VisIt use the simulation's copy of the mesh coordinates. */
        mesh->cmesh->xcoords = VisIt_CreateDataArrayFromFloat(
           VISIT_OWNER_SIM, (float *)cmesh_x);
        mesh->cmesh->ycoords = VisIt_CreateDataArrayFromFloat(
           VISIT_OWNER_SIM, (float *)cmesh_y);
        mesh->cmesh->zcoords = VisIt_CreateDataArrayFromFloat(
           VISIT_OWNER_SIM, (float *)cmesh_z);
    }

    return mesh;
}

VisIt_ScalarData *VisItGetScalar(int domain, const char *name)
{
    size_t sz = sizeof(VisIt_ScalarData);
    VisIt_ScalarData *scalar = (VisIt_ScalarData*)malloc(sz);
    memset(scalar, 0, sz);

    if(strcmp(name, "zonal") == 0)
    {
        scalar->len = (rmesh_dims[0]-1) * (rmesh_dims[1]-1);
        scalar->data = VisIt_CreateDataArrayFromFloat(
            VISIT_OWNER_SIM, zonal);
    }
    else if(strcmp(name, "nodal") == 0)
    {
        scalar->len = cmesh_dims[0] * cmesh_dims[1] *
            cmesh_dims[2];
        scalar->data = VisIt_CreateDataArrayFromDouble(
            VISIT_OWNER_SIM, (double*)nodal);
    }
    else 
    {
        free(scalar);
        scalar = NULL;
    }

    return scalar;
}

VisIt_SimulationCallback visitCallbacks =
{
    &VisItGetMetaData,
    &VisItGetMesh,
    NULL, /* GetMaterial */
    NULL, /* GetSpecies */
    &VisItGetScalar, /* GetScalar */
    NULL, /* GetCurve */
    NULL, /* GetMixedScalar */
    NULL /* GetDomainList */
};
