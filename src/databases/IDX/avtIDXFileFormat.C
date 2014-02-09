/*****************************************************************************
 *
 * Copyright (c) 2000 - 2012, Lawrence Livermore National Security, LLC
 * Produced at the Lawrence Livermore National Laboratory
 * LLNL-CODE-442911
 * All rights reserved.
 *
 * This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
 * full copyright notice is contained in the file COPYRIGHT located at the root
 * of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
 *
 * Redistribution  and  use  in  source  and  binary  forms,  with  or  without
 * modification, are permitted provided that the following conditions are met:
 *
 *  - Redistributions of  source code must  retain the above  copyright notice,
 *    this list of conditions and the disclaimer below.
 *  - Redistributions in binary form must reproduce the above copyright notice,
 *    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
 *    documentation and/or other materials provided with the distribution.
 *  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
 * ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
 * LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
 * DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
 * CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
 * LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
 * OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 *
 *****************************************************************************/

// ************************************************************************* //
//                            avtIDXFileFormat.C                           //
// ************************************************************************* //

#include <avtIDXFileFormat.h>

#include <string>
#include <vtkShortArray.h>
#include <vtkFloatArray.h>
#include <vtkRectilinearGrid.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <avtResolutionSelection.h>
#include <avtDatabaseMetaData.h>

#include <DBOptionsAttributes.h>
#include <Expression.h>

#include <InvalidVariableException.h>

#include <visuscpp/kernel/visus_kernel.h>
#include <visuscpp/db/visus_dataset.h>

#ifdef PARALLEL
#include <mpi.h>
#include <avtParallel.h>
#endif

#include <visuscpp/idx/visus_idx.h>

#include <visuscpp/kernel/io/visus_file.h>
#include <visuscpp/kernel/io/visus_fileinputstream.h>
#include <visuscpp/kernel/io/visus_fileoutputstream.h>
#include <visuscpp/kernel/array/visus_array.h>
#include <visuscpp/db/visus_multiplexaccess.h>
#include <visuscpp/db/visus_filteraccess.h>
#include <visuscpp/db/visus_ramaccess.h>
#include <visuscpp/db/visus_networkaccess.h>
#include <visuscpp/idx/visus_idx_dataset.h>
#include <visuscpp/idx/visus_idx_diskaccess.h>
#include <visuscpp/idx/visus_idx_hzorder.h>
#include <visuscpp/idx/visus_idx_pointquery.h>
#include <visuscpp/idx/visus_idx_boxquery.h>

#include <visuscpp/kernel/geometry/visus_vector3d.h>
#include <visuscpp/kernel/net/visus_url.h>
#include <visuscpp/kernel/net/visus_net_request.h>
#include <visuscpp/db/visus_field.h>
#include <visuscpp/db/visus_dataset.h>
#include <visuscpp/db/visus_access.h>
#include <visuscpp/db/visus_ramaccess.h>
#include <visuscpp/db/visus_boxquery.h>
#include <visuscpp/db/visus_pointquery.h>
#include <visuscpp/db/visus_dataset.h>

#include <visuscpp/kernel/visus_kernel.h>
#include <visuscpp/kernel/geometry/visus_quaternion.h>
#include <visuscpp/db/visus_db.h>
#include <visuscpp/idx/visus_idx_modvisus.h>
#include <visuscpp/gui/visus_guiapplication.h>
#include <visuscpp/appkit/visus_appkit.h>
#include <visuscpp/appkit/viewer/visus_viewer.h>


#include <visuscpp/kernel/array/visus_array.h>
#include <visuscpp/idx/visus_idx.h>

#include <visuscpp/kernel/array/visus_array.h>
#include <visuscpp/kernel/visus_kernel.h>
#include <visuscpp/kernel/array/visus_dtype.h>
#include <visuscpp/kernel/array/visus_range.h>
#include <visuscpp/kernel/geometry/visus_position.h>
#include <visuscpp/kernel/datastructure/visus_content.h>


using std::string;
using namespace VISUS_NAMESPACE;

//Global Variable required by ViSUS
static DatasetPtr dataset;
static int bounds[3];
static double extents[6];
static AccessPtr accessdata;
// ****************************************************************************
//  Method: avtIDXFileFormat constructor
//
//  Programmer: bremer5 -- generated by xml2avt
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

avtIDXFileFormat::avtIDXFileFormat(const char *filename)
: avtMTMDFileFormat(filename) {


    ApplicationPtr app = Application::New();
    app->setCommandLine(0, NULL);
    ENABLE_VISUS_IDX();

    string name(filename);
    name = "file://" + name;
    //std::cout << " [1] Opening (A): " << name << "\n";
    dataset = Dataset::open(URL(name));
    VisusReleaseAssert(dataset);
    
    //Redundant Code
    NdBox world_box = dataset->logic_box;
    //std::cout << "[2] Box dimension: " << world_box.to.x << ", " << world_box.to.y  << ", " << world_box.to.z  << "\n";
    
    accessdata = dataset->createAccess();
    //std::cout << "[3] Access GRANTED to IDX file\n";
    
    bounds[0] = world_box.to.x - world_box.from.x + 1;
    bounds[1] = world_box.to.y - world_box.from.y + 1;
    bounds[2] = world_box.to.z - world_box.from.z + 1;
    
    //TODO : USE logic to physic from .idx
    extents[0] = world_box.from.x ;
    extents[1] = world_box.to.x;
    extents[2] = world_box.from.y;
    
    extents[3] = world_box.to.y;
    extents[4] = world_box.from.z;
    extents[5] = world_box.to.z;

#ifdef PARALLEL
    rank = PAR_Rank();
#else
    rank = 0;
#endif

#ifdef PARALLEL
    nblocks = PAR_Size();
#else
    nblocks = 1;
#endif


}


// ****************************************************************************
//  Method: avtEMSTDFileFormat::GetNTimesteps
//
//  Purpose:
//      Tells the rest of the code how many timesteps there are in this file.
//
//  Programmer: bremer5 -- generated by xml2avt
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

int
avtIDXFileFormat::GetNTimesteps(void) {
    int NTimesteps = 1 + dataset->time->to - dataset->time->from;
    //std::cout << "[TS] Requesting Time-step info : " << NTimesteps << " [" << dataset->time->to << ", " << dataset->time->from << "] \n";
    return NTimesteps;
}


// ****************************************************************************
//  Method: avtIDXFileFormat::FreeUpResources
//
//  Purpose:
//      When VisIt is done focusing on a particular timestep, it asks that
//      timestep to free up any resources (memory, file descriptors) that
//      it has associated with it.  This method is the mechanism for doing
//      that.
//
//  Programmer: bremer5 -- generated by xml2avt
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

void
avtIDXFileFormat::FreeUpResources(void) {
    //std::cout << "Freeing up resources \n";
}


// ****************************************************************************
//  Method: avtIDXFileFormat::PopulateDatabaseMetaData
//
//  Purpose:
//      This database meta-data object is like a table of contents for the
//      file.  By populating it, you are telling the rest of VisIt what
//      information it can request from you.
//
//  Programmer: bremer5 -- generated by xml2avt
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

void
avtIDXFileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md, int timeState) 
{
    ///
    /// CYRUS NOTE:
    ///

    // resolution is set according to the MultiresControl Operator slider
    cout << rank << " selected resolution = " << resolution << endl;

    //std::cout << "[PDM 1] Populate Database Meta-data for timeState: " << timeState << "\n";

    //
    // CODE TO ADD A MESH
    //
    //string meshname = ...
    //
    // AVT_RECTILINEAR_MESH, AVT_CURVILINEAR_MESH, AVT_UNSTRUCTURED_MESH,
    // AVT_POINT_MESH, AVT_SURFACE_MESH, AVT_UNKNOWN_MESH

    int i;
    int decider = 0;
    avtMeshType mt = AVT_RECTILINEAR_MESH;
    string mesh_for_this_var;
    mesh_for_this_var.assign("CC_Mesh");
    int block_origin = 0;
    int spatial_dimension = 3;
    int topological_dimension = 3;
    int ndtype;
    //double bounds[3];
    
    

    //std::cout << "[PDM 2] Extent: [" << bounds[0] << ", " << bounds[1] << ", " << bounds[2] << "]\n";

    // dynamic decomposition
    md->SetFormatCanDoDomainDecomposition(true);

    // Here's the call that tells the meta-data object that we have a mesh:

    ///
    /// CYRUS NOTE:
    ///
    // For formats that do their own decomposition, VisIt expects we only advertise a single block.
    int metadata_nblocks = 1;

    AddMeshToMetaData(md, mesh_for_this_var, mt, extents, metadata_nblocks, block_origin, spatial_dimension, topological_dimension, bounds);

    ///
    /// CYRUS NOTE:
    ///
    // we want to set the LODs property of  the Mesh Meta Data. Since we only have one mesh
    // we can assume that md->GetMeshes(0) points to the avtMeshMetaData Object created
    // with the `AddMeshToMetaData' helper

    // Advertise 3 levels for now
    md->GetMeshes(0).LODs = 3;

    //std::cout << "[PDM 3] Mesh Added\n";

    std::vector<string> fieldnames = dataset->getFieldNames();
    FieldPtr field;
    
    for (i = 0; i < (int) fieldnames.size(); i++) {
        field = dataset->getFieldByName(fieldnames[i]);
        avtCentering cent = AVT_ZONECENT; // AVT_NODECENT;
        ndtype=1;
        VectorDTypePtr vdtype=VectorDType::cast(field->dtype);
        if (vdtype)
                ndtype=vdtype->num;
        printf("Number of samples per variable %d\n", ndtype);
        if (ndtype == 1){
            AddScalarVarToMetaData(md, fieldnames[i], mesh_for_this_var, cent);
            printf("i am currently here\n");
        }
        else
            AddVectorVarToMetaData(md, fieldnames[i], mesh_for_this_var, cent, ndtype);
    }

    //
    // CODE TO ADD A SCALAR VARIABLE
    //
    // string mesh_for_this_var = meshname; // ??? -- could be multiple meshes
    // string varname = ...
    //
    // AVT_NODECENT, AVT_ZONECENT, AVT_UNKNOWN_CENT
    // avtCentering cent = AVT_NODECENT;
    //
    //
    // Here's the call that tells the meta-data object that we have a var:
    //
    // AddScalarVarToMetaData(md, varname, mesh_for_this_var, cent);
    //

    //
    // CODE TO ADD A VECTOR VARIABLE
    //
    // string mesh_for_this_var = meshname; // ??? -- could be multiple meshes
    // string varname = ...
    // int vector_dim = 2;
    //
    // AVT_NODECENT, AVT_ZONECENT, AVT_UNKNOWN_CENT
    // avtCentering cent = AVT_NODECENT;
    //
    //
    // Here's the call that tells the meta-data object that we have a var:
    //
    // AddVectorVarToMetaData(md, varname, mesh_for_this_var, cent,vector_dim);
    //

    //
    // CODE TO ADD A TENSOR VARIABLE
    //
    // string mesh_for_this_var = meshname; // ??? -- could be multiple meshes
    // string varname = ...
    // int tensor_dim = 9;
    //
    // AVT_NODECENT, AVT_ZONECENT, AVT_UNKNOWN_CENT
    // avtCentering cent = AVT_NODECENT;
    //
    //
    // Here's the call that tells the meta-data object that we have a var:
    //
    // AddTensorVarToMetaData(md, varname, mesh_for_this_var, cent,tensor_dim);
    //

    //
    // CODE TO ADD A MATERIAL
    //
    // string mesh_for_mat = meshname; // ??? -- could be multiple meshes
    // string matname = ...
    // int nmats = ...;
    // vector<string> mnames;
    // for (int i = 0 ; i < nmats ; i++)
    // {
    //     char str[32];
    //     sprintf(str, "mat%d", i);
    //     -- or -- 
    //     strcpy(str, "Aluminum");
    //     mnames.push_back(str);
    // }
    // 
    // Here's the call that tells the meta-data object that we have a mat:
    //
    // AddMaterialToMetaData(md, matname, mesh_for_mat, nmats, mnames);
    //
    //
    // Here's the way to add expressions:
    //Expression momentum_expr;
    //momentum_expr.SetName("momentum");
    //momentum_expr.SetDefinition("{u, v}");
    //momentum_expr.SetType(Expression::VectorMeshVar);
    //md->AddExpression(&momentum_expr);
    //Expression KineticEnergy_expr;
    //KineticEnergy_expr.SetName("KineticEnergy");
    //KineticEnergy_expr.SetDefinition("0.5*(momentum*momentum)/(rho*rho)");
    //KineticEnergy_expr.SetType(Expression::ScalarMeshVar);
    //md->AddExpression(&KineticEnergy_expr);
    //
}


// ****************************************************************************
//  Method: avtIDXFileFormat::GetMesh
//
//  Purpose:
//      Gets the mesh associated with this file.  The mesh is returned as a
//      derived type of vtkDataSet (ie vtkRectilinearGrid, vtkStructuredGrid,
//      vtkUnstructuredGrid, etc).
//
//  Arguments:
//      timestate   The index of the timestate.  If GetNTimesteps returned
//                  'N' time steps, this is guaranteed to be between 0 and N-1.
//      domain      The index of the domain.  If there are NDomains, this
//                  value is guaranteed to be between 0 and NDomains-1,
//                  regardless of block origin.
//      meshname    The name of the mesh of interest.  This can be ignored if
//                  there is only one mesh.
//
//  Programmer: bremer5 -- generated by xml2avt
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

vtkDataSet *
avtIDXFileFormat::GetMesh(int timestate, int domain, const char *meshname) {
    //std::cout << "[GetMesh] Time state: " << timestate << " Domain: " << domain << " Meshname: " << meshname << "\n";


cout << "GetMesh rank = " << rank << " " << " nblocks = " << nblocks << " domain # = "  << domain <<endl;

    NdBox slice_box = dataset->logic_box;
    int zWidth =  slice_box.to[2] / nblocks; // bounds[2]/nblocks;
    printf("M[%d %s] : %d\n", domain, meshname, slice_box.to[2]);
    if(slice_box.to[2] % 2 == 0)
    {
        //case when the z dimension is an odd number
        if(domain == nblocks -1)
        {
                slice_box.from[2]= zWidth * domain;//(bounds[2] * domain) / nblocks;// nslice;
                slice_box.to[2]=slice_box.from[2] + zWidth + 0;
        }
        else
        {       
                slice_box.from[2]= zWidth * domain;//(bounds[2] * domain) / nblocks;// nslice;
                slice_box.to[2]=slice_box.from[2] + zWidth - 1;
        }
    }
    else
    {
                slice_box.from[2]= zWidth * domain;//(bounds[2] * domain) / nblocks;// nslice;
                slice_box.to[2]=slice_box.from[2] + zWidth - 1;
    }
    

    
    vtkRectilinearGrid *rgrid = vtkRectilinearGrid::New();
    int dims[3];
    float *arrayX;
    float *arrayZ;
    float *arrayY;
    vtkFloatArray *coordsX;
    vtkFloatArray *coordsY;
    vtkFloatArray *coordsZ;
    
    dims[0] = (slice_box.to[0] - slice_box.from[0] + 2) ;//*/bounds[0];
    dims[1] = (slice_box.to[1] - slice_box.from[1] + 2) ;//*/bounds[1];
    dims[2] = (slice_box.to[2] - slice_box.from[2] + 2) ;//*/bounds[2];
    rgrid->SetDimensions(dims[0], dims[1], dims[2]);

    printf("[%d] DIMS : [%d %d %d] :: [%d]\n", domain, dims[0], dims[1], dims[2], dims[0]*dims[1]*dims[2]);
    //printf("[%d] : From : TO :: %d : %d\n", domain, slice_box.from[2], slice_box.to[2]);
   
        coordsX = vtkFloatArray::New();
        coordsX->SetNumberOfTuples(dims[0]);
        arrayX = (float *) coordsX->GetVoidPointer(0);
        for (int i = 0; i < dims[0]; i++) {
            arrayX[i] = i;
        }
        rgrid->SetXCoordinates(coordsX);
        
        coordsY = vtkFloatArray::New();
        coordsY->SetNumberOfTuples(dims[1]);
        arrayY = (float *) coordsY->GetVoidPointer(0);
        for (int i = 0; i < dims[1]; i++) {
            arrayY[i] = i;
        }
        rgrid->SetYCoordinates(coordsY);
        
        coordsZ = vtkFloatArray::New();
        coordsZ->SetNumberOfTuples(dims[2]);
        arrayZ = (float *) coordsZ->GetVoidPointer(0);
        for (int i = 0; i < dims[2]; i++) {
            arrayZ[i] = slice_box.from[2] + i;
        }
        rgrid->SetZCoordinates(coordsZ);

        
    
    
    /*
    for (int c = 0; c < 3; c++) {
        vtkFloatArray *coords = vtkFloatArray::New();
        coords->SetNumberOfTuples(dims[c]);
        float *array = (float *) coords->GetVoidPointer(0);
        for (int i = 0; i < dims[c]; i++) {
            array[i] = i;
        }

        switch (c) {
            case 0:
                rgrid->SetXCoordinates(coords);
                break;
            case 1:
                rgrid->SetYCoordinates(coords);
                break;
            case 2:
                rgrid->SetZCoordinates(coords);
                break;
        }

        coords->Delete();
    }
    */

    return rgrid;

}


// ****************************************************************************
//  Method: avtIDXFileFormat::GetVar
//
//  Purpose:
//      Gets a scalar variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      timestate  The index of the timestate.  If GetNTimesteps returned
//                 'N' time steps, this is guaranteed to be between 0 and N-1.
//      domain     The index of the domain.  If there are NDomains, this
//                 value is guaranteed to be between 0 and NDomains-1,
//                 regardless of block origin.
//      varname    The name of the variable requested.
//
//  Programmer: bremer5 -- generated by xml2avt
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

vtkDataArray *
avtIDXFileFormat::GetVar(int timestate, int domain, const char *varname) {
   // std::cout << "[GetVarNN] Time state: " << timestate << "Domain: " << domain << "variable name: " << varname << "\n";


    printf("V[%d %s]\n", domain, varname);


    string name(varname);
    NdBox slice_box = dataset->logic_box;
    int zWidth =  slice_box.to[2] / nblocks; 
    if(slice_box.to[2] % 2 == 0)
    {
        //case when the z dimension is an odd number
        if(domain == nblocks -1)
        {
                slice_box.from[2]= zWidth * domain;//(bounds[2] * domain) / nblocks;// nslice;
                slice_box.to[2]=slice_box.from[2] + zWidth + 0;
        }
        else
        {       
                slice_box.from[2]= zWidth * domain;//(bounds[2] * domain) / nblocks;// nslice;
                slice_box.to[2]=slice_box.from[2] + zWidth - 1;
        }
    }
    else
    {
                slice_box.from[2]= zWidth * domain;//(bounds[2] * domain) / nblocks;// nslice;
                slice_box.to[2]=slice_box.from[2] + zWidth - 1;
    }
    int ntuples =  (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1);//bounds[0] * bounds[1] * bounds[2] / nblocks; // this is the number of entries in the variable.
    
    FieldPtr field; 
    field = dataset->getFieldByName(varname);
    
   // vtkDoubleArray *rv1 = vtkFloatArray::New();
    
    
    
    if (field->dtype->equals(DType::UINT8))
    {
        vtkUnsignedCharArray*rv = vtkUnsignedCharArray::New();
        rv->SetNumberOfComponents(1);
        rv->SetNumberOfTuples(ntuples);
        BoxQueryPtr query = dataset->createBoxQuery(Position::LogicSpace, slice_box, dataset->getFieldByName(name), timestate, 0, dataset->maxh, dataset->maxh);
        VisusReleaseAssert(query && query->dims.innerProduct() == (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1));
        NdPoint dims((slice_box.to[0] - slice_box.from[0] + 1), (slice_box.to[1] - slice_box.from[1] + 1), (slice_box.to[2] - slice_box.from[2] + 1), 1, 1);
        
        VisusReleaseAssert(rv->GetPointer(0) != NULL);
        query->data = Array::New(dims, DType::UINT8, (unsigned char*) rv->GetPointer(0)); //);
        
        VisusReleaseAssert(query->data->c_size() == sizeof (unsigned char) * (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1));
 
        VisusReleaseAssert(query->execute(accessdata, 'r') == QuerySucceed);
        return rv;
    }
    
    else if (field->dtype->equals(DType::FLOAT32))
    {
        vtkFloatArray *rv = vtkFloatArray::New();
        rv->SetNumberOfComponents(1);
        rv->SetNumberOfTuples(ntuples);
        BoxQueryPtr query = dataset->createBoxQuery(Position::LogicSpace, slice_box, dataset->getFieldByName(name), timestate, 0, dataset->maxh, dataset->maxh);
        VisusReleaseAssert(query && query->dims.innerProduct() == (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1)/*bounds[0] * bounds[1] * bounds[2] / nblocks*/);
        NdPoint dims((slice_box.to[0] - slice_box.from[0] + 1), (slice_box.to[1] - slice_box.from[1] + 1), (slice_box.to[2] - slice_box.from[2] + 1), 1, 1);
        
        VisusReleaseAssert(rv->GetPointer(0) != NULL);
        query->data = Array::New(dims, DType::FLOAT32, (unsigned char*) rv->GetPointer(0)); //);
        
        VisusReleaseAssert(query->data->c_size() == sizeof (float) * (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1));
 
        VisusReleaseAssert(query->execute(accessdata, 'r') == QuerySucceed);
        return rv;
    }

    
    else if (field->dtype->equals(DType::FLOAT64))
    {
        //printf("[DOUBLE] 1\n");
        vtkDoubleArray *rv1 = vtkDoubleArray::New();
        //printf("[DOUBLE] 2\n");
        rv1->SetNumberOfComponents(1);
        //printf("[DOUBLE] 3\n");
        rv1->SetNumberOfTuples(ntuples);
        //printf("[DOUBLE] 4  : %d\n", ntuples);
        BoxQueryPtr query = dataset->createBoxQuery(Position::LogicSpace, slice_box, dataset->getFieldByName(name), timestate, 0, dataset->maxh, dataset->maxh);
        //printf("[DOUBLE] 5 : : %d : %d\n", query->dims.innerProduct()), (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1);
        VisusReleaseAssert(query && query->dims.innerProduct() == (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1));
        //printf("[DOUBLE] 6\n");
        NdPoint dims((slice_box.to[0] - slice_box.from[0] + 1), (slice_box.to[1] - slice_box.from[1] + 1), (slice_box.to[2] - slice_box.from[2] + 1), 1, 1);
        
        VisusReleaseAssert(rv1->GetPointer(0) != NULL);
        query->data = Array::New(dims, DType::FLOAT64, (unsigned char*) rv1->GetPointer(0)); //);
        
        VisusReleaseAssert(query->data->c_size() == sizeof (double) * (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1));
 
        VisusReleaseAssert(query->execute(accessdata, 'r') == QuerySucceed);
        
         return rv1;
    }
    
    else if (field->dtype->equals(DType::UINT16))
    {
        //printf("[DOUBLE] 1\n");
        vtkShortArray *rv1 = vtkShortArray::New();
        //printf("[DOUBLE] 2\n");
        rv1->SetNumberOfComponents(1);
        //printf("[DOUBLE] 3\n");
        rv1->SetNumberOfTuples(ntuples);
        //printf("[DOUBLE] 4  : %d\n", ntuples);
        BoxQueryPtr query = dataset->createBoxQuery(Position::LogicSpace, slice_box, dataset->getFieldByName(name), timestate, 0, dataset->maxh, dataset->maxh);
        //printf("[DOUBLE] 5 : : %d : %d\n", query->dims.innerProduct()), (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1);
        VisusReleaseAssert(query && query->dims.innerProduct() == (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1));
        //printf("[DOUBLE] 6\n");
        NdPoint dims((slice_box.to[0] - slice_box.from[0] + 1), (slice_box.to[1] - slice_box.from[1] + 1), (slice_box.to[2] - slice_box.from[2] + 1), 1, 1);
        
        VisusReleaseAssert(rv1->GetPointer(0) != NULL);
        query->data = Array::New(dims, DType::UINT16, (unsigned char*) rv1->GetPointer(0)); //);
        
        VisusReleaseAssert(query->data->c_size() == sizeof (short) * (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1));
 
        VisusReleaseAssert(query->execute(accessdata, 'r') == QuerySucceed);
        
         return rv1;
    }
    
    
    
    
    
   
    
    
    

    
    
  //  std::cout << "GetVar [2] " << bounds[0] << ", " << bounds[1] << ", " << bounds[2] << "Queried Data : "<< query->dims.innerProduct() << "Desired Queried data : " << (slice_box.to[0] - slice_box.from[0] + 1) * (slice_box.to[1] - slice_box.from[1] + 1) * (slice_box.to[2] - slice_box.from[2] + 1) <<  "\n";
    
    
    
    //(NdPoint dims,DTypePtr dtype,const unsigned char* buffer) 
    
    //dims.x = 
    
  //  std::cout << "HERE MM " << rv->GetPointer(0) << "\n";
    
   // std::cout << "HERE 3 " << query->data->c_size() << "\n";
    //read data from disk
 
    

   // std::cout << "HERE 4\n";

    //float* Src=(float*)query->data->c_ptr();

    //float* data_copy = (float*)malloc(sizeof(float) * bounds[0]*bounds[1]*bounds[2]);

    //for(int i = 0 ; i < bounds[0]*bounds[1]*bounds[2] ; i++)
    //std::cout<<"value at "<< i  << " = " << Src[i*10] <<"\n";
    //data_copy[i] = Src[i];


    //rv->SetNumberOfComponents(1);
    //rv->SetArray((float*)GetPointer(0), bounds[0]*bounds[1]*bounds[2], 0);



    //for (int i = 0 ; i < ntuples ; i++)
    //{
    //   rv->SetTuple1(i, data_copy[i]/*Src[i]*/);  // you must determine value for ith entry.
    //}

    //YOU MUST IMPLEMENT THIS

    //
    // If you have a file format where variables don't apply (for example a
    // strictly polygonal format like the STL (Stereo Lithography) format,
    // then uncomment the code below.
    //
    // EXCEPTION1(InvalidVariableException, varname);
    //

    //
    // If you do have a scalar variable, here is some code that may be helpful.
    //
    // int ntuples = XXX; // this is the number of entries in the variable.
    // vtkFloatArray *rv = vtkFloatArray::New();
    // rv->SetNumberOfTuples(ntuples);
    // for (int i = 0 ; i < ntuples ; i++)
    // {
    //      rv->SetTuple1(i, VAL);  // you must determine value for ith entry.
    // }
    //
    // return rv;
    //

    
}


// ****************************************************************************
//  Method: avtIDXFileFormat::GetVectorVar
//
//  Purpose:
//      Gets a vector variable associated with this file.  Although VTK has
//      support for many different types, the best bet is vtkFloatArray, since
//      that is supported everywhere through VisIt.
//
//  Arguments:
//      timestate  The index of the timestate.  If GetNTimesteps returned
//                 'N' time steps, this is guaranteed to be between 0 and N-1.
//      domain     The index of the domain.  If there are NDomains, this
//                 value is guaranteed to be between 0 and NDomains-1,
//                 regardless of block origin.
//      varname    The name of the variable requested.
//
//  Programmer: bremer5 -- generated by xml2avt
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

vtkDataArray *
avtIDXFileFormat::GetVectorVar(int timestate, int domain, const char *varname) {
    std::cout << "[GetVectorVar : B] Times tate: " << timestate << "Domain: " << domain << "variable name: " << varname << "\n";
    //YOU MUST IMPLEMENT THIS
    //
    // If you have a file format where variables don't apply (for example a
    // strictly polygonal format like the STL (Stereo Lithography) format,
    // then uncomment the code below.
    //
    // EXCEPTION1(InvalidVariableException, varname);
    //

    //
    // If you do have a vector variable, here is some code that may be helpful.
    //
    // int ncomps = YYY;  // This is the rank of the vector - typically 2 or 3.
    // int ntuples = XXX; // this is the number of entries in the variable.
    // vtkFloatArray *rv = vtkFloatArray::New();
    // int ucomps = (ncomps == 2 ? 3 : ncomps);
    // rv->SetNumberOfComponents(ucomps);
    // rv->SetNumberOfTuples(ntuples);
    // float *one_entry = new float[ucomps];
    // for (int i = 0 ; i < ntuples ; i++)
    // {
    //      int j;
    //      for (j = 0 ; j < ncomps ; j++)
    //           one_entry[j] = ...
    //      for (j = ncomps ; j < ucomps ; j++)
    //           one_entry[j] = 0.;
    //      rv->SetTuple(i, one_entry); 
    // }
    //
    // delete [] one_entry;
    // return rv;
    //

    return NULL;
}


// ****************************************************************************
//  Method: avtChomboFileFormat::RegisterDataSelections
//
//  Purpose:
//      Tries to read requests for specific resolutions.
//
//  Programmer: Tom Fogal
//  Creation:   August 5, 2010
//
//  Modifications:
//
// ****************************************************************************
void avtIDXFileFormat::RegisterDataSelections(
       const std::vector<avtDataSelection_p>& sels,
       std::vector<bool>* applied)
{
    for(size_t i=0; i < sels.size(); ++i)
    {
        if(strcmp(sels[i]->GetType(), "avtResolutionSelection") == 0)
        {
            const avtResolutionSelection* sel =
                static_cast<const avtResolutionSelection*>(*sels[i]);
            this->resolution = sel->resolution();
            (*applied)[i] = true;
        }
    }
}

