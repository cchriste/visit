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
#include <avtFrustumSelection.h>
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
: avtMTMDFileFormat(filename)
{
    resolution = 3;

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

    ApplicationPtr app = Application::New();
    app->setCommandLine(0, NULL);
    ENABLE_VISUS_IDX();

    string name(filename);
    name = "file://" + name;
    dataset = Dataset::open(URL(name));
    VisusReleaseAssert(dataset);
    
    selList = std::vector<avtDataSelection_p>();
    selsApplied = NULL;

    //Redundant Code
    NdBox world_box = dataset->logic_box;
    
    accessdata = dataset->createAccess();
    
    bounds[0] = world_box.to.x - world_box.from.x + 1;
    bounds[1] = world_box.to.y - world_box.from.y + 1;
    bounds[2] = world_box.to.z - world_box.from.z + 1;
    cerr << "dataset->logic_box=("
         << world_box.from.x << "," << world_box.to.x << "),("
         << world_box.from.y << "," << world_box.to.y << "),("
         << world_box.from.z << "," << world_box.to.z << ")" << endl;
    cerr << "dataset->maxh=" << dataset->maxh << endl;

    // Convert the logical extents into physical extents.
    extents[0] = world_box.from.x;
    extents[1] = world_box.to.x + 1;

    extents[2] = world_box.from.y;
    extents[3] = world_box.to.y + 1;

    extents[4] = world_box.from.z;
    extents[5] = world_box.to.z + 1;
}


// ****************************************************************************
//  Method: avtIDXFileFormat destructor
//
//  Programmer: Cameron Christensen
//  Creation:   Monday, November 04, 2013
//
// ****************************************************************************

avtIDXFileFormat::~avtIDXFileFormat()
{
    //delete selsApplied; //don't think we own this...
}


// ****************************************************************************
//  Method: avtIDXFileFormat::GetNTimesteps
//
//  Purpose:
//      Tells the rest of the code how many timesteps there are in this file.
//
//  Programmer: bremer5 -- generated by xml2avt
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

int
avtIDXFileFormat::GetNTimesteps(void)
{
    int NTimesteps = dataset->time->to - dataset->time->from + 1;
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
avtIDXFileFormat::FreeUpResources(void)
{
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
avtIDXFileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md,
    int timeState) 
{
    avtMeshType mt = AVT_RECTILINEAR_MESH;
    string mesh_for_this_var;
    mesh_for_this_var.assign("CC_Mesh");
    int block_origin = 0;
    int spatial_dimension = 3;
    int topological_dimension = 3;
    int ndtype;

    //
    // dynamic decomposition
    //
    md->SetFormatCanDoDomainDecomposition(true);

    //
    // Add the mesh. For formats that do their own decomposition, VisIt
    // expects we only advertise a single block.
    //
    int metadata_nblocks = 1;

    //This is MUCH faster, just provide the bounds that actually encompass the data we will be loading.
    //TODO: now how do we change them when we up the resolution??
    int res = resolution;
    int resReduction = 1 << res;
    int reducedBounds[3];

    reducedBounds[0] = bounds[0] % resReduction == 0 ?
        bounds[0] / resReduction : bounds[0] / resReduction + 1;
    reducedBounds[1] = bounds[1] % resReduction == 0 ?
        bounds[1] / resReduction : bounds[1] / resReduction + 1;
    reducedBounds[2] = bounds[2] % resReduction == 0 ?
        bounds[2] / resReduction : bounds[2] / resReduction + 1;

    AddMeshToMetaData(md, mesh_for_this_var, mt, extents, metadata_nblocks, block_origin, spatial_dimension, topological_dimension, reducedBounds);
    //AddMeshToMetaData(md, mesh_for_this_var, mt, extents, metadata_nblocks, block_origin, spatial_dimension, topological_dimension, bounds);

    //
    // We want to set the LOD property of the Mesh Meta Data. Since we only
    // have one mesh we can assume that md->GetMeshes(0) points to the
    // avtMeshMetaData Object created with the `AddMeshToMetaData' helper.
    //
    // Advertise 3 levels for now
    //
    md->GetMeshes(0).LODs = (dataset->maxh - 15) / 3;

    //
    // Add the scalars and vectors.
    //
    std::vector<string> fieldnames = dataset->getFieldNames();
    
    for (int i = 0; i < (int) fieldnames.size(); i++)
    {
        FieldPtr field = dataset->getFieldByName(fieldnames[i]);
        avtCentering cent = AVT_ZONECENT;
        ndtype=1;
        VectorDTypePtr vdtype=VectorDType::cast(field->dtype);
        if (vdtype)
            ndtype=vdtype->num;
        if (ndtype == 1)
            AddScalarVarToMetaData(md, fieldnames[i], mesh_for_this_var, cent);
        else
            AddVectorVarToMetaData(md, fieldnames[i], mesh_for_this_var, cent, ndtype);
    }
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
avtIDXFileFormat::GetMesh(int timestate, int domain, const char *meshname)
{
    NdBox slice_box = dataset->logic_box;

    cerr << "avtIDXFileFormat::GetMesh: timestate("<<timestate<<") domain("<<domain<<") meshname("<<meshname<<") resolution("<<resolution<<")\n";
    //
    // Determine the reduced grid size taking into account the multi-res
    // reduction.
    //
    int res = resolution;
    int resReduction = 1 << res;
    int reducedBounds[3];

    reducedBounds[0] = bounds[0] % resReduction == 0 ?
        bounds[0] / resReduction : bounds[0] / resReduction + 1;
    reducedBounds[1] = bounds[1] % resReduction == 0 ?
        bounds[1] / resReduction : bounds[1] / resReduction + 1;
    reducedBounds[2] = bounds[2] % resReduction == 0 ?
        bounds[2] / resReduction : bounds[2] / resReduction + 1;

    //
    // Partition into slices along the Z direction. Assumes cell centered
    // data. This may fail if there are not enough cells to populate the
    // last block.
    //
    int zWidth = reducedBounds[2] % nblocks == 0 ?
        reducedBounds[2] / nblocks : reducedBounds[2] / nblocks + 1;
    if(reducedBounds[2] % nblocks == 0)
    {
        slice_box.from[2] = zWidth * resReduction * domain;
        slice_box.to[2]   = slice_box.from[2] + zWidth * resReduction - 1;
    }
    else
    {
        if(domain == nblocks - 1)
        {
            slice_box.from[2] = zWidth * resReduction * domain;
            slice_box.to[2]   = slice_box.to[2];
        }
        else
        {       
            slice_box.from[2] = zWidth * resReduction * domain;
            slice_box.to[2]   = slice_box.from[2] + zWidth * resReduction - 1;
        }
    }
    
    //
    // Create the mesh.
    //
    vtkRectilinearGrid *rgrid = vtkRectilinearGrid::New();
    int dims[3];
    float *arrayX;
    float *arrayZ;
    float *arrayY;
    vtkFloatArray *coordsX;
    vtkFloatArray *coordsY;
    vtkFloatArray *coordsZ;
    
    dims[0] = (slice_box.to[0] - slice_box.from[0] + 1) / resReduction + 1;
    dims[1] = (slice_box.to[1] - slice_box.from[1] + 1) / resReduction + 1;
    dims[2] = (slice_box.to[2] - slice_box.from[2] + 1) / resReduction + 1;
    rgrid->SetDimensions(dims[0], dims[1], dims[2]);
    cerr<<"avtIDXFileFormat::GetMesh() returning <"<<dims[0]<<"x"<<dims[1]<<"x"<<dims[2]<<"> mesh for resolution "<<resolution<<" (resReduction="<<resReduction<<")\n";

    coordsX = vtkFloatArray::New();
    coordsX->SetNumberOfTuples(dims[0]);
    arrayX = (float *) coordsX->GetVoidPointer(0);
    for (int i = 0; i < dims[0]; i++)
        arrayX[i] = i * resReduction;
    arrayX[dims[0]-1] = slice_box.to[0] + 1.;
    rgrid->SetXCoordinates(coordsX);

    coordsY = vtkFloatArray::New();
    coordsY->SetNumberOfTuples(dims[1]);
    arrayY = (float *) coordsY->GetVoidPointer(0);
    for (int i = 0; i < dims[1]; i++)
        arrayY[i] = i * resReduction;
    arrayY[dims[1]-1] = slice_box.to[1] + 1.;
    rgrid->SetYCoordinates(coordsY);

    coordsZ = vtkFloatArray::New();
    coordsZ->SetNumberOfTuples(dims[2]);
    arrayZ = (float *) coordsZ->GetVoidPointer(0);
    for (int i = 0; i < dims[2]; i++)
        arrayZ[i] = slice_box.from[2] + i * resReduction;
    arrayZ[dims[2]-1] = slice_box.to[2] + 1.;
    rgrid->SetZCoordinates(coordsZ);

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
avtIDXFileFormat::GetVar(int timestate, int domain, const char *varname)
{
    cerr << "avtIDXFileFormat::GetVar: timestate("<<timestate<<") domain("<<domain<<") varname("<<varname<<") resolution("<<resolution<<")\n";
    string name(varname);
    NdBox slice_box = dataset->logic_box;

    //
    // Determine the reduced grid size taking into account the multi-res
    // reduction.
    //
    int res = resolution;
    int resReduction = 1 << res;
    int reducedBounds[3];

    reducedBounds[0] = bounds[0] % resReduction == 0 ?
        bounds[0] / resReduction : bounds[0] / resReduction + 1;
    reducedBounds[1] = bounds[1] % resReduction == 0 ?
        bounds[1] / resReduction : bounds[1] / resReduction + 1;
    reducedBounds[2] = bounds[2] % resReduction == 0 ?
        bounds[2] / resReduction : bounds[2] / resReduction + 1;

    //
    // Partition into slices along the Z direction. Assumes cell centered
    // data. This may fail if there are not enough cells to populate the
    // last block.
    //
    int zWidth = reducedBounds[2] % nblocks == 0 ?
        reducedBounds[2] / nblocks : reducedBounds[2] / nblocks + 1;
    if(reducedBounds[2] % nblocks == 0)
    {
        slice_box.from[2] = zWidth * resReduction * domain;
        slice_box.to[2]   = slice_box.from[2] + zWidth * resReduction - 1;
    }
    else
    {
        if(domain == nblocks - 1)
        {
            slice_box.from[2] = zWidth * resReduction * domain;
            slice_box.to[2]   = slice_box.to[2];
        }
        else
        {       
            slice_box.from[2] = zWidth * resReduction * domain;
            slice_box.to[2]   = slice_box.from[2] + zWidth * resReduction - 1;
        }
    }

    
    FieldPtr field; 
    field = dataset->getFieldByName(varname);
    
    int maxh_read = dataset->maxh - res * 3;
    BoxQueryPtr query = dataset->createBoxQuery(Position::LogicSpace,
                        slice_box, dataset->getFieldByName(name), timestate, 0,
                        maxh_read, dataset->maxh);
    int ntuples = query->dims.innerProduct();
    VisusReleaseAssert(query);


    //mark frustum/resolution selections applied
    for (int i=0;i<selList.size();i++)
    {
        if (strcmp(selList[i]->GetType(),"avtResolutionSelection") == 0)
        {
            cerr<<"applying new (forced) resolution...(todo)\n";
            (*selsApplied)[i]=true;
        }
        else if (strcmp(selList[i]->GetType(), "avtFrustumSelection") == 0)
        {
            cerr<<"applying new frustum...(todo)\n";
            (*selsApplied)[i]=true;
        }
    }    

    if (field->dtype->equals(DType::UINT8))
    {
        vtkUnsignedCharArray*rv = vtkUnsignedCharArray::New();
        rv->SetNumberOfComponents(1);
        rv->SetNumberOfTuples(ntuples);
        
        VisusReleaseAssert(rv->GetPointer(0) != NULL);

        query->data = Array::New(query->dims, DType::UINT8, (unsigned char*) rv->GetPointer(0));
        
        VisusReleaseAssert(query->data->c_size() == sizeof (unsigned char) * ntuples);
 
        VisusReleaseAssert(query->execute(accessdata, 'r') == QuerySucceed);
        return rv;
    }
    else if (field->dtype->equals(DType::FLOAT32))
    {
        vtkFloatArray *rv = vtkFloatArray::New();
        rv->SetNumberOfComponents(1);
        rv->SetNumberOfTuples(ntuples);

        VisusReleaseAssert(rv->GetPointer(0) != NULL);
        query->data = Array::New(query->dims, DType::FLOAT32, (unsigned char*) rv->GetPointer(0));
        
        cerr << "IDX setup the array's" << endl;

        VisusReleaseAssert(query->data->c_size() == sizeof (float) * ntuples);
 
        VisusReleaseAssert(query->execute(accessdata, 'r') == QuerySucceed);

        cerr << "IDX gotten the data" << endl;

        return rv;
    }
    else if (field->dtype->equals(DType::FLOAT64))
    {
        vtkDoubleArray *rv1 = vtkDoubleArray::New();
        rv1->SetNumberOfComponents(1);
        rv1->SetNumberOfTuples(ntuples);
        
        VisusReleaseAssert(rv1->GetPointer(0) != NULL);
        query->data = Array::New(query->dims, DType::FLOAT64, (unsigned char*) rv1->GetPointer(0));
        
        VisusReleaseAssert(query->data->c_size() == sizeof (double) * ntuples);
 
        VisusReleaseAssert(query->execute(accessdata, 'r') == QuerySucceed);
        
        return rv1;
    }
    else if (field->dtype->equals(DType::UINT16))
    {
        vtkShortArray *rv1 = vtkShortArray::New();
        rv1->SetNumberOfComponents(1);
        rv1->SetNumberOfTuples(ntuples);

        VisusReleaseAssert(rv1->GetPointer(0) != NULL);
        query->data = Array::New(query->dims, DType::UINT16, (unsigned char*) rv1->GetPointer(0));
        
        VisusReleaseAssert(query->data->c_size() == sizeof (short) * ntuples);
 
        VisusReleaseAssert(query->execute(accessdata, 'r') == QuerySucceed);
        
        return rv1;
    }

    return NULL;
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
avtIDXFileFormat::GetVectorVar(int timestate, int domain, const char *varname)
{
    return NULL;
}


// ****************************************************************************
//  Method: avtIDXFileFormat::RegisterDataSelections
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

void
avtIDXFileFormat::RegisterDataSelections(
    const std::vector<avtDataSelection_p>& sels, std::vector<bool>* applied)
{
    this->selList = sels;
    this->selsApplied = applied;

    cerr << "avtIDXFileFormat::RegisterDataSelections" << endl;
    for(size_t i=0; i < sels.size(); ++i)
    {
        if(strcmp(sels[i]->GetType(), "avtResolutionSelection") == 0)
        {
            const avtResolutionSelection* sel = static_cast<const avtResolutionSelection*>(*sels[i]);
            cerr<<"new resolution: "<<sel->resolution()<<", (old resolution: "<<this->resolution<<")\n";
            if (this->resolution!=sel->resolution())
                this->resolution = sel->resolution();
            else
                (*applied)[i] = true;
        }
        else if (strcmp(sels[i]->GetType(), "avtFrustumSelection") == 0)
        {
            const avtFrustumSelection* sel = static_cast<const avtFrustumSelection*>(*sels[i]);
            cerr<<"new frustum\n";
            //todo: calculate new query bounds based on view
            for (int i=0;i<6;i++)
            {
                frustum2d[i]=sel->frustum2d[i];
                frustum3d[i]=sel->frustum3d[i];
            }
        }
        else if (strcmp(sels[i]->GetType(), "Resample Data Selection") == 0)
        {
            ;//ignore
        }
        else
            cerr<<"Error: unhandled selection "<<sels[i]->GetType()<<"!\n";
    }
}
