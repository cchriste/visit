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

#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCharArray.h>
#include <vtkShortArray.h>
#include <vtkIntArray.h>
#include <vtkLongArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkRenderWindow.h>
#include <vtkMatrix4x4.h>

#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <avtResolutionSelection.h>
#include <avtFrustumSelection.h>
#include <avtDatabaseMetaData.h>
#include <avtCallback.h>
#include <avtView3D.h>
#include <avtWorldSpaceToImageSpaceTransform.h>

#include <DBOptionsAttributes.h>
#include <Expression.h>

#include <InvalidVariableException.h>

#include <visuscpp/db/visus_db.h>

#ifdef PARALLEL
#include <mpi.h>
#include <avtParallel.h>
#endif

using std::string;
VISUS_USE_NAMESPACE

/////////////////////////////////////////////
class DummyNode : public DataflowNode , public IContentHolder
{
 public:

    //constructor
    DummyNode(String name="") : DataflowNode(name)
    {
        addInputPort("data");
        addOutputPort("data"); //if you do not need the original data, simply do not connect it!
    }

    //destructor
    virtual ~DummyNode()
    {}

    //getContentPhysicPosition (from IContentHolder class)
    virtual Position getContentPhysicPosition()
    {
        //if (!this->data)
            return Position::invalid();
        // else
        //     return Position(data->logic_to_physic,Box(Point3d(0,0,0),Point3d(data->dims.x,data->dims.y,data->dims.z)));
    }

    //from dataflow interface
    virtual bool processInput()
    {
        VisusAssert(false); //this really shouldn't be called anymore!
        //VisusInfo()<<"DummyNode::processInput... I guess we have some data. Ignore it (should get copied from publish)";
        // VisusAssert(node);

        // SharedPtr<Array> data=DynamicPointerCast<Array>(readInput("data"));

        // if (!data)
        // {
        //     VisusWarning()<<"Error: unable to read data";
        //     return false;
        // }
        // else
        // {
        //     VisusInfo()<<"YAY!!! read data of resolution <"<<data->dims.x<<","<<data->dims.y<<","<<data->dims.z<<">";
        //     node->bounds[0]=data->dims.x;
        //     node->bounds[1]=data->dims.y;
        //     node->bounds[2]=data->dims.z;
        //     node->data=(data);
        // }

        return true;
    }

 private:

    VISUS_DECLARE_NON_COPYABLE(DummyNode);


};


///////////////////////////////////////////////////////////
bool avtIDXQueryNode::publish(const DictObject &evt)
{
    SharedPtr<Array> data  =DynamicPointerCast<Array>   (evt.getattr("data"));
    SharedPtr<Position>dims=DynamicPointerCast<Position>(evt.getattr("dims"));
    SharedPtr<Position> pos=DynamicPointerCast<Position>(evt.getattr("position"));

    if (pos)
    {
        VisusInfo()<<"avtIDXQueryNode::publish: got position: "<<pos->box.toString();
        VisusInfo()<<"                          got dims:     "<<dims->box.toString();
        node->bounds[0] = dims->box.p2.x - dims->box.p1.x;
        node->bounds[1] = dims->box.p2.y - dims->box.p1.y;
        node->bounds[2] = dims->box.p2.z - dims->box.p1.z;

        if (node->resolution==1) //<ctc> use this to control whether or not to read new position (first few positions are very bad due to camera not being setup properly)
        {
            node->extents[0] = pos->box.p1.x;
            node->extents[1] = pos->box.p2.x;
            node->extents[2] = pos->box.p1.y;
            node->extents[3] = pos->box.p2.y;
            node->extents[4] = pos->box.p1.z;
            node->extents[5] = pos->box.p2.z;
        }
    }
    else if (data)
    {
        VisusInfo()<<"avtIDXQueryNode::publish: got data of resolution <"<<data->dims.x<<","<<data->dims.y<<","<<data->dims.z<<">";
        VisusAssert(node->bounds[0]==data->dims.x);
        VisusAssert(node->bounds[1]==data->dims.y);
        VisusAssert(node->bounds[2]==data->dims.z);
        node->haveData=true;
        node->data=data;
    }
    else
        VisusWarning()<<"Error in avtIDXQueryNode::publish";

    //return QueryNode::publish(evt); //don't even need to send it down the line.
    return true;
}

// ****************************************************************************
//  Method: avtIDXFileFormat constructor
//
//  Programmer: bremer5 -- generated by xml2avt
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

int avtIDXFileFormat::num_instances=0;

avtIDXFileFormat::avtIDXFileFormat(const char *filename)
: avtMTMDFileFormat(filename)
{
    resolution = 3;
    haveData=false;

    //<ctc> not sure this does anything or if we need to keep them (probably not since manual control does not happen soon enough)
    selList = std::vector<avtDataSelection_p>();
    selsApplied = NULL; 

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

    int dim=3; //<ctc> how to tell it when to extract only a slice?

    if (num_instances++<1)
    {
        app.reset(new Application);
        app->setCommandLine(0,NULL);
        ENABLE_VISUS_DB();
    }

    //dataflow
    this->dataflow.reset(new Dataflow);
    this->dataflow->oninput.connect(bind(&avtIDXFileFormat::onDataflowInput,this));

    string name(filename);

    //name="http://atlantis.sci.utah.edu/mod_visus?dataset=MRI_0901"; //<ctc> works!

    //<ctc> todo: add http:// url to ".urlidx" file and read it, or just open selected file directly
    //try to open a dataset
    //name = "file://" + name;
    dataset.reset(Dataset::loadDataset(name));
    if (!dataset)
    {
        VisusError()<<"could not load "<<name;
        VisusAssert(false); //<ctc> this shouldn't be done in the constructor: no way to fail if there is a problem.
    }

    //connect dataset
    DatasetNode*    dataflow_dataset   = new DatasetNode   ();
    dataflow_dataset->setDataset(dataset);

    VisusInfo()<<"creating the query node...";

    query=new avtIDXQueryNode(this);

    //position
    {
        Position* position=new Position();
        position->box=dataflow_dataset->getContentPhysicPosition().toAABB();
        if (dim==3)
        {
            position->box=position->box.scaleAroundCenter(1.0);
        }
        else
        {
            const int ref=2;
            double Z=position->box.center()[ref];
            position->box.p1[ref]=Z;
            position->box.p2[ref]=Z;
        }
        query->getInputPort("position")->writeValue(SharedPtr<Object>(position));
    }

    query->setAccessIndex(0);//<ctc> I think default (0) is fine...

    VisusInfo()<<"adding the dataflow_dataset to the dataflow...";

    dataflow->addNode(dataflow_dataset);
    dataflow->addNode(query);

    this->dataflow->connectNodes(dataflow_dataset,"dataset","dataset",query);

    //only load one level (VisIt doesn't support streaming)
    query->getInputPort("progression")->writeValue(SharedPtr<IntObject>(new IntObject(0)));

    //enable view-dependent data loading
    query->getInputPort("enable_viewdep")->writeValue(SharedPtr<BoolObject>(new BoolObject(true)));

    //fieldname
    query->getInputPort("fieldname")->writeValue(SharedPtr<StringObject>(new StringObject(dataset->default_field.name)));

    VisusInfo()<<"querying the bounds...";

    bounds[0] = dataset->logic_box.to.x-dataset->logic_box.from.x+1;
    bounds[1] = dataset->logic_box.to.y-dataset->logic_box.from.y+1;
    bounds[2] = dataset->logic_box.to.z-dataset->logic_box.from.z+1;

    Box physic_box=dataflow_dataset->getContentPhysicPosition().toAABB();
    fullextents[0] = extents[0] = physic_box.p1.x;
    fullextents[1] = extents[1] = physic_box.p2.x;
    fullextents[2] = extents[2] = physic_box.p1.y;
    fullextents[3] = extents[3] = physic_box.p2.y;
    fullextents[4] = extents[4] = physic_box.p1.z;
    fullextents[5] = extents[5] = physic_box.p2.z;
    
    NdPoint dims(bounds[0],bounds[1],bounds[2],1,1);
    NdBox   exts(NdPoint(extents[0],extents[2],extents[4]),NdPoint(extents[1],extents[3],extents[5]));
    VisusInfo()<<"bounds:  "<<dims.toString();
    VisusInfo()<<"extents: "<<exts.toString();

    VisusInfo()<<"attaching a dummy node";

    DummyNode *dummy=new DummyNode;
    this->dataflow->addNode(dummy);
    this->dataflow->connectNodes(query,"data","data",dummy);

    VisusInfo()<<"end of constructor!";
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
    VisusInfo()<<"destructor?!";
    num_instances--;

    disableSlots(this);

    //call this as soon as possible!
    //emitDestroySignal(); 

    //delete selsApplied; //don't think we own this... (or need to keep it at all)
}

/////////////////////////////////////////////////////////////////////////////
void avtIDXFileFormat::onDataflowInput(DataflowNode *dnode)
{
    //VisusInfo()<<"avtIDXFileFormat::onDataflowInput...";
    if (!dnode)
    {
        VisusAssert(false);
        return;
    }

    //VisusInfo()<<"calling dnode->processInput()...";
    bool ret=dnode->processInput();
    VisusInfo()<<"avtIDXFileFormat::onDataflowInput: dnode->processInput() returned "<<ret;
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
    //VisusInfo()<<"avtIDXFileFormat::GetNTimesteps...";
    int NTimesteps = dataset->time_range.delta()/dataset->time_range.step;
    //VisusInfo()<<"range: "<<dataset->time_range.delta()<<",step: "<<dataset->time_range.step;
    //VisusInfo()<<"\tnum_timesteps="<<NTimesteps;
    return std::max(1,NTimesteps); //<ctc> needs to return at least 1!!
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
    VisusInfo()<<"avtIDXFileFormat::FreeUpResources...";
}

///////////////////////////////////////////////////////////
struct ViewInfo
{
    double   camera[3];
    double   focus[3];
    double   viewUp[3];
    double   viewAngle;
    double   eyeAngle;
    double   parallelScale;
    bool     setScale;
    double   nearPlane;
    double   farPlane;
    double   imagePan[2];
    double   imageZoom;
    bool     orthographic;
    double   shear[3];


    ///////////////////////////////////////////////////////////
    void SetFromView3D(const View3DAttributes &view)
    {
        double    distance;
        double    normal2[3];

        //
        // Calculate a unit length normal.
        //
        distance = sqrt(view.viewNormal[0] * view.viewNormal[0] + view.viewNormal[1] * view.viewNormal[1] +
                        view.viewNormal[2] * view.viewNormal[2]);
        distance = (distance != 0) ? distance : 1.;
        normal2[0] = view.viewNormal[0] / distance;
        normal2[1] = view.viewNormal[1] / distance;
        normal2[2] = view.viewNormal[2] / distance;

        //
        // The view up vector and focal point are the same.  The distance from the
        // camera to the focal point can be calculated from the parallel scale and
        // view angle.  The camera position is then found by moving along the view
        // plane normal from the focal point by the distance.
        //
        viewUp[0] = view.viewUp[0];
        viewUp[1] = view.viewUp[1];
        viewUp[2] = view.viewUp[2];

        focus[0] = view.focus[0];
        focus[1] = view.focus[1];
        focus[2] = view.focus[2];

        eyeAngle = view.eyeAngle;

        distance = view.parallelScale / tan(view.viewAngle * 3.1415926535 / 360.);
        camera[0] = view.focus[0] + normal2[0] * distance;
        camera[1] = view.focus[1] + normal2[1] * distance;
        camera[2] = view.focus[2] + normal2[2] * distance;

        //
        // Orthographic is the opposite of perspective, setScale is always true.
        // It forces vtk to use the parallel scale.
        //
        orthographic  = !view.perspective;
        setScale      = true;
        parallelScale = view.parallelScale;
        viewAngle     = view.viewAngle;

        //
        // The minimum near clipping distance must be adaptive to make good use
        // of the zbuffer.  The distance between the near and far planes seemed
        // like a good choice, another possibility could have been the distance
        // between the camera and focus.  The 5000. is a magic number.  The
        // number should be as large as possible.  10000 would probably also
        // work, but 100000 would start showing z buffering artifacts.
        //
        nearPlane = std::max(view.nearPlane + distance, (view.farPlane - view.nearPlane) / 5000.);
        farPlane = view.farPlane + distance;

        //
        // Set the image pan and image zoom.
        //
        imagePan[0] = -view.imagePan[0];
        imagePan[1] = -view.imagePan[1];
        imageZoom   = view.imageZoom;

        //
        // Set the vew shear.
        //
        shear[0] = view.shear[0];
        shear[1] = view.shear[1];
        shear[2] = view.shear[2];
    }
};

///////////////////////////////////////////////////////////
Frustum* calcFrustum()
{
    //<ctc> todo: handle 2d

    //get window attributes, then use avtViewInfo because it's a more direct mapping
    const View3DAttributes &atts3d=avtCallback::GetCurrentWindowAtts().GetView3D();
    avtView3D v3d;  v3d.SetFromView3DAttributes(&atts3d); //ugh, can't include this either
    avtViewInfo vi; v3d.SetViewInfoFromView(vi); //ugh, can't include it. 
    //ViewInfo vi; vi.SetFromView3D(atts3d);

    double scale[3] = {1,1,1};
    const int *sz=avtCallback::GetCurrentWindowAtts().GetSize();    
    float aspect=(float)sz[0]/(float)sz[1];
    vtkMatrix4x4 *transform = vtkMatrix4x4::New();
    avtWorldSpaceToImageSpaceTransform::CalculateTransform(vi, transform, scale, aspect);

    UniquePtr<Frustum> frustum(new Frustum);
    frustum->setViewport(Viewport(0,0,sz[0],sz[1]));
    frustum->loadProjection(Matrix((double*)(&transform->Element[0])));
    frustum->loadModelview(Matrix::identity());

    Point3d pos(vi.camera[0],vi.camera[2],vi.camera[2]);
    Point3d la(vi.focus[0],vi.focus[1],vi.focus[2]);
    Point3d vup(vi.viewUp[0],vi.viewUp[1],vi.viewUp[2]);
#if 0
    frustum->loadModelview(Matrix::lookAt(pos,la,vup));
    frustum->loadProjection(Matrix::perspective(vi.viewAngle,aspect,vi.nearPlane,vi.farPlane));
#endif
 
    VisusInfo()<<"viewport:   "<<frustum->getViewport().toString();
    VisusInfo()<<"modelview:  "<<frustum->getModelview().toString();
    VisusInfo()<<"projection: "<<frustum->getProjection().toString();

    VisusInfo()<<"created from...";
    VisusInfo()<<"\tcamera.pos: "<<pos;
    VisusInfo()<<"\tcamera.la:  "<<la;
    VisusInfo()<<"\tcamera.vup: "<<vup;
    VisusInfo()<<"\tcamera.fov  "<<vi.viewAngle;
    VisusInfo()<<"\twindow.w:   "<<sz[0];
    VisusInfo()<<"\twindow.h:   "<<sz[1];
    VisusInfo()<<"\tnearPlane:  "<<vi.nearPlane;
    VisusInfo()<<"\tfarPlane:   "<<vi.farPlane;
    
    return frustum.release();
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
    int timestate) 
{
    VisusInfo()<<"avtIDXFileFormat::PopulateDatabaseMetaData, timestate("<<timestate<<")";

    avtMeshType mt = AVT_RECTILINEAR_MESH;
    string mesh_for_this_var;
    mesh_for_this_var.assign("CC_Mesh");
    int block_origin = 0;
    int spatial_dimension = 3;
    int topological_dimension = 3;
    int ndtype;

    //    <ctc> see if I can break in this function when I create a
    //    pseudocolor plot. Maybe there's a way to do this without
    //    resampling. Currently pseudocolor doesn't work, never calls
    //    PopulateDatabaseMetaData again, but where does it call it
    //    the first time? Maybe I can force it to call again...


    //
    // dynamic decomposition
    //
    md->SetFormatCanDoDomainDecomposition(true);

    //
    // Add the mesh. For formats that do their own decomposition, VisIt
    // expects we only advertise a single block.
    //
    int metadata_nblocks = 1;

    //set frustum for view dependent read
    SharedPtr<Frustum> frustum(calcFrustum());
    query->getInputPort("viewdep")->writeValue(frustum);
    query->getInputPort("time")->writeValue(SharedPtr<IntObject>(new IntObject(timestate)));

    //need to get the size of the target volume here, but don't fetch data (we don't know the varname yet!)
    query->getInputPort("position_only")->writeValue(SharedPtr<BoolObject>(new BoolObject(true)));

    this->dataflow->dispatchPublishedMessages(); //<ctc> have to dispatch this time to propagate dataflow_dataset->query connection. It's basically a no-op after the first callde
    this->dataflow->oninput.emitSignal(query);   //<ctc> still need to do this since after the first propagation it will no longer call onInput for the downstream nodes.

    //VisusInfo()<<"now we should have the bounds and extents of the data!";
    NdPoint dims(bounds[0],bounds[1],bounds[2],1,1);
    NdBox   exts(NdPoint(extents[0],extents[2],extents[4]),NdPoint(extents[1],extents[3],extents[5]));
    VisusInfo()<<"PopulateDatabaseMetaData: bounds:  "<<dims.toString();
    VisusInfo()<<"PopulateDatabaseMetaData: extents: "<<exts.toString();

    // Add the mesh, scalars and vectors.
    AddMeshToMetaData(md, mesh_for_this_var, mt, extents, metadata_nblocks, block_origin, spatial_dimension, topological_dimension, bounds);
    std::vector<string> &fieldnames = dataset->fieldnames;
    
    for (int i = 0; i < (int) fieldnames.size(); i++)
    {
        Field field = dataset->getFieldByName(fieldnames[i]);
        avtCentering cent = AVT_ZONECENT;
        ndtype=1;
        if (field.dtype.isVector())
            ndtype=field.dtype.ncomponents();
        if (ndtype == 1)
            AddScalarVarToMetaData(md, fieldnames[i], mesh_for_this_var, cent);
        else
            AddVectorVarToMetaData(md, fieldnames[i], mesh_for_this_var, cent, ndtype);
    }

    //
    // We want to set the LOD property of the Mesh Meta Data. Since we only
    // have one mesh we can assume that md->GetMeshes(0) points to the
    // avtMeshMetaData Object created with the `AddMeshToMetaData' helper.
    //
    md->GetMeshes(0).LODs = 16;//(dataset->maxh - 15) / 3;
    //resolution <ctc> tie LOD (resolution) to "quality" port.
    //query->getInputPort("quality")->writeValue(resolution-8);
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

    //resolution=avtCallback::idx_get_resolution_hack();
    //cerr<<"resolution is at "<<resolution<<endl;

    VisusInfo() << "avtIDXFileFormat::GetMesh: timestate("<<timestate<<") domain("<<domain<<") meshname("<<meshname<<") resolution("<<resolution<<")";

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
    
    dims[0] = bounds[0]+1; //visit is so weird... if I have a x by y array, it wants me to say it's x+1 by y+1 :P
    dims[1] = bounds[1]+1;
    dims[2] = bounds[2]+1;

    {
        NdPoint dims(bounds[0],bounds[1],bounds[2],1,1);
        NdBox   exts(NdPoint(extents[0],extents[2],extents[4]),NdPoint(extents[1],extents[3],extents[5]));
        VisusInfo()<<"GetMesh: bounds:  "<<dims.toString();
        VisusInfo()<<"GetMesh: extents: "<<exts.toString();
    }

    // Point3d spacing((fullextents[1]-fullextents[0]+1)/bounds[0],(fullextents[3]-fullextents[2]+1)/bounds[1],(fullextents[5]-fullextents[4]+1)/bounds[2]);
    // if (resolution==1) //see note in publish above
    Point3d spacing((extents[1]-extents[0]+1)/bounds[0],(extents[3]-extents[2]+1)/bounds[1],(extents[5]-extents[4]+1)/bounds[2]);

    rgrid->SetDimensions(dims[0], dims[1], dims[2]);
    VisusInfo()<<"avtIDXFileFormat::GetMesh() returning <"<<dims[0]<<"x"<<dims[1]<<"x"<<dims[2]<<"> mesh";//for resolution "<<resolution<<" (resReduction="<<resReduction<<")";

    coordsX = vtkFloatArray::New();
    coordsX->SetNumberOfTuples(dims[0]);
    arrayX = (float *) coordsX->GetVoidPointer(0);
    for (int i = 0; i < dims[0]; i++)
        //arrayX[i] = i * resReduction;
        arrayX[i] = (extents[0]+i)*spacing.x;
    //arrayX[dims[0]-1] = slice_box.to[0] + 1.; //<ctc> maybe still need something like this
    rgrid->SetXCoordinates(coordsX);

    coordsY = vtkFloatArray::New();
    coordsY->SetNumberOfTuples(dims[1]);
    arrayY = (float *) coordsY->GetVoidPointer(0);
    for (int i = 0; i < dims[1]; i++)
        //arrayY[i] = i * resReduction;
        arrayY[i] = (extents[2]+i)*spacing.y;
    //arrayY[dims[1]-1] = slice_box.to[1] + 1.; //<ctc> maybe still need something like this
    rgrid->SetYCoordinates(coordsY);

    coordsZ = vtkFloatArray::New();
    coordsZ->SetNumberOfTuples(dims[2]);
    arrayZ = (float *) coordsZ->GetVoidPointer(0);
    for (int i = 0; i < dims[2]; i++)
        //arrayZ[i] = slice_box.from[2] + i * resReduction;
        arrayZ[i] = (extents[4]+i)*spacing.z;
    //arrayZ[dims[2]-1] = slice_box.to[2] + 1.; //<ctc> maybe still need something like this
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
    //resolution=avtCallback::idx_get_resolution_hack();
    //VisusInfo()<<"resolution is at "<<resolution<<endl;
    VisusInfo() << "avtIDXFileFormat::GetVar: timestate("<<timestate<<") domain("<<domain<<") varname("<<varname<<") resolution("<<resolution<<")";
    string name(varname);
    NdBox slice_box = dataset->logic_box;

    query->getInputPort("fieldname")->writeValue(SharedPtr<StringObject>(new StringObject(varname)));
    query->getInputPort("position_only")->writeValue(SharedPtr<BoolObject>(new BoolObject(false)));
    query->getInputPort("time")->writeValue(SharedPtr<IntObject>(new IntObject(timestate)));

    this->dataflow->oninput.emitSignal(query);
    //this->dataflow->dispatchPublishedMessages();

    VisusInfo()<<"query started, waiting for data...";
    Clock t0(Clock::now());

    //Ack! What if the onDataflowInput that calls query->processInput returns false?!
    // need to fail gracefully or... should that never happen?

    //wait for the data to arrive.
    haveData=false;
    while (!haveData) ;//wait

    Clock::timestamp_t msec=t0.msec();
    VisusInfo()<<msec<<"ms to fetch data, now send it to VisIt.";
    
    Field field = dataset->getFieldByName(varname);
    NdPoint dims(bounds[0],bounds[1],bounds[2],1,1);
    long ntuples = dims.innerProduct();
    int ncomponents=1;
    if (field.dtype==DType::UINT8)
    {
        vtkUnsignedCharArray*rv = vtkUnsignedCharArray::New();
        rv->SetNumberOfComponents(ncomponents); //<ctc> eventually handle vector data, since visit can actually render it!
        data->unmanaged=true; //giving the data to VisIt which will delete it when it's no longer needed
        rv->SetArray((unsigned char*)data->c_ptr(),ncomponents*ntuples,1/*delete when done*/,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        return rv;
    }
    if (field.dtype==DType::UINT16)
    {
        vtkUnsignedShortArray *rv = vtkUnsignedShortArray::New();
        rv->SetNumberOfComponents(ncomponents);
        data->unmanaged=true;
        rv->SetArray((unsigned short*)data->c_ptr(),ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        return rv;
    }
    if (field.dtype==DType::UINT32)
    {
        vtkUnsignedIntArray *rv = vtkUnsignedIntArray::New();
        rv->SetNumberOfComponents(ncomponents);
        data->unmanaged=true;
        rv->SetArray((unsigned int*)data->c_ptr(),ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        return rv;
    }
    if (field.dtype==DType::UINT32)
    {
        vtkUnsignedLongArray *rv = vtkUnsignedLongArray::New();
        rv->SetNumberOfComponents(ncomponents);
        data->unmanaged=true;
        rv->SetArray((unsigned long*)data->c_ptr(),ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        return rv;
    }
    if (field.dtype==DType::INT8)
    {
        vtkCharArray*rv = vtkCharArray::New();
        rv->SetNumberOfComponents(ncomponents);
        data->unmanaged=true;
        rv->SetArray((char*)data->c_ptr(),ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        return rv;
    }
    if (field.dtype==DType::INT16)
    {
        vtkShortArray *rv = vtkShortArray::New();
        rv->SetNumberOfComponents(ncomponents);
        data->unmanaged=true;
        rv->SetArray((short*)data->c_ptr(),ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        return rv;
    }
    if (field.dtype==DType::INT32)
    {
        vtkIntArray *rv = vtkIntArray::New();
        rv->SetNumberOfComponents(ncomponents);
        data->unmanaged=true;
        rv->SetArray((int*)data->c_ptr(),ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        return rv;
    }
    if (field.dtype==DType::INT64)
    {
        vtkLongArray *rv = vtkLongArray::New();
        rv->SetNumberOfComponents(ncomponents);
        data->unmanaged=true;
        rv->SetArray((long*)data->c_ptr(),ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        return rv;
    }
    if (field.dtype==DType::FLOAT32)
    {
        vtkFloatArray *rv = vtkFloatArray::New();
        rv->SetNumberOfComponents(ncomponents);
        data->unmanaged=true;
        rv->SetArray((float*)data->c_ptr(),ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        return rv;
    }
    if (field.dtype==DType::FLOAT64)
    {
        vtkDoubleArray *rv = vtkDoubleArray::New();
        rv->SetNumberOfComponents(ncomponents);
        data->unmanaged=true;
        rv->SetArray((double*)data->c_ptr(),ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        return rv;
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
    //<ctc> this function gets called later in the pipeline, not soon
    //enough to prevent incorrect dims from being set in
    //PopulateDatabaseMetaData. So we use another strategy:
    //avtCallback.

    this->selList = sels;
    this->selsApplied = applied;

    //VisusInfo() << "avtIDXFileFormat::RegisterDataSelections" << endl;
    for(size_t i=0; i < sels.size(); ++i)
    {
        if(strcmp(sels[i]->GetType(), "avtResolutionSelection") == 0)
        {
            const avtResolutionSelection* sel = static_cast<const avtResolutionSelection*>(*sels[i]);
            VisusInfo()<<"new resolution: "<<sel->resolution()<<", (old resolution: "<<resolution<<")";
            if (resolution!=sel->resolution())
            {
                ;
                resolution = sel->resolution(); //<ctc> with hack function above, just don't mess with it.
                //ClearCache();
            }
            else
            {
                (*applied)[i] = true; //<ctc> this doesn't really seem to do anything.
            }
        }
        else if (strcmp(sels[i]->GetType(), "avtFrustumSelection") == 0)
        {
            const avtFrustumSelection* sel = static_cast<const avtFrustumSelection*>(*sels[i]);
            //VisusInfo()<<"new frustum\n";
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
            VisusInfo()<<"avtIDXFileFormat::RegisterDataSelections: unhandled selection "<<sels[i]->GetType()<<"!\n";
    }
}
