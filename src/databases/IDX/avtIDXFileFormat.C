/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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

// *************************************************************************
//                           avtIDXFileFormat.C
// *************************************************************************

#include <avtIDXFileFormat.h>

#include <string>
#include <sstream>

#include <vtkFloatArray.h>
#include <vtkTypeFloat32Array.h>
#include <vtkDoubleArray.h>
#include <vtkCharArray.h>
#include <vtkShortArray.h>
#include <vtkIntArray.h>
#include <vtkCellType.h>
#include <vtkCellData.h>
#include <vtkLongArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnsignedShortArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnsignedLongArray.h>
#include <vtkSmartPointer.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkRenderWindow.h>
#include <vtkMatrix4x4.h>
#include <vtkXMLDataElement.h>
#include <vtkXMLDataParser.h>

#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <avtResolutionSelection.h>
#include <avtStructuredDomainBoundaries.h>
#include <avtVariableCache.h>
#include <avtParallel.h>
#include <avtDatabaseMetaData.h>
#include <avtMultiresSelection.h>
#include <avtCallback.h>
#include <avtView2D.h>
#include <avtView3D.h>

#include <DBOptionsAttributes.h>
#include <Expression.h>

#include <InvalidVariableException.h>
#include <dirent.h>

#ifdef PARALLEL
#include <avtParallel.h>
#endif

typedef std::string String;

using namespace VisusSimpleIO;

int        cint   (String s) {int    value;std::istringstream iss(s);iss>>value;return value;}
float      cfloat (String s) {float  value;std::istringstream iss(s);iss>>value;return value;}
double     cdouble(String s) {double value;std::istringstream iss(s);iss>>value;return value;}

void avtIDXFileFormat::loadBalance(){
    
  //std::cout << "Load balancing";
    
    int maxdir = 0; // largest extent axis
    int maxextent = 0;
    int maxbox = 0;
    
    for(int i=0; i < boxes.size(); i++){
        SimpleBox& box = boxes.at(i);
        
        for(int j=0; j < 3; j++){
            int extent = box.p2[j]-box.p1[j];
            if(extent > maxextent){
                maxdir = j;
                maxextent = extent;
                maxbox = i;
            }
        }
    }
    
    //std::cout << "max dir " << maxdir << " max extent " << maxextent << " box " << maxbox;
    
    int total_extent = 0;
    int avg_ext = 0;
    
    for(int i=0; i < boxes.size(); i++){
        SimpleBox& box = boxes.at(i);
        
        total_extent += box.p2[maxdir];
    }
    //total_extent++; // DIM
    avg_ext = total_extent / nprocs;
    int res_ext = (total_extent % nprocs);
    
    //  std::cout << "tot ext " << total_extent << " avg ext " << avg_ext << " res ext " << res_ext;
    
    std::vector<SimpleBox> newboxes;
    std::vector<SimpleBox> newphyboxes;
    
    for(int i=0; i < boxes.size(); i++){
        SimpleBox& box = boxes.at(i);
        
        int loc_avg_ext = box.p2[maxdir] - box.p1[maxdir];
        int loc_res = 0;
        
        if(loc_avg_ext > avg_ext){
            loc_res = loc_avg_ext % avg_ext;
            loc_avg_ext = avg_ext;
        }
        
        //        std::cout << "local avg ext " << loc_avg_ext << " local res " << loc_res;
        
        int part_p1 = box.p1[maxdir];
        int part_p2 = box.p1[maxdir] + loc_avg_ext;
        
        SimplePoint3d p1(box.p1);
        SimplePoint3d p2(box.p2);

        SimplePoint3d log2phy;
        SimpleBox& phybox = phyboxes.at(i);
      
      std::cout << "Old box p1: " << p1 << " p2: "<< p2 << " phy " << phybox.p1 << " p2 " << phybox.p2<< std::endl;
      
        for (int k = 0; k < dim; k++){
          log2phy[k] = (phybox.p2[k] - phybox.p1[k])/(box.p2[k] - box.p1[k] + 1);
       }
      
//      std::cout << "log2phy " <<log2phy << std::endl;
      
        while(part_p2 <= box.p2[maxdir]){
            
            p1[maxdir] = part_p1;
            p2[maxdir] = (part_p2 < box.p2[maxdir]) ? part_p2+1 : part_p2;
           
            SimpleBox newbox(p1,p2);
            newboxes.push_back(newbox);
          
            SimpleBox newphybox;
          
            for (int k = 0; k < dim; k++){
              newphybox.p1[k] = newbox.p1[k] * log2phy[k];
              newphybox.p2[k] = newbox.p2[k] * log2phy[k];
            }
          
            newphyboxes.push_back(newphybox);
            
            std::cout << "New box p1: " << p1 << " p2: "<< p2 << " phy " << newphybox.p1 << " p2 " << newphybox.p2<< std::endl;
          
            part_p1 += loc_avg_ext;
            part_p2 += loc_avg_ext;
           
        }
        
        if(loc_res > 0){
            SimpleBox& boxres = newboxes.at(newboxes.size()-1);
            boxres.p2[maxdir] += loc_res-1;
          
            SimpleBox& phyboxres = newphyboxes.at(newphyboxes.size()-1);
            phyboxres.p2[maxdir] += loc_res*log2phy[maxdir];
          
              std::cout << "Residual " << loc_res-1 <<" added to box "<< newboxes.size()-1 <<" p1: " << boxres.p1 << " p2: "<< boxres.p2;
        }

    }
    
    boxes.swap(newboxes);
    phyboxes.swap(newphyboxes);
  
    std::cout << "Total number of boxes/domains: " << boxes.size() << std::endl;
    std::cout << "----------Boxes----------" << std::endl;
    for(int i=0; i< boxes.size(); i++){
        std::cout << i << " = "<< boxes.at(i).p1 << " , " << boxes.at(i).p2 << " phy: "
          << phyboxes.at(i).p1 << " , " << phyboxes.at(i).p2 << std::endl;
    }
    std::cout << "-------------------------" << std::endl;
    
}

template <typename Type>
Type* avtIDXFileFormat::convertComponents(const unsigned char* src, int src_ncomponents, int dst_ncomponents, long long totsamples){
    int n=src_ncomponents;
    int m=dst_ncomponents;
    int ncomponents=std::min(m,n);
    
    std::cout << "n " << n << " m " << m << " ncomp " << ncomponents;
    
    Type* dst = (Type*)calloc(totsamples*m, sizeof(Type));
    std::cout << "new vector size " << totsamples << std::endl;
    
    //for each component...
    for (int C=0;C<ncomponents;C++)
    {
        Type* src_p=((Type*)src)+C;
        Type* dst_p=((Type*)dst)+C;
        for (long long I=0; I<totsamples; I++,src_p+=n,dst_p+=m)
        {
            *dst_p=*src_p;
            
            //std::cout << "c " << C << " I " <<I << std::endl;
        }
        
    }
    
    std::cout << "data converted " << std::endl;
    
    return dst;
}

// TODO consider the physical box
void avtIDXFileFormat::calculateBoundsAndExtents(){
    
    // TODO deallocate this stuff
    for(int i=0; i< boxes.size(); i++){
        SimpleBox& box = boxes.at(i);
        int* my_bounds = new int[3];
            
        my_bounds[0] = box.p2.x-box.p1.x;
        my_bounds[1] = box.p2.y-box.p1.y;
        my_bounds[2] = box.p2.z-box.p1.z;
        
        boxes_bounds.push_back(my_bounds);
        
    }
    
}

void avtIDXFileFormat::createBoxes(){
    
    size_t found = dataset_filename.find_last_of("/\\");
    String folder = dataset_filename.substr(0,found);
    
    String upsfilename = "noupsfile.ups";
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (folder.c_str())) != NULL) {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL) {
            String name(ent->d_name);
            if(name.substr(name.find_last_of(".") + 1) == "ups"){
                upsfilename = name;
                std::cout<< ".ups file found " << upsfilename << std::endl;
                upsfilename = folder + "/" +upsfilename;
                break;
            }
        }
        closedir (dir);
    } else {
        std::cout<< "No .ups file found" << std::endl;
    }
    
    vtkSmartPointer<vtkXMLDataParser> parser = vtkSmartPointer<vtkXMLDataParser>::New();
    
    upsfilename.replace(upsfilename.end()-3, upsfilename.end(),"ups");
    
    parser->SetFileName(upsfilename.c_str());
    if (!parser->Parse())
    {
        std::cout<< "No .ups file found" << std::endl;
        multibox = false;
        
        std::cout << "Single-box mode" << std::endl;
    }else{
        multibox = true;
        std::cout << "Multi-box mode" << std::endl;
    }
    
    if(multibox){
        vtkXMLDataElement *root = parser->GetRootElement();
        vtkXMLDataElement *level = root->FindNestedElementWithName("Grid")->FindNestedElementWithName("Level");
        int nboxes = level->GetNumberOfNestedElements();
        
        std::cout << "Found " << nboxes << " boxes" << std::endl;
        
        for(int i=0; i < nboxes; i++){
            
            vtkXMLDataElement *xmlbox = level->GetNestedElement(i);
            String lower(xmlbox->FindNestedElementWithName("lower")->GetCharacterData());
            String upper(xmlbox->FindNestedElementWithName("upper")->GetCharacterData());
            String extra_cells(xmlbox->FindNestedElementWithName("extraCells")->GetCharacterData());
            String resolution(xmlbox->FindNestedElementWithName("resolution")->GetCharacterData());
            
            lower = lower.substr(1,lower.length()-2);
            upper = upper.substr(1,upper.length()-2);
            
            // TODO Extra cells managements (add flag on the plugin)
            extra_cells = extra_cells.substr(1,extra_cells.length()-2);
            resolution = resolution.substr(1,resolution.length()-2);
            
            //std::cout<< "lower " << lower << " upper " << upper;
            
            SimplePoint3d p1phy, p1log;
            SimplePoint3d p2phy, p2log;
            int eCells[3];
            int resdata[3];
            double phy2log[3];
            
            std::stringstream ress(resolution);
            std::stringstream ss1(lower);
            std::stringstream ss2(upper);
            std::stringstream ssSpace(extra_cells);
            std::string p1s, p2s, espace, res;
            for (int k=0; k < 3; k++){
                std::getline(ss1, p1s, ',');
                std::getline(ss2, p2s, ',');
                std::getline(ssSpace, espace, ',');
                std::getline(ress, res, ',');
                
                eCells[k] = cint(espace);
                resdata[k] = cint(res);
                
                p1phy[k] = cfloat(p1s);
                p2phy[k] = cfloat(p2s);
                
                phy2log[k] = (p2phy[k]-p1phy[k])/resdata[k];
                p2phy[k] += phy2log[k];
                
                p1log[k] = p1phy[k] / phy2log[k];
                p2log[k] = p1log[k] + resdata[k] +1;
                
                if (use_extracells)
                    p2log[k] += eCells[k];
            }
        
            std::cout <<"Read box phy: p1 " << p1phy << " p2 "<< p2phy << std::endl;
            std::cout <<"     box log: p1 " << p1log << " p2 "<< p2log << std::endl;
            
            phyboxes.push_back(SimpleBox(p1phy, p2phy));
            physicalBox = physicalBox.getUnion((const SimpleBox)phyboxes.at(i));
            boxes.push_back(SimpleBox(p1log,p2log));
            
        }
        
    }
    else{
        boxes.push_back(reader->getLogicBox());
        physicalBox = reader->getLogicBox();
        phyboxes.push_back(physicalBox);
    }

}

void avtIDXFileFormat::createTimeIndex(){
    vtkSmartPointer<vtkXMLDataParser> parser = vtkSmartPointer<vtkXMLDataParser>::New();
    size_t found = dataset_filename.find_last_of("/\\");
    String folder = dataset_filename.substr(0,found);
    
    String udafilename = folder + "/index.xml";
    
    parser->SetFileName(udafilename.c_str());
    if (!parser->Parse()){
        std::cout<< "No index.xml file found" << udafilename << std::endl;
        
        std::vector<double> times = reader->getTimes();
        
        for(int i=0; i< times.size(); i++)
            timeIndex.push_back(times.at(i));
        
        return;
    }
    else{
        
        std::cout << "Found index.xml file" << std::endl;

        vtkXMLDataElement *root = parser->GetRootElement();
        vtkXMLDataElement *level = root->FindNestedElementWithName("timesteps");
        int ntimesteps = level->GetNumberOfNestedElements();
        
        std::cout << "Found " << ntimesteps << " timesteps" << std::endl;
        
        for(int i=0; i < ntimesteps; i++){
            
            vtkXMLDataElement *xmltime = level->GetNestedElement(i);
            String timestr(xmltime->GetAttribute("time"));
            std::cout << "time " << timestr << std::endl;
            
            double time = cdouble(timestr);
            
            timeIndex.push_back(time);
        }
    }
    
}

// ****************************************************************************
//  Method: avtIDXFileFormat constructor
//
//  Programmer: bremer5 -- generated by xml2avt
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

avtIDXFileFormat::avtIDXFileFormat(const char *filename, DBOptionsAttributes* attrs)
: avtMTMDFileFormat(filename)
{
    for (int i=0; attrs!=0 && i<attrs->GetNumberOfOptions(); ++i) {
        if (attrs->GetName(i) == "Big Endian") {
            reverse_endian = attrs->GetBool("Big Endian");
        }
        else if (attrs->GetName(i) == "Use extra cells") {
            use_extracells = attrs->GetBool("Use extra cells");
        }
    }
    
    if(use_extracells)
        std::cout << "Using extra cells" << std::endl;
    else
        std::cout << "Not using extra cells" << std::endl;
        
    if(reverse_endian)
        std::cout << "Using Big Endian";
    else
        std::cout << "Using Little Endian";
    
#ifdef PARALLEL
    rank = PAR_Rank();
    nprocs = PAR_Size();
#else
    rank = 0;
    nprocs = 1;
#endif
    
    std::cout << "~~~PROC " << rank << " / " << nprocs << std::endl;
  
    reader = new SimpleIO();
  
    if (!reader->openDataset(filename))
    {
        std::cout <<"could not load "<<filename << std::endl;
        return;
    }
    
    dataset_filename = filename;
    
//    std::cout <<"dataset loaded";
    dim = reader->getDimension(); //<ctc> //NOTE: it doesn't work like we want. Instead, when a slice (or box) is added, the full data is read from disk then cropped to the desired subregion. Thus, I/O is never avoided.
    
    // TODO (if necessary) read only with rank 0 and then broadcast to the other processors
    createBoxes();
    createTimeIndex();
    loadBalance();
    calculateBoundsAndExtents();
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
    std::cout<<"(avtIDXFileFormat destructor)" << std::endl;

    for(int i=0; i < boxes_bounds.size(); i++)
        if(boxes_bounds.at(i) != NULL)
            delete [] boxes_bounds.at(i);
  
    if(reader != NULL)
      delete reader;
  
}

// ****************************************************************************
//  Method: avtIDXFileFormat::GetNTimesteps
//
//  Purpose:
//      Tells the rest of the code how many timesteps there are in this file.
//
//  Programmer: spetruzza, bremer5
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

int
avtIDXFileFormat::GetNTimesteps(void)
{
    return reader->getNTimesteps();
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
//  Programmer: spetruzza, bremer5
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

void
avtIDXFileFormat::FreeUpResources(void)
{
    std::cout<<"avtIDXFileFormat::FreeUpResources..." << std::endl;
    //<ctc> todo... something (is destructor also called?)
}


//bool avtIDXFileFormat::CanCacheVariable(const char *var)
//{
//    return false;
//}

// ****************************************************************************
//  Method: avtIDXFileFormat::PopulateDatabaseMetaData
//
//  Purpose:
//      This database meta-data object is like a table of contents for the
//      file.  By populating it, you are telling the rest of VisIt what
//      information it can request from you.
//
//  Programmer: spetruzza, bremer5
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

void
avtIDXFileFormat::PopulateDatabaseMetaData(avtDatabaseMetaData *md,
    int timestate) 
{
    std::cout << rank << ": Meta data" << std::endl;

    avtMeshMetaData *mesh = new avtMeshMetaData;
    mesh->name = "Mesh";
    
    mesh->meshType = AVT_RECTILINEAR_MESH;
    
    mesh->numBlocks = boxes.size();
    mesh->blockOrigin = 0;
    mesh->LODs = reader->getMaxResolution();
    mesh->spatialDimension = dim;
    mesh->topologicalDimension = dim;
    
    mesh->blockTitle = "box";
    mesh->blockPieceName = "box";
    
    // Set bounds and extents for SLIVR rendering
    // TODO use the physical box (logic_to_physic)
    mesh->hasSpatialExtents = true;

    mesh->minSpatialExtents[0] = physicalBox.p1.x;
    mesh->maxSpatialExtents[0] = physicalBox.p2.x;
    mesh->minSpatialExtents[1] = physicalBox.p1.y;
    mesh->maxSpatialExtents[1] = physicalBox.p2.y;
    mesh->minSpatialExtents[2] = physicalBox.p1.z;
    mesh->maxSpatialExtents[2] = physicalBox.p2.z;

    SimpleBox logicBox = reader->getLogicBox();
    mesh->hasLogicalBounds = true;
    mesh->logicalBounds[0] = logicBox.p2.x - logicBox.p1.x;
    mesh->logicalBounds[1] = logicBox.p2.y - logicBox.p1.y;
    mesh->logicalBounds[2] = logicBox.p2.z - logicBox.p1.z;
    
    md->Add(mesh);
    
    std::cout << rank << ": Added mesh";

    const std::vector<SimpleField>& fields = reader->getFields();
    
    int ndtype;
    for (int i = 0; i < (int) fields.size(); i++)
    {
        const SimpleField& field = fields[i];
        
        if (!field.isVector)
            md->Add(new avtScalarMetaData(field.name,mesh->name,AVT_ZONECENT));
        else
            md->Add(new avtVectorMetaData(field.name,mesh->name,AVT_ZONECENT, field.ncomponents));
    }
    
    std::cout << rank << ": Added fields";
        
    avtRectilinearDomainBoundaries *rdb =
    new avtRectilinearDomainBoundaries(true);
    rdb->SetNumDomains(phyboxes.size());
    
    for (long long i = 0 ; i < phyboxes.size() ; i++)
    {
        int extents[6];
        extents[0] = phyboxes.at(i).p1.x;
        extents[1] = phyboxes.at(i).p2.x;
        extents[2] = phyboxes.at(i).p1.y;
        extents[3] = phyboxes.at(i).p2.y;
        extents[4] = phyboxes.at(i).p1.z;
        extents[5] = phyboxes.at(i).p2.z;
        
        rdb->SetIndicesForRectGrid(i, extents);
    }
    rdb->CalculateBoundaries();
    
    std::cout << rank << ": Calculated boundaries";
    
    void_ref_ptr vr = void_ref_ptr(rdb,
                                   avtStructuredDomainBoundaries::Destruct);
    cache->CacheVoidRef("any_mesh",                  AUXILIARY_DATA_DOMAIN_BOUNDARY_INFORMATION, -1, -1, vr);

    return;
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
//  Programmer: spetruzza
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

vtkDataSet *
avtIDXFileFormat::GetMesh(int timestate, int domain, const char *meshname)
{
    //std::cout<< rank << ": start getMesh "<< meshname << " domain " << domain << std::endl;
    
    SimpleBox slice_box;
    
    int* my_bounds = NULL;
    int* my_extents = NULL;

    slice_box = phyboxes.at(domain);
    my_bounds = boxes_bounds.at(domain);
    
    vtkRectilinearGrid *rgrid = vtkRectilinearGrid::New();
    int dims[3];
    float *arrayX;
    float *arrayZ;
    float *arrayY;
    vtkFloatArray *coordsX;
    vtkFloatArray *coordsY;
    vtkFloatArray *coordsZ;
    
    int my_dims[3];
    
    my_dims[0] = my_bounds[0]+1;
    my_dims[1] = my_bounds[1]+1;
    my_dims[2] = my_bounds[2]+1;
    
//    std::cout << rank << ": dims " << my_dims[0] << " " << my_dims[1] << " " << my_dims[2] << std::endl;
    
    rgrid->SetDimensions(my_dims[0], my_dims[1], my_dims[2]);
    
    coordsX = vtkFloatArray::New();
    coordsX->SetNumberOfTuples(my_dims[0]);
    arrayX = (float *) coordsX->GetVoidPointer(0);
    
    float steps[3];
    
    for (int i = 0; i < dim; i++){
        steps[i] = (slice_box.p2[i] - slice_box.p1[i])/my_dims[i];
    }
    
    for (int i = 0; i < my_dims[0]; i++)
        arrayX[i] = slice_box.p1.x + i*steps[0];//+i;
    rgrid->SetXCoordinates(coordsX);
    
    coordsY = vtkFloatArray::New();
    coordsY->SetNumberOfTuples(my_dims[1]);
    arrayY = (float *) coordsY->GetVoidPointer(0);
    for (int i = 0; i < my_dims[1]; i++)
        arrayY[i] = slice_box.p1.y + i*steps[1];//+i;
    rgrid->SetYCoordinates(coordsY);
    
    coordsZ = vtkFloatArray::New();
    coordsZ->SetNumberOfTuples(my_dims[2]);
    arrayZ = (float *) coordsZ->GetVoidPointer(0);
    for (int i = 0; i < my_dims[2]; i++)
        arrayZ[i] = slice_box.p1.z + i*steps[2];//+i;
    rgrid->SetZCoordinates(coordsZ);
    
    std::cout << "end mesh";
    
    return rgrid;
    
}

void
avtIDXFileFormat::GetCycles(std::vector<int> &cycles)
{
    for(int i = 0; i < reader->getNTimesteps(); ++i)
        cycles.push_back(i);
}

void
avtIDXFileFormat::GetTimes(std::vector<double> &times)
{
    times.swap(timeIndex);
//    std::map<int, double> m;
//    std::vector<double> v;
//    for(std::map<int,double>::iterator it = m.begin(); it != m.end(); ++it) {
//        times.push_back(it->second);
//    }
    
//    std::vector<double> tsteps = reader->getTimes();
//    times.swap(tsteps);
}

vtkDataArray* avtIDXFileFormat::queryToVtk(int timestate, int domain, const char *varname){
    
    const SimpleBox& my_box = boxes.at(domain);
    
    unsigned char* data = reader->getData(my_box, timestate, varname, reverse_endian);
    
    if(data == NULL){
        std::cout << " NO DATA " << std::endl;
        return NULL;
    }
    
    SimpleField field = reader->getCurrField();
    SimpleDTypes type = field.type;
    
    int* my_bounds = boxes_bounds.at(domain);
    int ztuples = (dim == 2) ? 1 : (my_bounds[2]);
    long long ntuples = (my_bounds[0])*(my_bounds[1])*ztuples;
    
    int ncomponents = 1;
    
    bool isVector = field.isVector;
    
    if(isVector)
        ncomponents = 3; // Visit wants 3 components vectors
    
    long long ntotal = ntuples * ncomponents;
    
//    if( data != NULL)
//         std::cout<< rank << ": size data bytes " << data->c_size() << std::endl;
//    
//    std::cout << rank << ": size array " << ncomponents*ntuples << std::endl;
  
    if(type == VisusSimpleIO::UINT8){
        vtkUnsignedCharArray*rv = vtkUnsignedCharArray::New();
        rv->SetNumberOfComponents(ncomponents); //<ctc> eventually handle vector data, since visit can actually render it!
        
        if(isVector && dim < 3){
            
            unsigned char* newdata = convertComponents<unsigned char>(data, field.ncomponents, 3, ntuples);
            rv->SetArray((unsigned char*)newdata,ncomponents*ntuples,1,vtkDataArrayTemplate<unsigned char>::VTK_DATA_ARRAY_FREE);
            
            delete data;
        }
        else
            rv->SetArray((unsigned char*)data,ncomponents*ntuples,1/*delete when done*/,vtkDataArrayTemplate<unsigned char>::VTK_DATA_ARRAY_FREE);
        return rv;
    }
    else if(type == VisusSimpleIO::UINT16){
    
        vtkUnsignedShortArray *rv = vtkUnsignedShortArray::New();
        rv->SetNumberOfComponents(ncomponents);
        
        if(isVector && dim < 3){
            
            unsigned short* newdata = convertComponents<unsigned short>(data, field.ncomponents, 3, ntuples);
            rv->SetArray((unsigned short*)newdata,ncomponents*ntuples,1,vtkDataArrayTemplate<unsigned short>::VTK_DATA_ARRAY_FREE);
            
            delete data;
        }
        else
            rv->SetArray((unsigned short*)data,ncomponents*ntuples,1,vtkDataArrayTemplate<unsigned short>::VTK_DATA_ARRAY_FREE);
        
//        if(reverse_endian){
//            unsigned short *buff = (unsigned short *) rv->GetVoidPointer(0);
//            for (long long i = 0 ; i < ntotal ; i++)
//            {
//                int tmp;
//                int16_Reverse_Endian(buff[i], (unsigned char *) &tmp);
//                buff[i] = tmp;
//            }
//        }
    
        return rv;
    }
    else if(type == VisusSimpleIO::UINT32){
        vtkUnsignedIntArray *rv = vtkUnsignedIntArray::New();
        rv->SetNumberOfComponents(ncomponents);
        
        if(isVector && dim < 3){
            
            unsigned int* newdata = convertComponents<unsigned int>(data, field.ncomponents, 3, ntuples);
            rv->SetArray((unsigned int*)newdata,ncomponents*ntuples,1,vtkDataArrayTemplate<unsigned int>::VTK_DATA_ARRAY_FREE);
            
            delete data;
        }
        else
            rv->SetArray((unsigned int*)data,ncomponents*ntuples,1,vtkDataArrayTemplate<unsigned int>::VTK_DATA_ARRAY_FREE);
        
//        if(reverse_endian){
//            unsigned int *buff = (unsigned int *) rv->GetVoidPointer(0);
//            for (long long i = 0 ; i < ntotal ; i++)
//            {
//                int tmp;
//                int32_Reverse_Endian(buff[i], (unsigned char *) &tmp);
//                buff[i] = tmp;
//            }
//        }
      
        return rv;
    }
    else if(type == VisusSimpleIO::INT8){
        vtkCharArray*rv = vtkCharArray::New();
        rv->SetNumberOfComponents(ncomponents);
        
        if(isVector && dim < 3){
            
            char* newdata = convertComponents<char>(data, field.ncomponents, 3, ntuples);
            rv->SetArray((char*)newdata,ncomponents*ntuples,1,vtkDataArrayTemplate<char>::VTK_DATA_ARRAY_FREE);
            
            delete data;
        }
        else
            rv->SetArray((char*)data,ncomponents*ntuples,1,vtkDataArrayTemplate<char>::VTK_DATA_ARRAY_FREE);
        
        return rv;
    }
    else if(type == VisusSimpleIO::INT16){
        vtkShortArray *rv = vtkShortArray::New();
        rv->SetNumberOfComponents(ncomponents);
        
        if(isVector && dim < 3){
            
            short* newdata = convertComponents<short>(data, field.ncomponents, 3, ntuples);
            rv->SetArray((short*)newdata,ncomponents*ntuples,1,vtkDataArrayTemplate<short>::VTK_DATA_ARRAY_FREE);
            
            delete data;
        }
        else
            rv->SetArray((short*)data,ncomponents*ntuples,1,vtkDataArrayTemplate<short>::VTK_DATA_ARRAY_FREE);
        
//        if(reverse_endian){
//            short *buff = (short *) rv->GetVoidPointer(0);
//            for (long long i = 0 ; i < ntotal ; i++)
//            {
//                int tmp;
//                int16_Reverse_Endian(buff[i], (unsigned char *) &tmp);
//                buff[i] = tmp;
//            }
//        }
      
        return rv;
    }
    else if(type == VisusSimpleIO::INT32){
        vtkIntArray *rv = vtkIntArray::New();
        rv->SetNumberOfComponents(ncomponents);
        
        if(isVector && dim < 3){
            
            int* newdata = convertComponents<int>(data, field.ncomponents, 3, ntuples);
            rv->SetArray((int*)newdata,ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
            
            delete data;
        }
        else
            rv->SetArray((int*)data,ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        
//        if(reverse_endian){
//            int *buff = (int *) rv->GetVoidPointer(0);
//            for (long long i = 0 ; i < ntotal ; i++)
//            {
//                int tmp;
//                int32_Reverse_Endian(buff[i], (unsigned char *) &tmp);
//                buff[i] = tmp;
//            }
//        }
      
        return rv;
    }
    else if(type == VisusSimpleIO::INT64){
        vtkLongArray *rv = vtkLongArray::New();
        rv->SetNumberOfComponents(ncomponents);
    
        // ?? is it correct to use long here ??
        if(isVector && dim < 3){
            
            long* newdata = convertComponents<long>(data, field.ncomponents, 3, ntuples);
            rv->SetArray((long*)newdata,ncomponents*ntuples,1,vtkDataArrayTemplate<long>::VTK_DATA_ARRAY_FREE);
            
            delete data;
        }
        else
            rv->SetArray((long*)data,ncomponents*ntuples,1,vtkDataArrayTemplate<long>::VTK_DATA_ARRAY_FREE);
        
//        if(reverse_endian){
//            long *buff = (long *) rv->GetVoidPointer(0);
//            for (long long i = 0 ; i < ntotal ; i++)
//            {
//                long tmp;
//                double64_Reverse_Endian(buff[i], (unsigned char *) &tmp);
//                buff[i] = tmp;
//            }
//        }
      
        return rv;
    }
    else if(type == VisusSimpleIO::FLOAT32){

        vtkFloatArray *rv = vtkFloatArray::New();
        rv->SetNumberOfComponents(ncomponents);
        
        if(isVector && dim < 3){
            float* newdata = convertComponents<float>(data, field.ncomponents, 3, ntuples);
            rv->SetArray((float*)newdata,ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
            
            delete data;
        }
        else
            rv->SetArray((float*)data,ncomponents*ntuples,1,vtkDataArrayTemplate<int>::VTK_DATA_ARRAY_FREE);
        
//        if(reverse_endian){
//            float *buff = (float *) rv->GetVoidPointer(0);
//            for (long long i = 0 ; i < ntotal ; i++)
//            {
//                float tmp;
//                float32_Reverse_Endian(buff[i], (unsigned char *) &tmp);
//                buff[i] = tmp;
//            }
//        }
    
        return rv;
    }
    else if(type == VisusSimpleIO::FLOAT64){

        vtkDoubleArray *rv = vtkDoubleArray::New();
        rv->SetNumberOfComponents(ncomponents);
        
        if(isVector && dim < 3){
            
            double* newdata = convertComponents<double>(data, field.ncomponents, 3, ntuples);
            rv->SetArray((double*)newdata,ncomponents*ntuples,1,vtkDataArrayTemplate<double>::VTK_DATA_ARRAY_FREE);
            
            delete data;
        }
        else{
          
          std::cout << "components " <<ncomponents << " tuples " <<ntuples;
            rv->SetArray((double*)data,ncomponents*ntuples,1,vtkDataArrayTemplate<double>::VTK_DATA_ARRAY_FREE);
        }
      
//        if(reverse_endian){
//            double *buff = (double *) rv->GetVoidPointer(0);
//            for (unsigned long long i = 0 ; i < ntotal ; i++)
//            {
//                double tmp;
//                double64_Reverse_Endian(buff[i], (unsigned char *) &tmp);
//                buff[i] = tmp;
//            }
//        }
    
        return rv;
    }
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
//  Programmer: spetruzza, bremer5
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

vtkDataArray *
avtIDXFileFormat::GetVar(int timestate, int domain, const char *varname)
{
    //std::cout<< rank << ": start getvar " << varname << " domain "<< domain;
    
    return queryToVtk(timestate, domain, varname);
    
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
//  Programmer: spetruzza, bremer5
//  Creation:   Mon Dec 10 15:06:44 PST 2012
//
// ****************************************************************************

vtkDataArray *
avtIDXFileFormat::GetVectorVar(int timestate, int domain, const char *varname)
{

  //std::cout<< rank << ": start getVectorVar " << varname << " domain "<< domain;
    
    return queryToVtk(timestate, domain, varname);
    
}
