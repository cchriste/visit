// ************************************************************************* //
//                            avtFluentFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_Fluent_FILE_FORMAT_H
#define AVT_Fluent_FILE_FORMAT_H
#include <avtSTMDFileFormat.h>
#include <vector>
#include <string>
#include <map>

// Fluent plugin
#include <visitstream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <sstream>

#include "vtkPoints.h"
#include "vtkTriangle.h"
#include "vtkTetra.h"
#include "vtkQuad.h"
#include "vtkHexahedron.h"
#include "vtkPyramid.h"
#include "vtkWedge.h"
#include "vtkConvexPointSet.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDoubleArray.h"

using namespace std;


// ****************************************************************************
//  Class: avtFluentFileFormat
//
//  Purpose:
//      Reads in Fluent files as a plugin to VisIt.
//
//  Programmer: bdotson -- generated by xml2avt
//  Creation:   Fri Jun 30 15:02:33 PST 2006
//
//  Modifications:
//
//    Hank Childs, Fri Sep  8 14:30:04 PDT 2006
//    Added methods HasInvariantMetaData and HasInvariantSIL on the advice
//    of Terry Jordan.
//
// ****************************************************************************

class avtFluentFileFormat : public avtSTMDFileFormat
{
  public:
                       avtFluentFileFormat(const char *);
    virtual           ~avtFluentFileFormat() {;};

    virtual const char    *GetType(void)   { return "Fluent"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

    //Number of variables is dynamic
    virtual bool          HasInvariantMetaData(void) const 
                                { return false; }; 
    //Number of Domains is dynamic
    virtual bool          HasInvariantSIL(void) const
                                { return false; };

  protected:

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);
    int                    OpenCaseFile(const char *filename);
    int                    OpenDataFile(const char *filename);
    int                    GetCaseChunk ();
    void                   GetNumberOfCellZones();
    int                    GetCaseIndex();
    void                   LoadVariableNames();
    int                    GetDataIndex();
    int                    GetDataChunk();

    void                   ParseCaseFile();
    int                    GetDimension();
    void                   GetLittleEndianFlag();
    void                   GetNodesAscii();
    void                   GetNodesSinglePrecision();
    void                   GetNodesDoublePrecision();
    void                   GetCellsAscii();
    void                   GetCellsBinary();
    void                   GetFacesAscii();
    void                   GetFacesBinary();
    void                   GetPeriodicShadowFacesAscii();
    void                   GetPeriodicShadowFacesBinary();
    void                   GetCellTreeAscii();
    void                   GetCellTreeBinary();
    void                   GetFaceTreeAscii();
    void                   GetFaceTreeBinary();
    void                   GetInterfaceFaceParentsAscii();
    void                   GetInterfaceFaceParentsBinary();
    void                   GetNonconformalGridInterfaceFaceInformationAscii();
    void                   GetNonconformalGridInterfaceFaceInformationBinary();
    void                   CleanCells();
    void                   PopulateCellNodes();
    int                    GetCaseBufferInt(int ptr);
    float                  GetCaseBufferFloat(int ptr);
    double                 GetCaseBufferDouble(int ptr);
    void                   PopulateTriangleCell(int i);
    void                   PopulateTetraCell(int i);
    void                   PopulateQuadCell(int i);
    void                   PopulateHexahedronCell(int i);
    void                   PopulatePyramidCell(int i);
    void                   PopulateWedgeCell(int i);
    void                   PopulatePolyhedronCell(int i);
    void                   ParseDataFile();
    int                    GetDataBufferInt(int ptr);
    float                  GetDataBufferFloat(int ptr);
    double                 GetDataBufferDouble(int ptr);
    void                   GetData(int dataType);

//
//  Variables
//

ifstream FluentCaseFile;
ifstream FluentDataFile;
string CaseBuffer;
string DataBuffer;

vtkPoints           *Points;
vtkTriangle         *Triangle;
vtkTetra            *Tetra;
vtkQuad             *Quad;
vtkHexahedron       *Hexahedron;
vtkPyramid          *Pyramid;
vtkWedge            *Wedge;
vtkConvexPointSet   *ConvexPointSet;

struct Cell {
  int type;
  int zone;
  vector<int> faces;
  int parent;
  int child;
  vector<int> nodes;
};

struct Face {
  int type;
  int zone;
  vector<int> nodes;
  int c0;
  int c1;
  int periodicShadow;
  int parent;
  int child;
  int interfaceFaceParent;
  int interfaceFaceChild;
  int ncgParent;
  int ncgChild;
};

struct ScalarDataChunk {
  int subsectionId;
  int zoneId;
  vector<double> scalarData;
};

struct VectorDataChunk {
  int subsectionId;
  int zoneId;
  vector<double> iComponentData;
  vector<double> jComponentData;
  vector<double> kComponentData;
};


vector< Cell > Cells;
vector< Face > Faces;
map< int, string > VariableNames;
vector< int >  CellZones;
vector< ScalarDataChunk > ScalarDataChunks;
vector< VectorDataChunk > VectorDataChunks;

vector< vector<int> > SubSectionZones;
vector< int > SubSectionIds;
vector< int > SubSectionSize;

vector< string > ScalarVariableNames;
vector< int > ScalarSubSectionIds;
vector< string > VectorVariableNames;
vector< int > VectorSubSectionIds;

int LittleEndianFlag;
int GridDimension;
int DataPass;
int NumberOfScalars;
int NumberOfVectors;

};


#endif
