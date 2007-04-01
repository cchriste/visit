
//========================================================================
//
//  Program:   VisIt
//  Module:    $RCSfile: vtkOnionPeelFilter.h,v $
//  Language:  C++
//  Date:      $Date: 2000/09/20 18:11:07 $
//  Version:   $Revision: 1.24 $
//
//
//
//=========================================================================
// .NAME vtkOnionPeelFilter - creates layers of cells around a given 
// seed cell 
// .SECTION Description
// vtkOnionPeelFilter is a filter object that, given an initial seed 
// cell in a rectilinear grid will create layers of "neighbor" cells around 
// the seed cell.  Neighbors are determined by either node adjacency (default) 
// or face adjacency.  Adjacency type can be controlled by the user.  
//
// To use this filter you should specifiy a starting seed cell and the number
// of layers. 
//
// Output is UnstructuredGrid
//
// .SECTION See Also
//
// .SECTION Caveats
// First element in layerCellIds is always SeedCellId.
// Layer zero offset is always zero, so first element in layerOffsets
// will be the AdjacencyType (in case of updates to adjacency
// while OnionPeelFilter is active)
//
//

#ifndef __vtkOnionPeelFilter_h
#define __vtkOnionPeelFilter_h
#include <visit_vtk_exports.h>

#include <vtkDataSetToUnstructuredGridFilter.h>
#include <vtkSetGet.h>

#define VTK_NODE_ADJACENCY 0
#define VTK_FACE_ADJACENCY 1

class vtkIdList;



//*****************************************************************************
//  Modifications:
//    Kathleen Bonnell, Thu Aug 15 18:37:59 PDT 2002  
//    Added a bool return for Intialize method.  Changed SetSeedCell method 
//    from a Macro to a regular method. Added logicalIndex and its associated
//    Set/Get methods, added useLogicalIndex.
//
//    Kathleen Bonnell, Tue Jan 18 19:37:46 PST 2005 
//    Added  bool 'ReconstructOriginalCells' and Set/Get methods.  Added
//    FindCellsCorrespondingToOriginal.
// 
//    Kathleen Bonnell, Wed Jan 19 15:54:38 PST 2005 
//    Took 'Cell' out of callback name, renamed 'SeedCellId' to 'SeedId'. 
//    Added 'SeedIdIsForCell'.  Added 'FindNodesCorrespondingToOriginal'.
//
//*****************************************************************************



typedef void (*BadSeedCallback)(void *, int, int, bool);

class VISIT_VTK_API  
vtkOnionPeelFilter : public vtkDataSetToUnstructuredGridFilter
{
public:
  static vtkOnionPeelFilter *New();
  vtkTypeMacro(vtkOnionPeelFilter,vtkDataSetToUnstructuredGridFilter);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the current Seed value.
  void SetSeedId(const int);
  vtkGetMacro(SeedId, int);

  // Description:
  // Turn on/off scaling of source geometry.
  vtkSetMacro(SeedIdIsForCell,int);
  vtkBooleanMacro(SeedIdIsForCell,int);
  vtkGetMacro(SeedIdIsForCell,int);

  // Description:
  // Set the current LogicalIndex value.
  void SetLogicalIndex(const int, const int, const int k = 0);
  int *GetLogicalIndex(void) { return logicalIndex; };

  // Description:
  // Set the current layer value.
  vtkSetClampMacro(RequestedLayer,int, 0, VTK_LARGE_INTEGER);
  vtkGetMacro(RequestedLayer,int);

  // Description:
  // Turn on/off scaling of source geometry.
  vtkSetMacro(ReconstructOriginalCells,int);
  vtkBooleanMacro(ReconstructOriginalCells,int);
  vtkGetMacro(ReconstructOriginalCells,int);

 

  // Description:
  // Specify which type of adjacency to use when determining neighbor cells.
  // There are two choices:  Face Adjacency and Node Adjacency.
  vtkSetClampMacro(AdjacencyType, int, VTK_NODE_ADJACENCY, VTK_FACE_ADJACENCY);
  vtkGetMacro(AdjacencyType, int);
  void SetAdjacencyTypeToFace()
       { this->SetAdjacencyType(VTK_FACE_ADJACENCY); };
  void SetAdjacencyTypeToNode()
       { this->SetAdjacencyType(VTK_NODE_ADJACENCY); };
  char *GetAdjacencyTypeAsString();

  bool Initialize(const int = VTK_LARGE_INTEGER);

  void SetBadSeedCallback(BadSeedCallback, void *);
 
protected:
// Protected Methods

  vtkOnionPeelFilter();
  ~vtkOnionPeelFilter();

  void Execute();
  void Grow();
  void GenerateOutputGrid();

  void FindCellNeighborsByNodeAdjacency(vtkIdList *, vtkIdList*);
  void FindCellNeighborsByFaceAdjacency(vtkIdList *, vtkIdList*);
  void FindCellsCorrespondingToOriginal(int, vtkIdList*);
  void FindCellsCorrespondingToOriginal(vtkIdList *, vtkIdList*);
  void FindNodesCorrespondingToOriginal(int, vtkIdList*);

// Protected Data Members

  vtkIdList *layerCellIds;
  vtkIdList *cellOffsets;

  int maxLayersReached;
  int maxLayerNum;
  int RequestedLayer;
  int AdjacencyType;
  int SeedId;
  int ReconstructOriginalCells;
  int SeedIdIsForCell;

  int logicalIndex[3];
  bool useLogicalIndex;
  BadSeedCallback  bsc_callback;
  void                *bsc_args;
  
private:
  vtkOnionPeelFilter(const vtkOnionPeelFilter&);
  void operator=(const vtkOnionPeelFilter&);

};

inline char *vtkOnionPeelFilter::GetAdjacencyTypeAsString(void)
{
  if ( this->AdjacencyType == VTK_FACE_ADJACENCY )
    {
    return "FaceAdjacency";
    }
  else
    {
    return "NodeAdjacency";
    }
}

#endif


