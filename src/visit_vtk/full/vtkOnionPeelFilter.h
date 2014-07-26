/*****************************************************************************
*
* Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
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

#include <vtkUnstructuredGridAlgorithm.h>

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
//    Jeremy Meredith, Thu Aug  7 14:21:23 EDT 2008
//    Adjacency type string should have been const char*.
//
//    Kathleen Biagas, Wed Jul 23 16:43:42 MST 2014
//    Added 'maxId' arg to FindCellsCorrespondingToOriginal, used in
//    error reporting.
//
//*****************************************************************************

typedef void (*BadSeedCallback)(void *, int, int, bool);

class VISIT_VTK_API  
vtkOnionPeelFilter : public vtkUnstructuredGridAlgorithm
{
public:
  vtkTypeMacro(vtkOnionPeelFilter,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  static vtkOnionPeelFilter *New();

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
  vtkSetClampMacro(RequestedLayer,int, 0, VTK_INT_MAX);
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
  const char *GetAdjacencyTypeAsString();

  bool Initialize(vtkDataSet *);

  void SetBadSeedCallback(BadSeedCallback, void *);
 
protected:
// Protected Methods

  vtkOnionPeelFilter();
  ~vtkOnionPeelFilter();

  virtual int RequestData(vtkInformation *,
                          vtkInformationVector **,
                          vtkInformationVector *);
  virtual int FillInputPortInformation(int port, vtkInformation *info);

  void Grow(vtkDataSet *);
  void GenerateOutputGrid(vtkDataSet *, vtkUnstructuredGrid *);

  void FindCellNeighborsByNodeAdjacency(vtkDataSet *, vtkIdList *, vtkIdList*);
  void FindCellNeighborsByFaceAdjacency(vtkDataSet *, vtkIdList *, vtkIdList*);
  void FindCellsCorrespondingToOriginal(vtkDataSet *, int, vtkIdList*, int &);
  void FindCellsCorrespondingToOriginal(vtkDataSet *, vtkIdList *, vtkIdList*, int &);
  void FindNodesCorrespondingToOriginal(vtkDataSet *, int, vtkIdList*);

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

inline const char *vtkOnionPeelFilter::GetAdjacencyTypeAsString(void)
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
