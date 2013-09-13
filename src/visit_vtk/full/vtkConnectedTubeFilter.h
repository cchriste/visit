/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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

#ifndef __vtkConnectedTubeFilter_h
#define __vtkConnectedTubeFilter_h
#include <visit_vtk_exports.h>

#include "vtkPolyDataAlgorithm.h"

class vtkCellArray;
class vtkPoints;

// ****************************************************************************
//  Class:  vtkConnectedTubeFilter
//
//  Purpose:
//    A more suitable implementation of the vtkTubeFilter.  It assumes the
//    lines are all connected (it will return an error otherwise).
//
//  Programmer:  Jeremy Meredith
//  Creation:    November  1, 2002
//
//  Modifications:
//    Eric Brugger, Wed Jan  9 11:29:49 PST 2013
//    Modified to inherit from vtkPolyDataAlgorithm.
//
// ****************************************************************************
class VISIT_VTK_API vtkConnectedTubeFilter : public vtkPolyDataAlgorithm
{
  public:
    vtkTypeMacro(vtkConnectedTubeFilter,vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);
    bool BuildConnectivityArrays(vtkPolyData *);

    // Description:
    // Construct object with radius 0.5, radius variation turned off, the number 
    // of sides set to 3, and radius factor of 10.
    static vtkConnectedTubeFilter *New();

    // Description:
    // Set the minimum tube radius (minimum because the tube radius may vary).
    vtkSetClampMacro(Radius,float,0.0,VTK_FLOAT_MAX);
    vtkGetMacro(Radius,float);

    // Description:
    // Set the number of sides for the tube. At a minimum, number of sides is 3.
    vtkSetClampMacro(NumberOfSides,int,3,VTK_INT_MAX);
    vtkGetMacro(NumberOfSides,int);

    // Description:
    // Set a boolean to control whether to create normals.
    // DefaultNormalOn is set.
    vtkSetMacro(CreateNormals,int);
    vtkGetMacro(CreateNormals,int);
    vtkBooleanMacro(CreateNormals,int);

    // Description:
    // Turn on/off whether to cap the ends with polygons.
    vtkSetMacro(Capping,int);
    vtkGetMacro(Capping,int);
    vtkBooleanMacro(Capping,int);

  protected:
    // ************************************************************************
    //  Class:  PointSequence
    //
    //  Purpose:
    //    Encapsulates a single doubly connected point sequence.
    //
    //    Jean Favre, Tue May  7 16:38:37 CEST 2013
    //    Used vtkIdType where needed
    // ************************************************************************
    struct PointSequence
    {
        int length;
        vtkIdType *index;
        vtkIdType *cellindex;
      public:
        PointSequence();
        ~PointSequence();
        void Init(int maxlen);
        void Add(vtkIdType i, vtkIdType ci);
    };

    // ************************************************************************
    //  Class:  PointSequenceList
    //
    //  Purpose:
    //    Encapsulates a list of separate point sequences.
    //
    //  Modifications:
    //    Rich Cook and Hank Childs, Thu Oct  2 16:31:45 PDT 2008
    //    Added data member to support tubing over loops.
    //
    //    Jean Favre, Tue May  7 16:38:37 CEST 2013
    //    Used vtkIdType where needed
    // ************************************************************************
    class PointSequenceList
    {
      private:
        // connectivity data
        int          len;
        vtkIdType         *numneighbors;
        vtkIdType         *connectivity[2];
        vtkIdType         *cellindex;
        vtkPoints   *pts;

        // traversal variables
        bool  *visited;
        vtkIdType    index;
        bool   lookforloops;
      public:
        PointSequenceList();
        ~PointSequenceList();
        bool Build(vtkPoints *points, vtkCellArray *lines);
        void InitTraversal();
        bool GetNextSequence(PointSequence &seq);
    };

  protected:
    vtkConnectedTubeFilter();
    ~vtkConnectedTubeFilter();

    virtual int RequestData(vtkInformation *,
                            vtkInformationVector **,
                            vtkInformationVector *);

    float Radius;      // minimum radius of tube
    int NumberOfSides; // number of sides to create tube
    int CreateNormals; // true to create normals
    int Capping;       // true to cap ends

    PointSequenceList *pseqlist;

  private:
    vtkConnectedTubeFilter(const vtkConnectedTubeFilter&);  // Not implemented.
    void operator=(const vtkConnectedTubeFilter&);  // Not implemented.
};

#endif
