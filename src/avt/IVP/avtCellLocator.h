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

// ************************************************************************* //
//                              avtCellLocator.h                             //
// ************************************************************************* //

#ifndef AVT_CELLLOCATOR_H
#define AVT_CELLLOCATOR_H

#include <ivp_exports.h>
#include <avtVector.h>
#include <vtkType.h>
#include <vector>

#include <ref_ptr.h>

class vtkDataSet;

struct avtInterpolationWeight
{
    vtkIdType i;    // point id 
    double    w;    // point weight
};

typedef std::vector<avtInterpolationWeight> avtInterpolationWeights;

// ************************************************************************* //
//  Modifications:
//
//    Hank Childs, Fri Oct 29 12:13:07 PDT 2010
//    Add new data members for efficient curvilinear location.
//
//    Hank Childs, Fri Nov 19 14:45:53 PST 2010
//    Add support for voxels.
//
//    Hank Childs, Sun Nov 28 11:34:04 PST 2010
//    Add support for caching cell locators via void_ref_ptr.
//
//    David Camp, Tue Sep 13 08:16:35 PDT 2011
//    Changed the SetDataSet function to virtual. You may need to reset
//    pointer to the new dataset.
//
// ************************************************************************* //

class IVP_API avtCellLocator
{
  public:
                    avtCellLocator( vtkDataSet* );
    virtual        ~avtCellLocator();

    vtkDataSet     *GetDataSet() { return dataSet; }
    virtual void    SetDataSet(vtkDataSet *ds);
    void            ReleaseDataSet();

    virtual vtkIdType FindCell( const double pos[3], 
                                avtInterpolationWeights* iw,
                                bool ignoreGhostCells ) const = 0;
    static void     Destruct(void *);

  protected:

    void CopyCell( vtkIdType cellid, vtkIdType* ids, 
                   double pts[][3] ) const;

    bool TestCell( vtkIdType id, const double pos[3],
                   avtInterpolationWeights* iw,
                   bool ignoreGhostCells ) const;

    bool TestTet( vtkIdType id, const double pos[3], 
                  avtInterpolationWeights* iw ) const; 
    bool TestHex( vtkIdType id, const double pos[3], 
                  avtInterpolationWeights* iw ) const; 
    bool TestPrism( vtkIdType id, const double pos[3], 
                    avtInterpolationWeights* iw ) const;
    bool TestVoxel( vtkIdType id, const double pos[3], 
                    avtInterpolationWeights* iw ) const;

    vtkDataSet*    dataSet;
    vtkIdType*     cellIdxPtr;
    vtkIdType*     cellLocPtr;
    int*           strDimPtr;
    bool           normal2D;
    bool           normal3D;
    float*         fCoordPtr;
    double*        dCoordPtr;
    unsigned char* ghostPtr;
};

typedef ref_ptr<avtCellLocator> avtCellLocator_p;

#endif

