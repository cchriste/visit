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

// ************************************************************************* //
//                             avtScatterFilter.h                            //
// ************************************************************************* //

#ifndef AVT_SCATTER_FILTER_H
#define AVT_SCATTER_FILTER_H

#include <avtDataTreeIterator.h>

#include <ScatterAttributes.h>


// ****************************************************************************
//  Class: avtScatterFilter
//
//  Purpose:
//      A filter that combines multiple scalar fields into a point mesh.
//
//  Programmer: Brad Whitlock
//  Creation:   Tue Nov 2 22:31:23 PST 2004
//
//  Modifications:
//    Brad Whitlock, Mon Jul 18 11:53:27 PDT 2005
//    Added overrides of PreExecute, ModifyContract, and added new 
//    method PopulateDataInputs. Also added extents members.
//
//    Cyrus Harrison, Tue Aug 17 11:52:24 PDT 2010
//    Added PostExecute method to:
//      1) Set legend if a color var is selected.
//      (This logic was moved out of ExecuteData to prevent a parallel hang
//       when there are more procs than chunks to process.)
//      2) Set proper spatial extents.
//
//    Cyrus Harrison, Thu Aug 19 13:35:08 PDT 2010
//    Changes to support using var1 from atts.
//
//    Hank Childs, Thu Aug 26 13:47:30 PDT 2010
//    Change extents name.
//
//    Kathleen Biagas, Thu Mar  1 14:49:50 MST 2012
//    Add keepNodeZone and dataArray (origNodes) arg to PointMeshFromVariables.
//
//    Eric Brugger, Tue Aug 19 11:13:13 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class avtScatterFilter : public avtDataTreeIterator
{
  public:
                               avtScatterFilter(const ScatterAttributes &);
    virtual                   ~avtScatterFilter();

    virtual const char        *GetType(void)  { return "avtScatterFilter"; };
    virtual const char        *GetDescription(void)
                                   { return "Creating point mesh"; };

protected:
    struct DataInput
    {
        vtkDataArray *data;
        bool          useMin;
        bool          useMax;
        float         min;
        float         max;
        int           scale;
        float         skew;
    };

    ScatterAttributes          atts;
    double                     xExtents[2];
    double                     yExtents[2];
    double                     zExtents[2];
    double                     colorExtents[2];
    bool                       needXExtents;
    bool                       needYExtents;
    bool                       needZExtents;
    bool                       needColorExtents;
    bool                       keepNodeZone;

    doubleVector               thisProcsSpatialExtents;

    virtual void               PreExecute(void);
    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    virtual void               PostExecute(void);
    virtual void               UpdateDataObjectInfo(void);
    virtual avtContract_p     
                               ModifyContract(avtContract_p spec);
    vtkDataArray              *GetDataArray(vtkDataSet *inDS,
                                            const std::string &name,
                                            avtCentering targetCentering,
                                            bool &deleteArray);
    vtkDataArray              *Recenter(vtkDataSet *ds, vtkDataArray *arr, 
                                        avtCentering cent) const;
    vtkDataSet                *PointMeshFromVariables(DataInput *d1,
                                                      DataInput *d2,
                                                      DataInput *d3,
                                                      DataInput *d4, bool &,
                                                      vtkDataArray *);
    int                        CountSpatialDimensions() const;
    void                       PopulateDataInputs(DataInput *orderedArrays,
                                                  vtkDataArray **arr) const;
    void                       PopulateNames(const char **) const;
    bool                       NeedSpatialExtents() const;
};


#endif


