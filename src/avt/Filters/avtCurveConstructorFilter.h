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
//                        avtCurveConstructorFilter.h                        //
// ************************************************************************* //

#ifndef AVT_CURVE_CONSTRUCTOR_FILTER_H
#define AVT_CURVE_CONSTRUCTOR_FILTER_H
#include <filters_exports.h>

#include <avtDatasetToDatasetFilter.h>
#include <vectortypes.h>


// ****************************************************************************
//  Class: avtCurveConstructorFilter
//
//  Purpose:
//      A filter that will construct one uniform curve from fragments of
//      curves.
//
//  Programmer: Kathleen Bonnell
//  Creation:   Sat Apr 20 13:01:58 PST 2002
//
//  Modifications:
//    Kathleen Bonnell, Fri Jul 12 16:53:11 PDT 2002  
//    Removed vtk filters associated with label-creation.  Now handled by
//    the plot.
//
//    Kathleen Bonnell, Mon Dec 23 08:23:26 PST 2002
//    Added UpdateDataObjectInfo. 
//    
//    Hank Childs, Fri Oct  3 11:10:29 PDT 2003
//    Moved from /plots/Curve.  Renamed to CurveConstructorFilter.
//
//    Kathleen Bonnell, Tue Jun 20 16:02:38 PDT 2006 
//    Add PostExecute and outputArray.
//
//    Kathleen Bonnell, Thu Mar 19 17:42:14 PDT 2009
//    Added 'ForceConstruction', needed by curve queries.
//
//    Kathleen Bonnell, Mon Mar 23 09:53:17 PDT 2009
//    Removed 'ForceConstruction'.
//
//    Kathleen Bonnell, Thu Feb 17 09:18:43 PST 2011
//    Added CreateSingeOutput method.
//
// ****************************************************************************

class AVTFILTERS_API avtCurveConstructorFilter : public avtDatasetToDatasetFilter
{
  public:
                              avtCurveConstructorFilter();
    virtual                  ~avtCurveConstructorFilter();

    virtual const char       *GetType(void)  
                                       { return "avtCurveConstructorFilter"; };
    virtual const char       *GetDescription(void)
                                  { return "Constructing Curve"; };

  protected:
    doubleVector              outputArray;
    vtkDataSet               *CreateSingleOutput(avtDataTree_p inTree);

    virtual void              Execute(void);
    virtual void              PostExecute(void);
    virtual void              VerifyInput(void);
    avtContract_p
                           ModifyContract(avtContract_p spec);
    virtual void              UpdateDataObjectInfo(void);
};


#endif


