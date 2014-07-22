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
//  File: avtLineoutFilter.h
// ************************************************************************* //

#ifndef AVT_LINEOUT_FILTER_H
#define AVT_LINEOUT_FILTER_H
#include <filters_exports.h>

#include <avtPluginDataTreeIterator.h>

class vtkDataSet;
class vtkIdList;
class vtkPoints;
class vtkRectilinearGrid;


// ****************************************************************************
//  Class: avtLineoutFilter
//
//  Purpose:
//      A plugin operator for Lineout.
//
//  Programmer: kbonnell -- generated by xml2info
//  Creation:   Thu Apr 25 16:01:28 PST 2002
//
//  Modifications:
//
//    Hank Childs, Tue Sep 10 16:46:57 PDT 2002
//    Re-worked memory management paradigm.
//
//    Kathleen Bonnell, Tue Dec 23 10:18:06 PST 2003 
//    Added ModifyContract. 
//
//    Kathleen Bonnell, Wed Jan 14 12:02:38 PST 2004 
//    Added PostExecute. 
//
//    Kathleen Bonnell, Thu Jul 29 09:55:49 PDT 2004 
//    Added Sampling, NoSampling, and CreatePolys methods.
//
//    Kathleen Bonnell, Wed Oct 20 17:35:10 PDT 2004 
//    Added arg to CreatePolys, added method CreatePolysFromOriginalCells, 
//    added var useOriginalCells. 
//
//    Kathleen Bonnell, Mon Jul 31 10:15:00 PDT 2006 
//    Curves represented as 1D RectilinearGrids instead of PolyData.
//
//    Kathleen Bonnell, Mon Aug 14 18:00:43 PDT 2006
//    Added ndims. 
//
//    Hank Childs, Thu Jan 24 09:44:45 PST 2008
//    Moved to /components/Filters and divorced from LineoutAttributes.
//
//    Hank Childs, Fri Jan 25 09:59:29 PST 2008
//    Remove ignoreGlobal, which was unused.
//
//    Kathleen Bonnell, Thu Mar  6 09:07:33 PST 2008 
//    Add AVT_FILTERS_API for build on windows.
//
//    Eric Brugger, Mon Jul 21 14:08:32 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class AVTFILTERS_API avtLineoutFilter : public avtDataTreeIterator
{
  public:
                             avtLineoutFilter();
    virtual                 ~avtLineoutFilter() {;};

    virtual const char      *GetType(void)  
                             { return "avtLineoutFilter"; };
    virtual const char      *GetDescription(void)
                             { return "Lineout"; };

    void                     SetPoint1(double *p)
                               { point1[0] = p[0]; point1[1] = p[1]; point1[2] = p[2]; };
    void                     SetPoint2(double *p)
                               { point2[0] = p[0]; point2[1] = p[1]; point2[2] = p[2]; };
    void                     SetSamplingOn(bool so)
                               { samplingOn = so; };
    void                     SetNumberOfSamplePoints(int nosp)
                               { numberOfSamplePoints = nosp; };

  protected:
    double                   point1[3];
    double                   point2[3];
    bool                     samplingOn;
    int                      numberOfSamplePoints;

    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    virtual void              PostExecute(void);
    virtual void              VerifyInput(void);
    virtual void              UpdateDataObjectInfo(void);
    virtual avtContract_p
                              ModifyContract(avtContract_p);

  private:
    bool                      useOriginalCells;
    vtkDataSet               *Sampling(vtkDataSet *, int);
    vtkDataSet               *NoSampling(vtkDataSet *, int);
    vtkRectilinearGrid       *CreateRGrid(vtkDataSet *, double *, double *,
                                          vtkPoints *, vtkIdList *);
    vtkRectilinearGrid       *CreateRGridFromOrigCells(vtkDataSet *, double *, 
                                          double *, vtkPoints *, vtkIdList *);

   int ndims;
};

#endif
