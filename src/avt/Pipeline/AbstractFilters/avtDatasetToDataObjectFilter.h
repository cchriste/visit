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
//                       avtDatasetToDataObjectFilter.h                      //
// ************************************************************************* //

#ifndef AVT_DATASET_TO_DATA_OBJECT_FILTER_H
#define AVT_DATASET_TO_DATA_OBJECT_FILTER_H

#include <pipeline_exports.h>

class     vtkObject;

#include <avtFilter.h>
#include <avtDatasetSink.h>


// ****************************************************************************
//  Class: avtDatasetToDataObjectFilter
//
//  Purpose:
//      A filter that takes in a dataset as input and has a data object as
//      output.
//
//  Programmer: Hank Childs
//  Creation:   May 31, 2001
//
//  Modifications:
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined constructor and destructor definitions to .C files
//    because certain compilers have problems with them.
//
//    Hank Childs, Fri Dec  3 14:28:02 PST 2004
//    Added variable name argument to SearchDataForDataExtents.
//
//    Jeremy Meredith, Thu Feb 15 11:44:28 EST 2007
//    Added support for rectilinear grids with an inherent transform.
//
//    Hank Childs, Sun Nov 28 06:19:25 PST 2010
//    Add methods for caching VTK objects in the database.
//
//    Hank Childs, Tue Nov 30 20:38:36 PST 2010
//    Add method SearchDataForSpatialExtents.
//
// ****************************************************************************

class PIPELINE_API avtDatasetToDataObjectFilter
    : virtual public avtFilter, virtual public avtDatasetSink
{
  public:
                       avtDatasetToDataObjectFilter();
    virtual           ~avtDatasetToDataObjectFilter();

  protected:
    //                 Note that these variables are only used when an active
    //                 variable is set.
    bool               activeVariableIsPointData;
    bool               hasPointVars;
    bool               hasCellVars;

    void               InputSetActiveVariable(const char *);
    virtual void       SearchDataForDataExtents(double *, const char *);
    virtual void       SearchDataForSpatialExtents(double *);
    virtual void       PreExecute(void);

    vtkObject         *FetchArbitraryVTKObject(int dependencies, const char *name, int dom, 
                                               int ts, const char *type);
    void               StoreArbitraryVTKObject(int dependencies, const char *name, int dom, 
                                               int ts, const char *type, vtkObject *);
};


#endif


