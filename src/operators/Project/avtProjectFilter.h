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
//  File: avtProjectFilter.h
// ************************************************************************* //

#ifndef AVT_Project_FILTER_H
#define AVT_Project_FILTER_H

#include <avtPluginDataTreeIterator.h>

#include <ProjectAttributes.h>

class vtkDataSet;
class vtkPointSet;
class vtkRectilinearGrid;


// ****************************************************************************
//  Class: avtProjectFilter
//
//  Purpose:
//      A plugin operator for Project.
//
//  Programmer: Jeremy Meredith
//  Creation:   Tue May 18 14:35:38 PST 2004
//
//  Modifications:
//    Jeremy Meredith, Thu Sep  9 16:44:38 PDT 2004
//    Added ModifyContract so pick could get the node/zone numbers if
//    needed.  Force a rectilinear dataset to always project into a
//    curvilinear one.  Added projection of vectors.
//
//    Hank Childs, Thu Jan 20 10:27:29 PST 2005
//    Added extents calculation in PostExecute.
//
//    Jeremy Meredith, Thu Apr  1 14:41:33 EDT 2010
//    Added double precision version of projection for increased
//    accuracy.
//
//    Eric Brugger, Thu Jul 31 14:42:52 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class avtProjectFilter : public avtPluginDataTreeIterator
{
  public:
                        avtProjectFilter();
    virtual            ~avtProjectFilter();

    static avtFilter   *Create();

    virtual const char *GetType(void)  { return "avtProjectFilter"; }
    virtual const char *GetDescription(void) { return "Project"; }

    virtual void        SetAtts(const AttributeGroup*);
    virtual bool        Equivalent(const AttributeGroup*);

  protected:
    ProjectAttributes   atts;

    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    virtual void                PostExecute(void);
    virtual void                UpdateDataObjectInfo(void);
    avtContract_p  ModifyContract(avtContract_p);

  private:
    void                ProjectPoint(float &x, float &y, float &z);
    void                ProjectPoint(double &x, double &y, double &z);
    vtkPointSet        *ProjectPointSet(vtkPointSet*);
    vtkPointSet        *ProjectRectilinearGrid(vtkRectilinearGrid*);

    template <class T> 
    void                ProjectVectors(vtkDataSet*,vtkDataSet*,
                                       vtkDataArray*,vtkDataArray*,bool);
};


#endif
