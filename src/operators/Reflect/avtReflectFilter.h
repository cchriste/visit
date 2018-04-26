/*****************************************************************************
*
* Copyright (c) 2000 - 2018, Lawrence Livermore National Security, LLC
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
//  File: avtReflectFilter.h
// ************************************************************************* //

#ifndef AVT_Reflect_FILTER_H
#define AVT_Reflect_FILTER_H

#include <avtSIMODataTreeIterator.h>
#include <avtPluginFilter.h>

#include <ReflectAttributes.h>

class vtkDataArray;
class vtkDataSet;
class vtkPointSet;
class vtkRectilinearGrid;


// ****************************************************************************
//  Class: avtReflectFilter
//
//  Purpose:
//      A plugin operator for Reflect.
//
//  Programmer: childs -- generated by xml2info
//  Creation:   Thu Mar 7 10:35:24 PDT 2002
//
//  Modifications:
//
//    Hank Childs, Fri Sep  3 12:10:47 PDT 2010
//    Add Boolean to zero out velocities on the boundary.
//
//    Eric Brugger, Thu Jul 31 19:15:11 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
//    Alister Maguire, Tue Nov 15 11:51:26 PST 2016
//    Added ThreadSafe method.     
//
//    Alister Maguire, Tue Apr 10 13:22:11 PDT 2018
//    Added PlaneReflect for reflecting about an arbitrary plane, 
//    ghostsCreated to signify when ghost zones have been created, 
//    and changed names of previous reflect functions to AxisReflect*.
//
// ****************************************************************************

class avtReflectFilter : public virtual avtSIMODataTreeIterator,
                         public virtual avtPluginFilter
{
  public:
                         avtReflectFilter();
    virtual             ~avtReflectFilter();

    static avtFilter    *Create();

    virtual const char  *GetType(void)  { return "avtReflectFilter"; };
    virtual const char  *GetDescription(void)
                             { return "Reflecting the data"; };

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    ReflectAttributes     atts;
    double                xReflect;
    double                yReflect;
    double                zReflect;
    bool                  zeroOutVelocitiesOnBoundary;
    bool                  ghostsCreated;

    virtual bool          ThreadSafe(void) { return(true); };
    virtual void          PreExecute(void);
    virtual void          PostExecute(void);
    virtual avtDataTree_p ExecuteDataTree(avtDataRepresentation *);
    virtual avtContract_p
                          ModifyContract(avtContract_p);
    virtual void          UpdateDataObjectInfo(void);

    vtkDataSet           *PlaneReflect(vtkDataSet *, bool);
    vtkDataSet           *AxisReflect(vtkDataSet *, int);
    vtkDataSet           *AxisReflectRectilinear(vtkRectilinearGrid *, int);
    vtkDataSet           *AxisReflectPointSet(vtkPointSet *, int);
    vtkDataArray         *ReflectDataArray(vtkDataArray *, double);
    void                  HasNeighbor(int, bool &, bool &, bool &);
};


#endif
