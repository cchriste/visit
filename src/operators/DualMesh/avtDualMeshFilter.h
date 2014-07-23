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
//  File: avtDualMeshFilter.h
// ************************************************************************* //

#ifndef AVT_DUALMESH_FILTER_H
#define AVT_DUALMESH_FILTER_H

#include <avtPluginDataTreeIterator.h>

#include <DualMeshAttributes.h>

class vtkDataArray;


// ****************************************************************************
//  Class: avtDualMeshFilter
//
//  Purpose:
//      Converts rectilinear mesh data between dual representations.
//
//        Nodes to Zones: Creates output zones centered at input nodes and 
//        converts point data to cell data.
//        
//        Zones to Nodes: Creates output nodes centered at input zone centers 
//        and converts cell data to point data.
//        
//        Auto: Looks at the primary varaible to determine conversion mode.
//        If there is no primary var (which is a valid case for a mesh plot)
//        will default to "Nodes to Zones". 
//
//  Programmer: Cyrus Harrison
//  Creation:   Wed May 7 15:59:34 PST 2008
//
//  Modifications:
//    Brad Whitlock, Wed Aug 15 12:21:41 PDT 2012
//    Override ModifyContract.
//
//    Eric Brugger, Wed Jul 23 12:06:25 PDT 2014
//    Modified the class to work with avtDataRepresentation.
//
// ****************************************************************************

class avtDualMeshFilter : public avtPluginDataTreeIterator
{
  public:
                         avtDualMeshFilter();
    virtual             ~avtDualMeshFilter();

    static avtFilter    *Create();

    virtual const char  *GetType(void)  { return "avtDualMeshFilter"; };
    virtual const char  *GetDescription(void)
                             { return "Dual Mesh"; };

    virtual void         SetAtts(const AttributeGroup*);
    virtual bool         Equivalent(const AttributeGroup*);

  protected:
    DualMeshAttributes   atts;
    std::string          actualVar;
    int                  actualMode;

    virtual void         UpdateDataObjectInfo(void);
    virtual avtDataRepresentation *ExecuteData(avtDataRepresentation *);
    virtual void         ExamineContract(avtContract_p contract);
    virtual avtContract_p ModifyContract(avtContract_p);

    vtkDataArray        *ExpandDual(vtkDataArray *coords);
    vtkDataArray        *ContractDual(vtkDataArray *coords);
    
};


#endif
