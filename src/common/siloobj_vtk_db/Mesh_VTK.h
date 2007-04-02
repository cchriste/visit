/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
* Produced at the Lawrence Livermore National Laboratory
* All rights reserved.
*
* This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
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
*    documentation and/or materials provided with the distribution.
*  - Neither the name of the UC/LLNL nor  the names of its contributors may be
*    used to  endorse or  promote products derived from  this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
* CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
* ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/

// ************************************************************************* //
//                                Mesh_VTK.h                                 //
// ************************************************************************* //

#ifndef MESH_VTK_H
#define MESH_VTK_H
#include <siloobj_vtk_exports.h>

#include <visitstream.h>

#include <silo.h>

#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>

#include <Field_VTK.h>
#include <IntervalTree_VTK.h>
#include <Mesh.h>


//
//  Forward declaration of classes.
//

class  TableOfContents;


// ****************************************************************************
//  Class: Mesh_VTK
//
//  Purpose:
//      A derived type of Mesh that constructs VTK objects.
//
//  Programmer: Hank Childs
//  Creation:   February 3, 2000
//
//  Modifications:
//    
//    Hank Childs, Sat Mar  4 10:36:28 PST 2000
//    Added GetCoords routine for rectilinear meshes.
//
//    Hank Childs, Mon Apr  3 15:06:55 PDT 2000
//    Added GetCoords routine for curvilinear meshes.
//
//    Hank Childs, Wed Apr 12 21:47:40 PDT 2000
//    Removed method GetDomainList, added method GetMetaData.
//
// ****************************************************************************

class SILOOBJ_VTK_API Mesh_VTK : public Mesh
{
  public:
                             Mesh_VTK();
    virtual                 ~Mesh_VTK();

    void                     GetCoords(const int *, int, 
                                       vtkUnstructuredGrid **);
    void                     GetCoords(const int *, int, 
                                       vtkRectilinearGrid **);
    void                     GetCoords(const int *, int, vtkStructuredGrid **);
    const IntervalTree_VTK  *GetMetaData(void);

    void                     UpdateReferences(TableOfContents *);

  protected:
    // The table of contents that holds this mesh.
    TableOfContents         *toc;

    // The coordinate field (as a pointed to the _VTK type).
    Field_VTK               *coordsField;
};


#endif


