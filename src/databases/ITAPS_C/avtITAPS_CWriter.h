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
//                           avtITAPS_CWriter.h                              //
// ************************************************************************* //

#ifndef AVT_ITAPS_C_WRITER_H
#define AVT_ITAPS_C_WRITER_H

#include <avtDatabaseWriter.h>

#include <iMesh.h>

#include <string>
#include <vector>

class avtMeshMetaData;
class DBOptionsAttributes;
class vtkCellData;

// ****************************************************************************
//  Class: avtITAPS_CWriter
//
//  Purpose:
//      A module that writes out ITAPS_C files.
//
//  Programmer: Mark C. Miller 
//  Creation:   November 20, 2008 
//
//  Modifications:
//    Mark C. Miller, Wed Jan 14 17:54:21 PST 2009
//    Added some bools to control behavior of output to iMesh 
//
//    Mark C. Miller, Tue Apr 21 15:53:42 PDT 2009
//    Added storage for spatial and topoligical dimension of mesh.
//
//    Mark C. Miller, Wed May 20 09:56:36 PDT 2009
//    Added WriteMaterial method
// ****************************************************************************

class avtITAPS_CWriter : public virtual avtDatabaseWriter
{
  public:
                   avtITAPS_CWriter(DBOptionsAttributes *);
    virtual       ~avtITAPS_CWriter();

  protected:

    virtual bool   CanHandleMaterials(void) { return true; };

    virtual void   OpenFile(const std::string &, int);
    virtual void   WriteHeaders(const avtDatabaseMetaData *,
                                std::vector<std::string> &, 
                                std::vector<std::string> &,
                                std::vector<std::string> &);
    virtual void   WriteChunk(vtkDataSet *, int);
    virtual void   CloseFile(void);

    virtual void   WriteMaterial(vtkCellData *cd, int chunk,
                       iMesh_Instance itapsMesh,
                       iBase_EntitySetHandle rootSet,
                       iBase_EntitySetHandle chunkSet,
                       iBase_EntityHandle *clHdls);

    enum iMesh_EntityTopology VTKZoneTypeToITAPS_MOABZoneType(int);

  private:
    std::string    stem;
    std::string    dir;
    std::string    saveOptions;
    std::string    formatExtension;
    bool           simplexify;
    bool           addFacesFor3DEnts;
    bool           preventDupsToiMesh;
    int            nblocks;
    int            spatialDim;
    int            topoDim;
};

#endif
