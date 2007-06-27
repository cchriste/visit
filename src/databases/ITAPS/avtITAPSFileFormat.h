/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
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
//                            avtITAPSFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_ITAPS_FILE_FORMAT_H
#define AVT_ITAPS_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>

#include <iBase.h>
#include <iMesh.h>

#include <map>
#include <string>
#include <vector>

using std::map;
using std::string;
using std::vector;

typedef struct _VisItEntitySetInfo {
    int domainId;
    unsigned int groupId;
    int matId;
    vector<string> varNames;
} VisItEntitySetInfo_t;

// ****************************************************************************
//  Class: avtITAPSFileFormat
//
//  Purpose:
//      Reads in ITAPS files as a plugin to VisIt.
//
//  Programmer: Mark C. Miller 
//  Creation:   Wed Mar 7 17:15:33 PST 2007
//
//  Modifications:
//
//    Mark C. Miller, Thu Mar 22 09:37:55 PDT 2007
//    Added handles for primitiveEntities so we don't attempt to get them
//    more than once.
//    Added stuff to start processing the entity set hierarchy. It is not
//    yet in use though.
//
// ****************************************************************************

class avtITAPSFileFormat : public avtSTMDFileFormat
{
  public:
                       avtITAPSFileFormat(const char *);
    virtual           ~avtITAPSFileFormat();

    virtual const char    *GetType(void)   { return "ITAPS"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

  protected:

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);

  private:
      string               vmeshFileName;
      iMesh_Instance       itapsMesh;
      iBase_EntitySetHandle rootSet;
      iBase_EntityHandle  *vertEnts;
      int                  vertEnts_allocated;
      iBase_EntityHandle  *edgeEnts;
      int                  edgeEnts_allocated;
      iBase_EntityHandle  *faceEnts;
      int                  faceEnts_allocated;
      iBase_EntityHandle  *regnEnts;
      int                  regnEnts_allocated;
      int                  geomDim;
      int                  topoDim;
      int                  numVerts;
      int                  numEdges;
      int                  numFaces;
      int                  numRegns;
      iBase_EntityType     domToEntType[4];
      bool                 haveMixedElementMesh;
      vector<iBase_TagHandle>    primitiveTagHandles[4];
      map<iBase_EntitySetHandle,VisItEntitySetInfo_t> esMap;

};


#endif
