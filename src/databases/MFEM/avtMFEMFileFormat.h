/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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
//                            avtMFEMFileFormat.h                            //
// ************************************************************************* //

#ifndef AVT_MFEM_FILE_FORMAT_H
#define AVT_MFEM_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>
#include <avtDataSelection.h>
#include <vector>
#include <map>

#include "mfem.hpp"
class JSONRoot;


// ****************************************************************************
//  Class: avtMFEMFileFormat
//
//  Purpose:
//       MFEM database reader with LOD / refinement controls.
//
//  Programmer: harrison37 -- generated by xml2avt
//  Creation:   Fri May 23 15:16:20 PST 2014
//
//  Modifications:
//   Cyrus Harrison, Wed Sep 24 10:47:00 PDT 2014
//   Enable time varying metadata and SIL.
//
// ****************************************************************************

class avtMFEMFileFormat : public avtSTMDFileFormat
{
  public:
                       avtMFEMFileFormat(const char *);
    virtual           ~avtMFEMFileFormat();
    
    // VisIt can't cache for us b/c we need to implement LOD support. 
    virtual bool           CanCacheVariable(const char *var) {return false;}
    // Used to enable support for avtResolutionSelection
    virtual void        RegisterDataSelections(
                            const std::vector<avtDataSelection_p> &selList,
                            std::vector<bool> *selectionsApplied);


    virtual const char    *GetType(void)   { return "MFEM"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

    virtual bool           ReturnsValidCycle() const;
    virtual int            GetCycle();
    virtual bool           ReturnsValidTime() const;
    virtual double         GetTime();

    virtual void           ActivateTimestep(void);
  protected:

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);
    bool                   HasInvariantMetaData(void) const { return false; };
    bool                   HasInvariantSIL(void) const      { return false; };
  
  private:
    int                              selectedLOD;

    Mesh                            *FetchMesh(const std::string &mesh_name,
                                              int chunk);
                                              
    vtkDataSet                      *GetRefinedMesh(const std::string &mesh_name,
                                                    int chunk,
                                                    int lod);
                                                    
    vtkDataArray                    *GetRefinedVar(const std::string &mesh_name,
                                                   int chunk,
                                                   int lod);
                                                   
    vtkDataArray                    *GetRefinedElementColoring(const std::string &mesh_name, 
                                                                int domain, 
                                                                int lod);
                                                                
    vtkDataArray                    *GetRefinedElementAttribute(const std::string &mesh_name, 
                                                                int domain, 
                                                                int lod);
                                             
    JSONRoot                        *root;  
    
};


#endif
