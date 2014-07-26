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
//                            avtVLIFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_VLI_FILE_FORMAT_H
#define AVT_VLI_FILE_FORMAT_H

#include <vtkObject.h>
#include <avtDataSelection.h>
#include <avtMTMDFileFormat.h>
#include <VLIFileManager.h>

#include <vector>


// ****************************************************************************
//  Class: avtVLIFileFormat
//
//  Purpose:
//      Reads in VLI files as a plugin to VisIt.
//
//  Programmer: Markus Glatter <glatter@cs.utk.edu> -- generated by xml2avt
//  Creation:   Mon May 7 13:54:06 PST 2007
//
// ****************************************************************************

class avtVLIFileFormat : public avtMTMDFileFormat
{
  public:
                           avtVLIFileFormat(const char *);
    virtual               ~avtVLIFileFormat();

    virtual const char    *GetType(void)   { return "VLI"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, int, const char *);
    virtual int            GetNTimesteps(void);
    virtual vtkDataArray  *GetVar(int, int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, int, const char *);
    virtual void           RegisterDataSelections(const std::vector<avtDataSelection_p>&, std::vector<bool> *);
    virtual void           RegisterVariableList(const char *, const std::vector<CharStrRef>&);
    virtual void          *threadCommServer(void *in);
    virtual void          *threadPacer(void *); 

   protected:

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *, int);
    virtual void           startServers(void);
    virtual int            startCommServer(std::string);
    virtual float          ConvertToFloat(unsigned short, int);
    virtual vtkObject     *Query(int, int, const char *);
    virtual void           ClearCache();
    virtual bool           CanCacheVariable(const char *v);

  private:

    int                             initialized;  
    int                             procNum;
    int                             procCount;
    VLIFileManager                 *config;
    volatile bool                   error;
    volatile bool                   loaded;
    volatile bool                   info;
    std::string                    *shostname;
    int                            *sport;
    std::string                    *ehostname;
    int                            *eport;
    void                           *ptid;
    volatile int                    socket;
    std::vector<avtDataSelection_p> selList;
    std::vector<bool>              *selsApplied;
    std::vector<int>                registeredVars;
    int                             queryTimestate;
    int                             queryDomain;
    std::vector<vtkObject *>        queryObjects;
    void                           *pacertid;
    int                             pport;
    volatile bool                   inQuery;
};

#endif
