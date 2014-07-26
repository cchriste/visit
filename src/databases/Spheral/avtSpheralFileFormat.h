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
//                           avtSpheralFileFormat.h                          //
// ************************************************************************* //

#ifndef AVT_SPHERAL_FILE_FORMAT_H
#define AVT_SPHERAL_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>

#include <vector>
#include <string>
#include <visitstream.h>
#include <avtTypes.h>


class vtkDataArray;
class vtkPolyData;


// ****************************************************************************
//  Class: avtSpheralFileFormat
//
//  Purpose:
//      A file format reader for the Spheral++ code.
//
//  Programmer: Hank Childs
//  Creation:   March 12, 2003
//
//  Modifications:
//    Brad Whitlock, Tue Apr 15 11:59:39 PDT 2003
//    Made it work on Windows.
//
// ****************************************************************************

struct AllOfOneDomain
{
    std::vector< vtkPolyData * >                  meshes;
    std::vector<std::vector< vtkDataArray * > >   fields;
};


class avtSpheralFileFormat : public avtSTMDFileFormat
{
  public:
                          avtSpheralFileFormat(const char *);
    virtual              ~avtSpheralFileFormat();
    
    virtual const char   *GetType(void) { return "Spheral++ File Format"; };
    
    virtual vtkDataSet   *GetMesh(int, const char *);
    virtual vtkDataArray *GetVar(int, const char *);
    virtual vtkDataArray *GetVectorVar(int, const char *);

    virtual void          FreeUpResources(void);
    virtual void          PopulateDatabaseMetaData(avtDatabaseMetaData *);
    virtual bool          PerformsMaterialSelection(void) { return true; };
    virtual bool          HasVarsDefinedOnSubMeshes(void) { return true; };

    virtual void          RegisterVariableList(const char *,
                                          const std::vector<CharStrRef> &);

    virtual int           GetCycle(void) { return cycle; };
    virtual double        GetTime(void)  { return dtime; };

    virtual double        GetTimeFromFilename(const char *) const;
    virtual int           GetCycleFromFilename(const char *) const;


  protected:
    std::string           rootfile;
    int                   ndomains;
    int                   cycle;
    bool                  gotCycle;
    double                dtime;
    bool                  gotTime;
    bool                  readInMetaData;

    std::vector<std::string>                nodeLists;
    std::vector<int>                        nodeListSizes;
    std::vector<std::string>                positionField;
    std::vector<int>                        positionDimension;
    std::vector<std::string>                fields;
    std::vector< std::vector<bool> >        fieldDefinedOnNodeList;
    std::vector<avtVarType>                 fieldType;
    std::vector<int>                        fieldDim1;
    std::vector<int>                        fieldDim2;
    std::vector<std::string>                domain_files;
    std::vector<bool>                       read_domain;
    std::vector<bool>                       validNodeLists;

    std::vector<AllOfOneDomain>             cache;

    std::string                             current_file;

    void                  DetermineSubFiles(istream &, int);
    int                   GetLine(istream &, char *, std::vector<int> &);
    void                  ParseHeader(istream &);
    void                  ParseField(char *, int, std::vector<int> &,bool,int);
    void                  ParseNodeList(char *, int, std::vector<int> &);

    void                  ReadDomain(int);
    vtkPolyData          *ReadNodeList(istream &, int);
    vtkDataArray         *ReadField(istream &, int, int &);
    int                   GetNodeListIndexFromName(const char *);
    int                   GetFieldIndexFromName(const char *);
    void                  ReadInMetaData(void);
};


#endif


