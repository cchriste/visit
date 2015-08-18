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
//                           avtPLOT3DFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_PLOT3D_FILE_FORMAT_H
#define AVT_PLOT3D_FILE_FORMAT_H

#include <avtMTMDFileFormat.h>

#include <vector>
#include <string>


class     DBOptionsAttributes;
class     vtkPLOT3DReader;


// ****************************************************************************
//  Class: avtPLOT3DFileFormat
//
//  Purpose:
//      A file format reader for PLOT3D files.
//
//  Programmer: Hank Childs
//  Creation:   May 3, 2002
//
//  Modifications:
//    Kathleen Biagas, Thu Apr 23 10:36:09 PDT 2015
//    Added 'haveSolutionFile' flag.
//
//    Kathleen Biagas, Fri Jun 26 10:24:26 PDT 2015
//    Change this from type STMD to MTMD.
//    Add solutionFiles, times, haveReadMetaFile, haveProcessedQ, previousTS.
//
// ****************************************************************************

class avtPLOT3DFileFormat : public avtMTMDFileFormat
{
  public:
                          avtPLOT3DFileFormat(const char *, DBOptionsAttributes *);
    virtual              ~avtPLOT3DFileFormat();
    
    virtual void           GetTimes(std::vector<double> &);
    virtual int            GetNTimesteps(void);

    virtual const char   *GetType(void) { return "PLOT3D File Format"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, int, const char *);
    virtual vtkDataArray  *GetVar(int, int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, int, const char *);

    virtual void           ActivateTimestep(int ts);
  protected:
    vtkPLOT3DReader *reader;
    std::string           visitMetaFile;
    std::string           xFileName;
    std::string           qFileName;
    std::string           solutionRoot;
    std::vector<std::string> solutionFiles;
    std::vector<double>   times;
    bool                  haveSolutionFile;
    bool                  haveReadMetaFile;
    bool                  haveProcessedQ;
    int                   previousTS;

    virtual void          PopulateDatabaseMetaData(avtDatabaseMetaData *, int);

  private:
    bool                  ReadVisItMetaFile(void);
    bool                  ProcessQForTimeSeries(void);
    void                  SetTimeStep(int timeState);
};


#endif


