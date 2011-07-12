/*****************************************************************************
*
* Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
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
//                            avtViSUSFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_ViSUS_FILE_FORMAT_H
#define AVT_ViSUS_FILE_FORMAT_H

#include <idx_input.h>

#include <database_exports.h>

#include <avtMTSDFileFormat.h>

#include <map>
#include <string>
#include <vector>


// ****************************************************************************
//  Class: avtViSUSFileFormat
//
//  Purpose:
//      Reads in ViSUS files as a plugin to VisIt.
//
//  Programmer: mcmiller -- generated by xml2avt
//  Creation:   Tue Sep 14 20:37:24 PST 2004
//
//  Modifications:
//
//    Mark C. Miller, Tue May 17 18:48:38 PDT 2005
//    Added timeState arg to PopulateDatabaseMetaData to satisfy new interface
//
//    Mark C. Miller, Tue Aug 16 13:56:55 PDT 2005
//    Added GetFile/GetTimeInfo
//
// ****************************************************************************

class avtViSUSFileFormat : public avtMTSDFileFormat
{
  public:
                       avtViSUSFileFormat(const char *);
    virtual           ~avtViSUSFileFormat();

    virtual const char    *GetType(void)   { return "ViSUS"; };
    virtual void           FreeUpResources(void); 

    virtual void           GetCycles(std::vector<int> &);
    virtual void           GetTimes(std::vector<double> &);
    virtual int            GetNTimesteps(void);

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);

    virtual bool           CanCacheVariable(const char *var);

    virtual void           RegisterDataSelections(
                               const std::vector<avtDataSelection_p> &selList,
                               std::vector<bool> *selectionsApplied);

  protected:

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *, int);

  private:

    void                   GetFile();
    void                   GetTimeInfo();
    void                   SetupRectilinearCoordinates();
    void                   SetupDomainAndZoneIndexing(int *zoneCounts,
                                                      int *stepSizes,
                                                      int *startZones,
                                                      int *domCounts,
                                                      int *domIndices);

    std::string fileName;
    IDX_file_descriptor idxFile;
    bool haveOpenedFile;

    bool useGetData3D;
    bool haveSetupCoordinates;
    bool ignoreDataSelections;

    int globalMin[3], globalMax[3];
    double origin3D[3], gridSpacing3D[3];

    int minTimeIndex, maxTimeIndex;
    std::vector<int> cycleVals;
    std::vector<double> timeVals;

    int dataDimension;

    int numFields;

    int globalZoneTotal;

    int globalZoneCount[3];

    double *globalNodeCoords[3];

    int *fieldSampleSize;

    std::map<std::string, int> fieldMap;

    int procNum;
    int procCount;

    std::vector<avtDataSelection_p> selList;
    std::vector<bool>              *selsApplied;

};


#endif
