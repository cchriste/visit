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
//                             avtMTSDFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_MTSD_FILE_FORMAT_H
#define AVT_MTSD_FILE_FORMAT_H

#include <database_exports.h>

#include <avtFileFormat.h>

#include <vector>

#include <void_ref_ptr.h>


class     vtkDataSet;
class     vtkDataArray;

class     avtIOInformation;


// ****************************************************************************
//  Class: avtMTSDFileFormat
//
//  Purpose:
//      This is an implementation of avtFileFormat for file formats that have
//      Multiple Timesteps and a Single Domain.
//
//  Programmer: Hank Childs
//  Creation:   February 23, 2001
//
//  Modifications:
//
//    Hank Childs, Mon Mar 11 08:48:41 PST 2002
//    Removed SetDatabaseMetaData since it is defined in the base class.
//
//    Kathleen Bonnell, Mon Mar 18 17:15:28 PST 2002  
//    vtkScalars and vtkVectors have been deprecated in VTK 4.0, 
//    use vtkDataArray instead.
//
//    Kathleen Bonnell, Mon Mar 18 17:15:28 PST 2002  
//    vtkScalars and vtkVectors have been deprecated in VTK 4.0, 
//    use vtkDataArray instead.
//
//    Brad Whitlock, Mon Oct 13 14:14:21 PST 2003
//    Added a method to get the times.
//
//    Mark C. Miller, Mon Feb 23 20:38:47 PST 2004
//    Added method, ActivateTimestep
//
//    Mark C. Miller, Tue Mar 16 14:28:42 PST 2004
//    Added method, PopulateIOInformation
//
//    Hank Childs, Mon Aug 16 16:22:56 PDT 2004
//    Allow for the domain to be set.
//
//    Mark C. Miller, Tue May 17 18:48:38 PDT 2005
//    Removed GetCycles/Times. Added SetDatabaseMetaData and both
//    time-qualified and non-time-qualified PopulateDatabaseMetaData methods
//    See note below.
//
//    Hank Childs, Sun May  9 18:47:06 CDT 2010
//    Add support for time slice offsets (used when group MT files with .visit
//    files).
//
//    Mark C. Miller, Fri Oct 29 09:58:43 PDT 2010
//    Moved implementation of SetDatabaseMetaData to the .C file.
//
//    Hank Childs, Tue Apr 10 15:12:58 PDT 2012
//    Add method SetReadAllCyclesAndTimes.
//
//    Brad Whitlock, Thu Jun 19 10:50:25 PDT 2014
//    Pass mesh name to PopulateIOInformation.
//
// ****************************************************************************

class DATABASE_API avtMTSDFileFormat : public avtFileFormat
{
  public:
                           avtMTSDFileFormat(const char * const *, int);
    virtual               ~avtMTSDFileFormat();

    virtual void          *GetAuxiliaryData(const char *var, int,
                                            const char *type, void *args,
                                            DestructorFunction &);

    void                   SetDomain(int d) { myDomain = d; };
    void                   SetTimeSliceOffset(int ts) { timeSliceOffset = ts; };

    virtual int            GetNTimesteps(void);

    virtual const char    *GetFilename(void) { return filenames[0]; };

    virtual vtkDataSet    *GetMesh(int, const char *) = 0;
    virtual vtkDataArray  *GetVar(int, const char *) = 0;
    virtual vtkDataArray  *GetVectorVar(int, const char *);

    virtual void           ActivateTimestep(int ts)
                               { avtFileFormat::ActivateTimestep(); };
    virtual bool           PopulateIOInformation(int ts, const std::string &meshname,
                                                 avtIOInformation& ioInfo) { return false; }
    virtual void           SetDatabaseMetaData(avtDatabaseMetaData *md, int ts = 0);

    void                   SetReadAllCyclesAndTimes(bool b) 
                                             { readAllCyclesAndTimes = b; };
    bool                   GetReadAllCyclesAndTimes(void)
                                             { return readAllCyclesAndTimes; };

  protected:
    char                 **filenames;
    int                    nFiles;
    int                    myDomain;
    int                    timeSliceOffset;
    bool                   readAllCyclesAndTimes;

    // The second of these should really be pure virtual and the first
    // non-existant. However, both are just virtual to maintain 
    // backward compatibility with older MTXX plugins and to allow 
    // MTXX plugins to implement a time-qualified request to populate
    // database metadata.
    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *md);
    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *md, int);

    int                    AddFile(const char *);
    static const int       MAX_FILES;
};


#endif


