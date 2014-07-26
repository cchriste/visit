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
//                        avtMTSDFileFormatInterface.h                       //
// ************************************************************************* //

#ifndef AVT_MTSD_FILE_FORMAT_INTERFACE_H
#define AVT_MTSD_FILE_FORMAT_INTERFACE_H

#include <database_exports.h>

#include <avtFileFormatInterface.h>
#include <avtMTSDFileFormat.h>
#include <vector>

class avtIOInformation;

// ****************************************************************************
//  Class: avtMTSDFileFormatInterface
//
//  Purpose:
//      This is an implementation of avtFileFormatInterface for file formats
//      that have Multiple Timesteps and a Single Domain.  It enables such
//      file formats to be treated by the database as if there were multiple
//      domains, but to be written as if there was only one domain.
//
//  Programmer: Hank Childs
//  Creation:   October 8, 2001
//
//  Modifications:
//    Kathleen Bonnell, Fri Feb  8 11:03:49 PST 2002
//    vtkScalars and vtkVectors havebeen deprecated in VTK 4.0, 
//    use vtkDataArray instead.
//
//    Hank Childs, Fri Mar 14 17:28:55 PST 2003
//    Removed calls that could be put in the base class.
//
//    Brad Whitlock, Wed May 14 09:24:22 PDT 2003
//    Added an optional timeState argument to SetDatabaseMetaData.
//
//    Mark C. Miller, Mon Feb 23 20:38:47 PST 2004
//    Added method, ActivateTimestep
//
//    Mark C. Miller, Tue May 17 18:48:38 PDT 2005
//    Added bool, forceReadAllCyclesTimes, to SetDatabaseMetaData
//
//    Mark C. Miller, Tue May 31 20:12:42 PDT 2005
//    Added method SetCycleTimeInDatabaseMetaData
//
//    Jeremy Meredith, Thu Jan 28 13:11:07 EST 2010
//    MTSD files can now be grouped not just into a faux MD format by having
//    more than one block, but also into a longer sequence of MT files,
//    each chunk with one or more timesteps.
//
//    Hank Childs, Tue Dec 20 14:15:33 CST 2011
//    Add method CreateCacheNameIncludingSelections.
//
//    Brad Whitlock, Thu Jun 19 10:50:25 PDT 2014
//    Pass mesh name to PopulateIOInformation.
//
// ****************************************************************************

class DATABASE_API avtMTSDFileFormatInterface : public avtFileFormatInterface
{
  public:
                    avtMTSDFileFormatInterface(avtMTSDFileFormat ***,
                                               int ntsgroups, int nblocks);
    virtual        ~avtMTSDFileFormatInterface();

    virtual vtkDataSet     *GetMesh(int, int, const char *);
    virtual vtkDataArray   *GetVar(int, int, const char *);
    virtual vtkDataArray   *GetVectorVar(int, int, const char *);
    virtual void           *GetAuxiliaryData(const char *var, int, int,
                                             const char *type, void *args,
                                             DestructorFunction &);
    virtual std::string     CreateCacheNameIncludingSelections(std::string,
                                                               int, int);

    virtual const char     *GetFilename(int);
    virtual void            SetDatabaseMetaData(avtDatabaseMetaData *md,
                                int timeState = 0,
                                bool forceReadAllCyclesTimes = false);
    virtual void            SetCycleTimeInDatabaseMetaData(
                                avtDatabaseMetaData *md, int ts);
    virtual void            FreeUpResources(int, int);

    virtual void            ActivateTimestep(int ts);

    virtual bool            PopulateIOInformation(int ts, const std::string &meshname,
                                                  avtIOInformation& ioInfo);

  protected:
    avtMTSDFileFormat    ***chunks;
    int                     nTimestepGroups;
    int                     nBlocks;
    std::vector<int>        tsPerGroup;
    int                     nTotalTimesteps;

    virtual int             GetNumberOfFileFormats(void)
                              { return nTimestepGroups*nBlocks; };
    virtual avtFileFormat  *GetFormat(int n) const
                              { int block = n % nBlocks;
                                int tsg   = n / nBlocks;
                                return chunks[tsg][block]; };
    void                    GenerateTimestepCounts();
    int                     GetTimestepGroupForTimestep(int ts);
    int                     GetTimestepWithinGroup(int ts);
};


#endif


