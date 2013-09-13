/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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
//                            avtXYZFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_XYZ_FILE_FORMAT_H
#define AVT_XYZ_FILE_FORMAT_H

#include <avtMTSDFileFormat.h>

#include <vector>

#define MAX_XYZ_VARS 6

// ****************************************************************************
//  Class: avtXYZFileFormat
//
//  Purpose:
//      Reads in XYZ molecular files.  Supports up to MAX_XYZ_VARS extra
//      data fields.  Supports multiple timesteps in one file.
//
//  Programmer: Jeremy Meredith
//  Creation:   June 14, 2007
//
//  Modifications:
//    Jeremy Meredith, Thu Jul 26 14:33:16 EDT 2007
//    Changed to read only the requested time steps on demand instead
//    of reading the whole file at once.  Keep track of file positions
//    for all time steps to accessing even late time steps is instantanous.
//
//    Jeremy Meredith, Wed May 20 11:17:19 EDT 2009
//    Added support for CrystalMakers more "human-friendly" (but
//    less computer friendly) style of XYZ file.
//    nAtoms is now stored per-timestep, not globally.
//
// ****************************************************************************

class avtXYZFileFormat : public avtMTSDFileFormat
{
  public:
                       avtXYZFileFormat(const char *);
    virtual           ~avtXYZFileFormat() {;};

    virtual int            GetNTimesteps(void);

    virtual const char    *GetType(void)   { return "XYZ"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

  protected:
    enum XYZStyle { Normal, CrystalMaker };

    ifstream                       in;
    std::vector<istream::pos_type> file_positions;
    std::string                    filename;
    bool                           metaDataRead;
    int                            nTimeSteps;
    int                            nVars;
    std::vector<int>               nAtoms;
    XYZStyle                       fileStyle;
    bool                           corruptElement;

    std::vector< std::vector<int> >   e;
    std::vector< std::vector<float> > x;
    std::vector< std::vector<float> > y;
    std::vector< std::vector<float> > z;
    std::vector< std::vector<float> > v[MAX_XYZ_VARS];

    virtual void    PopulateDatabaseMetaData(avtDatabaseMetaData *, int);
    void            OpenFileAtBeginning();

    void ReadTimeStep(int);
    void ReadAllMetaData();
};


#endif
