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

#ifndef AVT_OUTCAR_FILE_FORMAT_H
#define AVT_OUTCAR_FILE_FORMAT_H

#include <avtFileFormatInterface.h>
#include <avtMTSDFileFormat.h>
#include <visitstream.h>
#include <string>
#include <vector>

struct Atom
{
    int elementtype_index;
    float x;
    float y;
    float z; 
    float fx;
    float fy;
    float fz;
    float vx;
    float vy;
    float vz;
};

// ****************************************************************************
//  Class: avtOUTCARFileFormat
//
//  Purpose:
//      Reads in OUTCAR files as a plugin to VisIt.
//
//  Programmer: Jeremy Meredith
//  Creation:   August 29, 2006
//
//  Modifications:
//    Jeremy Meredith, Fri Feb 23 15:22:37 EST 2007
//    Added support for seeking directly to preset timesteps.
//
//    Jeremy Meredith, Fri Apr 20 15:01:36 EDT 2007
//    Added support for magnetization fields.
//
//    Jeremy Meredith, Tue Mar 10 17:42:20 EDT 2009
//    Added support for POTIM field to get time values.
//
//    Jeremy Meredith, Mon May 10 18:01:50 EDT 2010
//    Changed the way cycles and times are generated.
//
//    Jeremy Meredith, Thu Aug 12 16:26:24 EDT 2010
//    Allowed per-cycle changes in unit cell vectors.
//
//    Jeremy Meredith, Tue Oct 19 12:59:24 EDT 2010
//    Added support for optional velocities.
//
// ****************************************************************************

class avtOUTCARFileFormat : public avtMTSDFileFormat
{
  public:
    static bool        Identify(const std::string&);
    static avtFileFormatInterface *CreateInterface(
                       const char *const *list, int nList, int nBlock);

                       avtOUTCARFileFormat(const char *filename);
    virtual           ~avtOUTCARFileFormat() {;};

    virtual int            GetNTimesteps(void);

    virtual const char    *GetType(void)   { return "OUTCAR"; };
    virtual void           FreeUpResources(void); 
    virtual void           GetCycles(std::vector<int>&);
    virtual void           GetTimes(std::vector<double>&);

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

    virtual bool          HasInvariantMetaData(void) const { return false; };

  protected:
    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *,int);


    void OpenFileAtBeginning();
    void ReadAllMetaData();
    void ReadAtomsForTimestep(int);

    ifstream in;
    std::string filename;
    bool metadata_read;

    int ntimesteps;
    int natoms;

    struct UCV
    {
        double v[3][3];
        double *operator[](int i) { return v[i]; }
    };

    std::vector<istream::pos_type>   file_positions;

    bool has_velocities;
    bool has_magnetization;

    std::vector<float>               mags,magp,magd,magtot;

    std::vector<float>               free_energy;
    std::vector< std::vector<Atom> > allatoms;

    std::vector<UCV> unitCell;
    double potim;

    std::vector<std::string> element_names;
    std::vector<int>         element_types;
    std::vector<int>         element_counts;
};


#endif
