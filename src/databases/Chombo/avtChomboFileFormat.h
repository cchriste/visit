/*****************************************************************************
*
* Copyright (c) 2000 - 2006, The Regents of the University of California
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
//                           avtChomboFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_Chombo_FILE_FORMAT_H
#define AVT_Chombo_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>

#include <vector>

#include <hdf5.h>

struct {
  int i;
  int j;
} typedef intvect2d;

struct {
  int i;
  int j;
  int k;
} typedef intvect3d;

struct {
  intvect2d lo;
  intvect2d hi;
} typedef box2d;

struct{
  intvect3d lo;
  intvect3d hi;
} typedef box3d;

union
{
  box2d b2;
  box3d b3;
} typedef box;


// ****************************************************************************
//  Class: avtChomboFileFormat
//
//  Purpose:
//      Reads in Chombo files as a plugin to VisIt.
//
//  Programmer: Hank Childs
//  Creation:   Thu Jan 19 11:17:14 PDT 2006
//
//  Modifications:
//
//    Hank Childs, Mon Jun 19 16:44:06 PDT 2006
//    Add support for ghost zones.
//
//    Brad Whitlock, Mon Sep 25 13:54:59 PST 2006
//    I added some fixes for getting cycle,time and for time-varying metadata.
//
// ****************************************************************************

class avtChomboFileFormat : public avtSTMDFileFormat
{
  public:
                       avtChomboFileFormat(const char *);
    virtual           ~avtChomboFileFormat();

    virtual const char    *GetType(void)   { return "Chombo"; };
    virtual void           FreeUpResources(void); 
    virtual void           ActivateTimestep(void);

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);
    virtual void          *GetAuxiliaryData(const char *var, int,
                                            const char *type, void *args,
                                            DestructorFunction &);
  
  protected:
    bool                   initializedReader;
    int                    dimension;
    hid_t                  file_handle;
    std::vector<std::string>  varnames;
    double                 dtime;
    int                    cycle;
    int                    max_level;
    int                    num_levels;
    std::vector<int>       numGhosts;
    std::vector<int>       patchesPerLevel;
    std::vector<int>       refinement_ratio;
    std::vector<double>    dx;

    std::vector<int>       lowI;
    std::vector<int>       hiI;
    std::vector<int>       lowJ;
    std::vector<int>       hiJ;
    std::vector<int>       lowK;
    std::vector<int>       hiK;

    void                   InitializeReader(void);
    void                   GetLevelAndLocalPatchNumber(int global_patch,
                                           int &level, int &local_patch) const;
    void                   CalculateDomainNesting(void);

    virtual int            GetCycle(void);
    virtual double         GetTime(void);
    virtual int            GetCycleFromFilename(const char *f) const;

    virtual bool           HasInvariantMetaData(void) const { return false; };
    virtual bool           HasInvariantSIL(void) const { return false; };
};


#endif


