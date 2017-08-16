/*****************************************************************************
*
* Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
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
//                            avtUintahFileFormat.h                          //
// ************************************************************************* //

#ifndef AVT_UINTAH_FILE_FORMAT_H
#define AVT_UINTAH_FILE_FORMAT_H

#include <StandAlone/tools/uda2vis/udaData.h>

#include <avtMTMDFileFormat.h>

#include <vector>
#include <map>
#include <string>

// ****************************************************************************
//  Class: avtUintahFileFormat
//
//  Purpose:
//      Reads in Uintah files as a plugin to VisIt.
//
//  Programmer: sshankar -- generated by xml2avt
//  Creation:   Tue May 13 19:02:26 PST 2008
//
// ****************************************************************************

// the DataArchive and GridP types are considered opaque from this library, we dlsym()
// functions from uintah libs to allocate, free, and perform operations on them.
class DataArchive;
class GridP;
class DBOptionsAttributes;

class avtUintahFileFormat : public avtMTMDFileFormat
{
public:
  avtUintahFileFormat( const char * filename, DBOptionsAttributes* attrs);
  virtual           ~avtUintahFileFormat();

  virtual double        GetTime( int timestep );
  virtual int           GetNTimesteps( void );

  virtual const char    *GetType( void )   { return "Uintah"; };
  virtual void          ActivateTimestep( int timestep ); 

  virtual vtkDataSet    *GetMesh( int timestate, int domain, const char * meshname );
  virtual vtkDataArray  *GetVar(  int timestate, int domain, const char * varname );
  virtual vtkDataArray  *GetVectorVar( int timestate, int domain, const char * varname );


protected:

  virtual void     PopulateDatabaseMetaData(avtDatabaseMetaData *, int);
  void             ReadMetaData(avtDatabaseMetaData *, int);

  virtual void     *GetAuxiliaryData(const char *var, int,
                                     const char *type, void *args,
                                     DestructorFunction &);

  void             GetLevelAndLocalPatchNumber(int, int&, int&);
  int              GetGlobalDomainNumber(int, int);
  void             CalculateDomainNesting(int, const std::string&);
        
  virtual bool     HasInvariantMetaData(void) const { return !dataVariesOverTime; };
  virtual bool     HasInvariantSIL(void) const { return !dataVariesOverTime; };

  void             AddExpressionsToMetadata(avtDatabaseMetaData *md);
  void             CheckNaNs(double *data, const int num,
                             const char* varname,
                             const int level, const int patch);

  // DATA MEMBERS
  bool useExtraCells;
  bool dataVariesOverTime;
  int currTimeStep;
  bool forceMeshReload;

  std::string mesh_for_patch_data;
  
  // VisIt meshes (see https://visitbugs.ornl.gov/issues/52)
  std::map<std::string, void_ref_ptr> mesh_domains;
  std::map<std::string, void_ref_ptr> mesh_boundaries;

  // Data that is not dependent on time
  DataArchive *archive;
  std::vector<double> cycleTimes;

  // Data that is dependent on time
  GridP *grid;
  TimeStepInfo *stepInfo;
        
  // Interface to the uda2vis library
  void  * libHandle;

  DataArchive*     (*openDataArchive)(const std::string&);
  void             (*closeDataArchive)(DataArchive*);

  GridP*           (*getGrid)(DataArchive*, int);
  void             (*releaseGrid)(GridP*);

  std::vector<double>   (*getCycleTimes)(DataArchive*);
  TimeStepInfo*    (*getTimeStepInfo)(DataArchive*, GridP*, int, bool);

  GridDataRaw*     (*getGridData)(DataArchive*, GridP*, int, int, std::string, int, int, int[3], int[3], bool);

  bool             (*variableExists)(DataArchive*, std::string);

  ParticleDataRaw* (*getParticleData)(DataArchive*, GridP*, int, int, std::string, int, int);

  std::string      (*getParticlePositionName)(DataArchive*);
};
#endif
