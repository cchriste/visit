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
//                            avtMFIXFileFormat.h                            //
// ************************************************************************* //

#ifndef AVT_MFIX_FILE_FORMAT_H
#define AVT_MFIX_FILE_FORMAT_H

#include <avtMTMDFileFormat.h>

#include <vectortypes.h>

#ifdef _WIN32
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
  #define MFIX_NAME_MAX MAX_PATH
#else
  #include <limits.h>
  #ifdef __APPLE__
    #include <sys/syslimits.h>
  #endif
  #define MFIX_NAME_MAX (NAME_MAX+1)
#endif

class DBOptionsAttributes;
class vtkDoubleArray;
class vtkFloatArray;
class vtkIntArray;
class vtkLongLongArray;
class vtkMFIXReader;
class vtkStringArray;
class vtkUnstructuredGrid;

// ****************************************************************************
//  Class: avtMFIXFileFormat
//
//  Purpose:
//      Reads in MFIX files as a plugin to VisIt.
//
//  Programmer: bdotson -- generated by xml2avt
//  Creation:   Fri May 26 08:59:22 PDT 2006
//
//  Modifications:
//    Jeremy Meredith, Thu Aug  7 15:46:38 EDT 2008
//    Use a C++ string for Version so comparisons are easier.
//
//    Kathleen Bonnell, Wed Oct 22 17:11:07 PDT 2008
//    Reworked to use vtkMFIXReader for bulk of work, and to have DBOptions.
//
// ****************************************************************************

class avtMFIXFileFormat : public avtMTMDFileFormat
{
 public:
     avtMFIXFileFormat(const char *, DBOptionsAttributes *);
     virtual               ~avtMFIXFileFormat();

     virtual int            GetNTimesteps(void);

     virtual void           GetTimes(doubleVector &_times);

     virtual const char    *GetType(void)   { return "MFIX"; };
     virtual void           FreeUpResources(void);

     virtual vtkDataSet    *GetMesh(int, int, const char *);
     virtual vtkDataArray  *GetVar(int, int, const char *);
     virtual vtkDataArray  *GetVectorVar(int, int, const char *);
     virtual void          *GetAuxiliaryData(const char *, int, int,
                                const char *, void *, DestructorFunction &);

     int                    get_var_index(const char *);
     void                   make_fn(char *, int);

 protected:

     enum DataType {
         DATATYPE_INT,
         DATATYPE_FLOAT,
         DATATYPE_DOUBLE,
         DATATYPE_CHAR
     };

     virtual void PopulateDatabaseMetaData(avtDatabaseMetaData *, int);

 private:
     void               ReadInformation(void);
     void               SetUp(void);
     void               ReadRestartFile(void);
     vtkDataSet*        BuildMesh(int xDomain, int yDomain, int zDomain);
     void               RestartVersionNumber(const char* buffer);
     void               GetInt(istream& in, int &val);
     void               GetDouble(istream& in, double& val);
     void               SkipBytes(istream& in, int n);
     void               SwapDouble(double &value);
     void               SwapFloat(float &value);
     void               SwapInt(int &value);
     void               GetBlockOfDoubles(istream& in, vtkDoubleArray *v, int n);
     void               GetBlockOfInts(istream& in, vtkIntArray *v, int n);
     void               SkipBlockOfInts(istream& in, int n);
     void               CreateVariableNames(void);
     void               GetTimeSteps(void);
     void               CalculateMaxTimeStep(void);
     void               MakeTimeStepTable(int nvars);
     void               GetNumberOfVariablesInSPXFiles(void);
     void               MakeSPXTimeStepIndexTable(int nvars);
     void               GetAllTimesTweaked(void);
     void               GetSubBlock(void * buf, const char *fname,
                            DataType dataType,
                            long long startOffsetBytes,
                            int startSkipVals,
                            int xVals, int xStrideVals,
                            int yVals, int yStrideVals,
                            int zVals);

     void               CalcDomainBreakdown2D(long targetDomains,
                            int cellsX, int cellsY, int* nX, int* nY);
     void               CalcDomainBreakdown3D(long targetDomains,
                            int cellsX, int cellsY, int cellsZ,
                            int* nX, int* nY, int* nZ);
     void               decompose_domains(int, int *, int *, int *);
     void               get_limit(int, int, int, vtkDoubleArray *,
                            int *, int *, double **);

     char*              GSB_currentFname;
     FILE*              GSB_file;
     int                SwapByteOrder;
     char               RestartFileName[MFIX_NAME_MAX];
     doubleVector       times;
     stringVector       variables;
     bool               readInData;
     bool               readInformation;
     bool               fileBigEndian;
     bool               setupCompleteFlag;
     int                numDomainsPerProc;
     int                numXDomains;
     int                numYDomains;
     int                numZDomains;

     int                par_rank;
     int                par_size;

     // The following have been incorporated from vtkMFIXReader.
     char               DataBuffer[513];
     char               Version[120];
     float              VersionNumber;
     int                IMinimum1;
     int                JMinimum1;
     int                KMinimum1;
     int                IMaximum;
     int                JMaximum;
     int                KMaximum;
     int                IMaximum1;
     int                JMaximum1;
     int                KMaximum1;
     int                IMaximum2;
     int                JMaximum2;
     int                KMaximum2;
     int                IJMaximum2;
     int                IJKMaximum2;
     int                MMAX;
     int                TimeStep;
     //int ActualTimeStep;
     int                CurrentTimeStep;
     int                NumberOfTimeSteps;
     int               *TimeSteps;
     int                TimeStepRange[2];
     int                TimeStepWasReadOnce;
     double             DeltaTime;
     double             XMinimum;
     double             XLength;
     double             YLength;
     double             ZLength;
     int                DimensionIc;
     int                DimensionBc;
     int                DimensionC;
     int                DimensionIs;
     double             Ce;
     double             Cf;
     double             Phi;
     double             PhiW;
     char               CoordinateSystem[17];
     char               Units[17];
     int                NumberOfSPXFilesUsed;
     int                NumberOfScalars;
     int                NumberOfReactionRates;
     bool               BkEpsilon;
     int                MaximumTimestep;       // maximum timesteps amongst
                                               // the variables
     int                SPXRecordsPerTimestep; // number of records in a single
                                               // timestep for a variable
     char               FlagFilename[MFIX_NAME_MAX];     // What file are the flags in?
     streampos          FlagFieldOffset;       // Where in the flag file is
                                               // that array?

     vtkIntArray       *NMax;                  // Array to hold number of
                                               // species per phase
     vtkDoubleArray    *C;                     // Array used to parse restart
                                               // file
     vtkDoubleArray    *Dx;                    // Cell widths in x axis
     vtkDoubleArray    *Lx;                    // Integral of Dx
     vtkDoubleArray    *Dy;                    // Cell widths in y axis
     vtkDoubleArray    *Ly;                    // Integral of Ly
     vtkDoubleArray    *Dz;                    // Cell widths in z axis
     vtkDoubleArray    *Lz;                    // Integral of Lz
     vtkDoubleArray    *TempD;                 // Array used to parse restart
                                               // file
     vtkIntArray       *TempI;                 // Array used to parse restart
                                               // file
     vtkIntArray       *SpxFileExists;         // Array for keeping track of
                                               // what spx files exist.
     vtkStringArray    *VariableNames;
     vtkIntArray       *VariableComponents;
     vtkIntArray       *VariableIndexToSPX;    // This gives the spx file number
                                               // for the particular variable.
     vtkIntArray       *VariableTimesteps;     // number of timesteps for each
                                               // variable
     vtkIntArray       *VariableTimestepTable; // Since the number of timesteps
                                               // vary between variables
                                               //  this is a table that looks
                                               //  up the appropriate timestep
                                               // for the particular variable.
     vtkIntArray       *VariableToSkipTable;   // skip value for each variable,
                                               // this is needed in spx files
                                               // with more than one variable.
     vtkIntArray       *SPXToNVarTable;        // number of variables in each
                                               // spx file
     vtkLongLongArray  *SPXTimestepIndexTable; // This a table look up for the
                                               //  index into a file for a
                                               //  certain variable.

#if 0
     void                 ReadInData(void);
     void                 GetBlockOfFloats(istream& in, vtkFloatArray *v,
         int n);
     vtkIntArray *Flag;            // Cell Flag array
     vtkFloatArray **CellDataArray; // Arrays for variables that will
     //attach to mesh
     vtkPoints *Points;            // Points array for building grid
     vtkUnstructuredGrid *FluidMesh;    // Unstructured Grid
     vtkUnstructuredGrid *InletMesh;    // Unstructured Grid
     vtkUnstructuredGrid *OutletMesh;    // Unstructured Grid
     vtkUnstructuredGrid *ObstructionMesh;    // Unstructured Grid
     vtkHexahedron *AHexahedron;   // Hexahedron type cell
     vtkWedge *AWedge;             // Wedge type cell
     vtkQuad *AQuad;               // Quad type cell

     char FileExtension[15];
     char RunName[256]; */
#endif
};


#endif
