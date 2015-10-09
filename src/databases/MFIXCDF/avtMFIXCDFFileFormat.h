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
//                            avtMFIXCDFFileFormat.h                         //
// ************************************************************************* //

#ifndef AVT_MFIXCDF_FILE_FORMAT_H
#define AVT_MFIXCDF_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>
#include <netcdfcpp.h>

#include <string>
#include <set>
#include <map>
#include <memory>

class DBOptionsAttributes;
class vtkDoubleArray;
class vtkFloatArray;
class vtkIntArray;
class vtkLongLongArray;
class vtkMFIXReader;
class vtkStringArray;
class vtkUnstructuredGrid;
class vtkUnsignedCharArray;

// ****************************************************************************
//  Class: avtMFIXCDFFileFormat
//
//  Purpose:
//      Reads in MFIXCDF files as a plugin to VisIt.
//
//  Programmer: welling -- generated by xml2avt
//  Creation:   Wed Aug 3 16:35:13 PST 2011
//
// ****************************************************************************

class avtMFIXCDFFileFormat : public avtSTMDFileFormat
{
public:
    avtMFIXCDFFileFormat(const char *, DBOptionsAttributes *);
    virtual           ~avtMFIXCDFFileFormat();

    virtual double GetTime(void);

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
    virtual void      *GetAuxiliaryData(const char *var, int domain,
                                        const char *type, void *args,
                                        DestructorFunction &);
    //

    //
    // If you know the cycle number, overload this function.
    // Otherwise, VisIt will make up a reasonable one for you.
    //
    // virtual int         GetCycle(void);
    //

    virtual const char    *GetType(void)
    {
        return "MFIXCDF";
    };
    virtual void           FreeUpResources(void);

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

protected:
    // DATA MEMBERS

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);

private:
    void               CalcDomainBreakdown2D(long targetDomains,
                                             int cellsX, int cellsY, int* nX, int* nY);
    void               CalcDomainBreakdown3D(long targetDomains,
                                             int cellsX, int cellsY, int cellsZ,
                                             int* nX, int* nY, int* nZ);
    void               decompose_domains(int, int *, int *, int *);
    void               get_limit(int, int, int, vtkDoubleArray *,
                                 int *, int *, double **);
    void               getBlockOfDoubles1D(NcFile* file,
                                           const char* varname,
                                           double* data,
                                           long offset, int n);
    void               getBlockOfFloats3D(NcFile* file,
                                          const char* varname,
                                          float* data,
                                          long iOffset, int iN,
                                          long jOffset, int jN,
                                          long kOffset, int kN);
    void               getBlockOfInts3D(NcFile* file,
                                        const char* varname,
                                        int* data,
                                        long iOffset, int iN,
                                        long jOffset, int jN,
                                        long kOffset, int kN);
    void               inferVectorVariableNames(avtDatabaseMetaData *md,
                                                std::set<std::string>& varNames);
    void               checkCoordArrays(void);

    // These members are initialized in the constructor
    std::string*       filePath;
    NcFile*            dataFile;
    double             timeNow;
    int                numXDomains;
    int                numYDomains;
    int                numZDomains;
    int                iSz;
    int                jSz;
    int                kSz;
    int                par_rank;
    int                par_size;
    int                coordCode;

    // These members only become valid after execution of
    // PopulateDatabaseMetaData()
    vtkDoubleArray    *Dx;                     // Cell widths in x axis
    vtkDoubleArray    *Lx;                     // Integral of Dx
    vtkDoubleArray    *Dy;                     // Cell widths in y axis
    vtkDoubleArray    *Ly;                     // Integral of Ly
    vtkDoubleArray    *Dz;                     // Cell widths in z axis
    vtkDoubleArray    *Lz;                     // Integral of Lz
};


#endif
