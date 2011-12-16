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
//                            avtXdmfFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_Xdmf_FILE_FORMAT_H
#define AVT_Xdmf_FILE_FORMAT_H

#include <avtMTMDFileFormat.h>

#include <vector>

#include <XdmfObject.h>

class XdmfArray;
class XdmfAttribute;
class XdmfDOM;
class XdmfElement;
class XdmfGrid;

class vtkDataArray;
class vtkImageData;
class vtkRectilinearGrid;
class vtkStructuredGrid;
class vtkUnstructuredGrid;

// ****************************************************************************
//  Class: avtXdmfFileFormat
//
//  Purpose:
//      Reads in Xdmf files as a plugin to VisIt.
//
//  Programmer: kleiter -- generated by xml2avt
//  Creation:   Mon Mar 29 15:43:05 PST 2010
//
//  Modifications:
//
//    Hank Childs, Thu Aug 18 17:22:19 PDT 2011
//    Change timesteps to double (from int).
//
// ****************************************************************************

class avtXdmfFileFormat: public avtMTMDFileFormat
{
    public:
        avtXdmfFileFormat(const char *);
        virtual ~avtXdmfFileFormat();

        //
        // This is used to return unconvention data -- ranging from material
        // information to information about block connectivity.
        //
        // virtual void      *GetAuxiliaryData(const char *var, int timestep,
        //                                     int domain, const char *type, void *args,
        //                                     DestructorFunction &);
        //

        //
        // If you know the times and cycle numbers, overload this function.
        // Otherwise, VisIt will make up some reasonable ones for you.
        //
        // virtual void        GetCycles(std::vector<int> &);
        virtual void        GetTimes(std::vector<double> &);

        virtual int GetNTimesteps(void);

        virtual const char *GetType(void)
        {
            return "Xdmf";
        }
        ;
        virtual void FreeUpResources(void);

        virtual vtkDataSet *GetMesh(int, int, const char *);
        virtual vtkDataArray *GetVar(int, int, const char *);
        virtual vtkDataArray *GetVectorVar(int, int, const char *);

    protected:
        // DATA MEMBERS
        std::string filename;
        std::string firstGrid;
        XdmfGrid * currentGrid;
        XdmfDOM *dom;
        int Stride[3];
        int numGrids;
        std::vector<double> timesteps;

        void AddArrayExpressions(avtDatabaseMetaData *, std::string, std::vector<std::string> &);
        void AddTensorExpressions(avtDatabaseMetaData *, std::string, int, int);
        vtkDataArray* CopyXdmfArray(XdmfArray *, int, int);
        template<typename T> void CopyXdmfArray(XdmfArray *, vtkDataArray *, int, int);
        vtkDataArray * CopyXdmfArrayByPointer(XdmfArray *, int);
        XdmfAttribute * GetAttributeFromName(XdmfGrid *, const char *);
        std::vector<std::string> GetComponentNames(std::string, XdmfInt32, int);
        void GetDims(int[6], int[3]);
        std::string GetFormattedExpressionName(std::string &);
        XdmfGrid * GetGrid(int);
        int GetMeshDataType(XdmfGrid *);
        int GetNumberOfComponents(XdmfGrid *, XdmfAttribute *);
        long GetNumberOfCellComponents(XdmfGrid *, XdmfAttribute *);
        long GetNumberOfNodeComponents(XdmfGrid *, XdmfAttribute *);
        long GetNumberOfPoints(XdmfGrid *);
        int GetNumberOfSymmetricalTensorComponents(int);
        long GetNumberOfValues(XdmfElement *);
        int GetSpatialDimensions(XdmfInt32);
        int GetTopologicalDimensions(XdmfInt32);
        int GetVTKCellType(XdmfInt32);
        bool GetWholeExtent(XdmfGrid *, int[6]);
        virtual void PopulateDatabaseMetaData(avtDatabaseMetaData *, int);
        vtkRectilinearGrid* ReadRectilinearGrid(XdmfGrid *);
        vtkStructuredGrid* ReadStructuredGrid(XdmfGrid *);
        vtkUnstructuredGrid* ReadUnstructuredGrid(XdmfGrid *);
        void ScaleExtents(int[6], int[6], int[3]);
        void SetCurrentGrid(int, const char *);
};

#endif
