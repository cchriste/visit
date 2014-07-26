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
//                             avtTecplotWriter.h                            //
// ************************************************************************* //

#ifndef AVT_TECPLOT_WRITER_H
#define AVT_TECPLOT_WRITER_H

#include <avtDatabaseWriter.h>
#include <DBOptionsAttributes.h>

#include <string>
#include <vector>
#include <visitstream.h>

class vtkPoints;
class vtkPolyData;
class vtkRectilinearGrid;
class vtkStructuredGrid;
class vtkUnstructuredGrid;

// ****************************************************************************
//  Class: avtTecplotWriter
//
//  Purpose:
//      A module that writes out Tecplot files.
//
//  Programmer: Jeremy Meredith
//  Creation:   February 15, 2005
//
//  Modifications:
//
//    Hank Childs, Tue Sep 27 10:21:36 PDT 2005
//    Use virtual inheritance.
//
//    Jeremy Meredith, Tue Mar 27 17:03:47 EDT 2007
//    Added numblocks (currently ignored) to the OpenFile interface.
//
//    Brad Whitlock, Wed Sep  2 14:16:43 PDT 2009
//    I added methods for writing rectilinear and polydata datasets.
//
// ****************************************************************************

class avtTecplotWriter : public virtual avtDatabaseWriter
{
  public:
                   avtTecplotWriter(DBOptionsAttributes *);
    virtual       ~avtTecplotWriter();

  protected:
    std::string    stem;
    ofstream       file;

    virtual void   OpenFile(const std::string &, int);
    virtual void   WriteHeaders(const avtDatabaseMetaData *,
                                std::vector<std::string> &, 
                                std::vector<std::string> &,
                                std::vector<std::string> &);
    virtual void   WriteChunk(vtkDataSet *, int);
    virtual void   CloseFile(void);

  private:
    void           WritePolyData(vtkPolyData *, int);
    void           WriteRectilinearMesh(vtkRectilinearGrid *, int);
    void           WriteCurvilinearMesh(vtkStructuredGrid *, int);
    void           WriteUnstructuredMesh(vtkUnstructuredGrid *, int);

    void           WritePoints(vtkPoints *pts, int dim);
    void           WriteDataArrays(vtkDataSet*);
    void           WriteVariables(const std::vector<std::string> &);

    bool           ReallyHasMaterials();

    std::vector<std::string> variableList;
    std::vector<std::string> materialList;
    bool variablesWritten;
};


#endif
