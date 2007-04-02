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
//                             avtTecplotWriter.h                              //
// ************************************************************************* //

#ifndef AVT_Tecplot_WRITER_H
#define AVT_Tecplot_WRITER_H

#include <avtDatabaseWriter.h>

#include <string>
#include <vector>
#include <visitstream.h>
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
// ****************************************************************************

class
avtTecplotWriter : public virtual avtDatabaseWriter
{
  public:
                   avtTecplotWriter();
    virtual       ~avtTecplotWriter();

  protected:
    std::string    stem;
    ofstream       file;

    virtual void   OpenFile(const std::string &);
    virtual void   WriteHeaders(const avtDatabaseMetaData *,
                                std::vector<std::string> &, 
                                std::vector<std::string> &,
                                std::vector<std::string> &);
    virtual void   WriteChunk(vtkDataSet *, int);
    virtual void   CloseFile(void);

  private:
    void           WriteCurvilinearMesh(vtkStructuredGrid *, int);
    void           WriteUnstructuredMesh(vtkUnstructuredGrid *, int);
    void           WriteDataArrays(vtkDataSet*);
    std::vector<std::string> variableList;
    bool hadMaterial;
};


#endif
