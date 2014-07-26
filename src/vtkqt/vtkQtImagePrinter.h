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

#ifndef __vtkQtImagePrinter_h
#define __vtkQtImagePrinter_h
#include <vtkqt_exports.h>
#include <QtCore>
#include <QPrinter>
#include <vtkImageWriter.h>

// ****************************************************************************
// Class: vtkQtImagePrinter
//
// Purpose:
//   This is an image file writer class that writes an image to a Qt printer
//   and thus prints the image. Printing was done in this manner to easily
//   fit into the rest of the pipeline while taking advantage of Qt's
//   printing capabilities.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Wed Feb 20 15:46:31 PST 2002
//
// Modifications:
//   
// ****************************************************************************

class VTKQT_API vtkQtImagePrinter : public vtkImageWriter
{
public:
  static vtkQtImagePrinter *New();
  vtkTypeMacro(vtkQtImagePrinter,vtkImageWriter);

  QPrinter &printer() { return print; };
protected:
  vtkQtImagePrinter();
  virtual ~vtkQtImagePrinter() {};
  vtkQtImagePrinter(const vtkQtImagePrinter&) {};
  void operator=(const vtkQtImagePrinter&) {};

  virtual void WriteFile(ofstream *file, vtkImageData *data, int ext[6], int wext[6]);
  virtual void WriteFileHeader(ofstream *, vtkImageData *, int [6]) { };
private:
  QPrinter print;
};

#endif
