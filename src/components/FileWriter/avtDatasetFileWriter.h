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
//                           avtDatasetFileWriter.h                          //
// ************************************************************************* //

#ifndef AVT_DATASET_FILE_WRITER_H
#define AVT_DATASET_FILE_WRITER_H

#include <file_writer_exports.h>

#include <visitstream.h>
#include <string>
#include <vector>

#include <avtOriginatingDatasetSink.h>


typedef enum
{
    CURVE                  = 0,
    OBJ,                  /* 1 */
    STL,                  /* 2 */
    VTK,                  /* 3 */
    ULTRA                 /* 4 */
} DatasetFileFormat;


// ****************************************************************************
//  Class: avtDatasetFileWriter
//
//  Purpose:
//      A type of dataset sink that writes the dataset to the specified file
//      format.
//
//  Programmer: Hank Childs
//  Creation:   May 24, 2002
//
//  Modifications:
//    Jeremy Meredith, Sat Apr 12 15:09:28 PDT 2003
//    Added the Ultra file format.
//
//    Jeremy Meredith, Tue Dec 30 09:16:12 PST 2003
//    Removed the obsolete Curve file format.  Renamed ULTRA to Curve.
//
//    Hank Childs, Thu Feb  5 17:11:06 PST 2004
//    Moved inlined destructor definition to .C file because certain compilers 
//    have problems with them.
//
//    Brad Whitlock, Mon Mar 6 17:35:28 PST 2006
//    I made it reset nFilesWritten if the nase changes.
//
// ****************************************************************************

class AVTFILEWRITER_API avtDatasetFileWriter : public avtOriginatingDatasetSink
{
  public:
                       avtDatasetFileWriter();
    virtual           ~avtDatasetFileWriter();

    void               Write(DatasetFileFormat, const char *filename, bool);

    char              *CreateFilename(const char *base, bool family,
                                      DatasetFileFormat);

  protected:
    static const char *extensions[];
    int                nFilesWritten;
    char              *oldFileBase;

    void               WriteSTLFile(const char *, bool);

    void               WriteCurveFile(const char *);

    void               WriteVTKFamily(const char *, bool);
    int                WriteVTKTree(avtDataTree_p, int, const char *, bool);
    void               WriteVTKFile(vtkDataSet *, const char *, bool);

    void               WriteOBJFamily(const char *);
    int                WriteOBJTree(avtDataTree_p, int, const char *);
    void               WriteOBJFile(vtkDataSet *, const char *, const char *);

    vtkDataSet        *GetSingleDataset(void);
    char              *GenerateName(const char *, const char *,
                                    std::vector<std::string> &);
};


#endif


