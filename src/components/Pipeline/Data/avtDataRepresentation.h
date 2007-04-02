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
//                       avtDataRepresentation.h                             //
// ************************************************************************* //

#ifndef AVT_DATAREPRESENTATION_H
#define AVT_DATAREPRESENTATION_H
#include <pipeline_exports.h>


#include <array_ref_ptr.h>
#include <string>

class  vtkDataSet;


typedef enum
{
    DATASET_TYPE_RECTILINEAR,
    DATASET_TYPE_CURVILINEAR,
    DATASET_TYPE_UNSTRUCTURED,
    DATASET_TYPE_POLYDATA,
    DATASET_TYPE_NULL,
    DATASET_TYPE_UNKNOWN
}  DataSetType;


// ****************************************************************************
//  Class: avtDataRepresentation
//
//  Purpose:
//      Buffers how data is being stored - as a vtkDataSet or as a char
//      string.
//
//  Programmer: Hank Childs
//  Creation:   September 22, 2000
//
//    Jeremy Meredith, Tue Sep 26 09:28:46 PDT 2000
//    Added constructors from VTK datasets and strings.  Added copy 
//    constructor.  Removed SetDomain/SetDomainString and Merge.
//    Added assignment operator.
//
//    Kathleen Bonnell, Thu Apr  5 13:27:02 PDT 2001 
//    Changed class name to avtDataRepresentation.
//
//    Hank Childs, Mon Jun  4 14:02:19 PDT 2001
//    Added domain identifiers with the data representation.
//
//    Hank Childs, Mon Sep 17 09:11:23 PDT 2001
//    Changed way the string representation is managed to prevent unnecessary
//    copies.  Removed unused method CreateDataString.
//
//    Kathleen Bonnell, Wed Sep 19 13:45:33 PDT 2001 
//    Added member 'label' and retrieval method. Modifed constructor to
//    include the label.
//    
//    Hank Childs, Wed Mar 17 19:16:48 PST 2004
//    Added a static null dataset to prevent unnecessary memory usage with SR.
//
//    Mark C. Miller, Wed Nov 16 13:31:19 PST 2005
//    Added compression
// ****************************************************************************

class PIPELINE_API avtDataRepresentation
{
  public:
                        avtDataRepresentation();
                        avtDataRepresentation(vtkDataSet *, int, std::string,
                                              bool dontCopyData = false);
                        avtDataRepresentation(char *, int, int, std::string,
                                              CharStrRef &, DataSetType);
                        avtDataRepresentation(const avtDataRepresentation &);
    virtual            ~avtDataRepresentation();

    const avtDataRepresentation    &operator=(const avtDataRepresentation &);

    vtkDataSet         *GetDataVTK(void);
    unsigned char      *GetDataString(int &, DataSetType &);
    unsigned char      *GetCompressedDataString(int &, DataSetType &);
    int                 GetNumberOfCells(int topoDim, bool polysOnly) const;

    float               GetCompressionRatio() const;
    float               GetTimeToCompress() const;
    float               GetTimeToDecompress() const;

    bool                Valid(void);

    int                 GetDomain(void) { return domain; };
    std::string         GetLabel(void)  { return label; };

  protected:
    vtkDataSet         *asVTK;
    unsigned char      *asChar;
    int                 asCharLength;
    CharStrRef          originalString;
    DataSetType         datasetType;

    float               compressionRatio;
    float               timeToCompress;
    float               timeToDecompress;

    int                 domain;
    std::string         label;

    static bool         initializedNullDataset;
    static vtkDataSet  *nullDataset;

    unsigned char      *GetDataString(int &, DataSetType &, bool);
    static void         InitializeNullDataset(void);
    static DataSetType  DatasetTypeForVTK(vtkDataSet *);
};


#endif


