/*****************************************************************************
*
* Copyright (c) 2000 - 2007, The Regents of the University of California
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
//                            avtDDCMDFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_DDCMD_FILE_FORMAT_H
#define AVT_DDCMD_FILE_FORMAT_H

#include <avtSTMDFileFormat.h>

#include <string>
#include <vector>

using std::string;
using std::vector;

// ****************************************************************************
//  Class: avtDDCMDFileFormat
//
//  Purpose:
//      Reads in DDCMD files as a plugin to VisIt.
//
//  Programmer: brugger -- generated by xml2avt
//  Creation:   Fri Aug 31 15:27:59 PST 2007
//
// ****************************************************************************

class avtDDCMDFileFormat : public avtSTMDFileFormat
{
  public:
                       avtDDCMDFileFormat(const char *);
    virtual           ~avtDDCMDFileFormat() {;};

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
    // virtual void      *GetAuxiliaryData(const char *var, int domain,
    //                                     const char *type, void *args, 
    //                                     DestructorFunction &);
    //

    virtual bool           ReturnsValidCycle() const { return true; };
    virtual int            GetCycle(void);
    virtual bool           ReturnsValidTime() const { return true; };
    virtual double         GetTime(void);

    void                   ActivateTimestep(void);

    virtual const char    *GetType(void)   { return "DDCMD"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

  protected:
    string                 fname;
    int                    nXFileBlocks, nYFileBlocks, nZFileBlocks, nBlocks;

    // Header information
    unsigned int           lRec, nRecord, nFiles, nFields, swap; 
    double                 hMatrix[9];
    int                    loop;
    double                 time;
    char                 **fieldNames, **fieldTypes;             
    unsigned int          *fieldSizes;
    unsigned int           nXFile, nYFile, nZFile;
    unsigned int           nSpecies;
    char                 **speciesNames;

    // Mesh information
    string                 coordsUnit;
    int                    nDims;
    int                    nXMesh, nYMesh, nZMesh;
    int                    nXMeshBlocks, nYMeshBlocks, nZMeshBlocks;
    float                  xMin, yMin, zMin;
    float                  dX, dY, dZ;

    // Block information
    int                    nXBlock, nYBlock, nZBlock;
    int                    nZonesBlock;
    float                **varsBlock;

    // Variable information
    int                    labelOffset, iSpeciesOffset;
    int                    nVars;
    vector<string>         varNames;
    int                   *varOffsets;
    int                   *varSizes;
    bool                  *varFloat; 

    // File information
    bool                   dataRead;
    int                    headerLength;
    long                  *nRecordsList;
    int                   *fileNumberList;
    off_t                 *fileOffsetList;

    char                  *readProcessorData;

    // Input buffer information
    int                    nInRecords;
    char                  *inProcessorData;

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);

    void                   Convert(void *ptr, int size);
    void                   DetermineBlockDecomposition();
    void                   CopyExchangeDataToBlocks();
    void                   ExchangeProcessorData();
    void                   ReadProcessorChunk();
    void                   DetermineProcessorReadOffset();
    void                   ReadData();
    void                   ReadHeader(const char *filename);
};


#endif
