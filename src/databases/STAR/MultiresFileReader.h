/*****************************************************************************
*
* Copyright (c) 2010, University of New Hampshire Computer Science Department
* All rights reserved.
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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  File:        MultiresFileReader.h                                        //
//  Programmer:  Andrew Foulks <rafoulks@cs.unh.edu>                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef _MULTI_RES_FILE_READER_H_
#define _MULTI_RES_FILE_READER_H_

#include <string>
#include <vector>

#define DEBUG_TO_STDERR

#include "DataManagerAPI.h"

/**
 *      The MultiresFileReader is a multiresolution file format
 *      being developed for a VisIt database plugin.  
 *      It is a type of an archive file which contains a 
 *      header in ascii followed by binary data.  The header 
 *      describes how many resolutions are in the file, and 
 *      where each resolution is located, as well as each 
 *      resolution's size.  All data is 32-bit floating point.  
 *      All data is assumed to be 3D.
 *
 *      Data Structure:
 *          mDataBuffer points to the data in the entire file, 
 *          which includes all resolutions.
 *
 *          mRawDataPtrs is an array of pointers, each pointer
 *          pointing to the start of its resolution within
 *          the DataBuffer.
 *
 *              mDataBuffer ------+
 *                                |
 *                                v
 *                               +-------------------------+
 *              mRawDataPtrs[0]=>|      resolution 0       |
 *                               |                         |
 *                               |                         |
 *                               +-------------------------+
 *              mRawDataPtrs[1]=>|      resolution 1       |
 *                               |                         |
 *                               +-------------------------+
 *              mRawDataPtrs[2]=>|      resolution 2       |
 *                               +-------------------------+
 *                   ...         |      ... etc ...        |
 *
 *
 *      @author Andrew Foulks
 **/

class MultiresFileReader : public DataManagerAPI
{
public:  // 'structors
            MultiresFileReader();
            MultiresFileReader(const char* filename);
    virtual ~MultiresFileReader();

public:  // api implemented from MetadataAPI
    virtual const char*     filename        () const {return mFilename.c_str();}
    virtual float           fileVersion     () const {return mFileVersion;}
    virtual int             numFiles        () const {return 1;}
    virtual int             numResolutions  () const {return mNumResolutions;}
    virtual std::string     gridFilename    () const;
    virtual std::vector<int> timeStepList    () const;

    // data rank: scalars, vectors, or tensors
    virtual int             numVariables        () const {return 1;}
    virtual std::string     variableNameAt      (int index) const;
    virtual int             indexOfVariableName (const std::string& name="") const
                                                    {return 0;}
    virtual bool            isScalar            (const std::string& name="") const;
    virtual char            isVectorComponent   (const std::string& name="") const;
    virtual bool            isVector            (const std::string& name="") const;
    virtual bool            isTensor            (const std::string& name="") const;
    virtual int             numVectorExpressions() const {return 0;}
    virtual std::string     vectorExpressionAt  (int index) const {return "";}

    // data size and structure
    virtual int             width               (int resolution) const;
    virtual int             height              (int resolution) const;
    virtual int             depth               (int resolution) const;
    virtual int             numXchunks          (int resolution) const;
    virtual int             numYchunks          (int resolution) const;
    virtual int             numZchunks          (int resolution) const;
    virtual bool            hasMinVal           () const {return mHasMinVal;}
    virtual bool            hasMaxVal           () const {return mHasMaxVal;}
    virtual float           minVal              () const {return mMinVal;}
    virtual float           maxVal              () const {return mMaxVal;}
    virtual const float*    rawData             (const std::string& varName,
                                                 int resolution,
                                                 int fileIndex,
                                                 int chunkIndex);
    virtual void            freeRawDataMemory   (const std::string& varName="",
                                                 int fileIndex=0);

public:
    virtual void            parseFile           (const char* filename);
    virtual const float*    rawData             (int res, int chunkIndex)
                        { return rawData(mVariableName, res, 0, chunkIndex); }

private: // helper functions
    void        loadDataFromFile        ();
    float       parseVersionNumber      (const char* line);
    std::string parseGridFilename       (const char* line);
    std::string parseDataFilename       (const char* line);
    std::string parseDataType           (const char* line);
    std::string parseDataRank           (const char* line);
    std::string parseVariableName       (const char* line);
    std::vector<int> parseNumXchunks    (const char* line);
    std::vector<int> parseNumYchunks    (const char* line);
    std::vector<int> parseNumZchunks    (const char* line);
    float       parseMinVal             (const char* line);
    float       parseMaxVal             (const char* line);
    int         parseNumResolutions     (const char* line);
    int         parseNumErrorDataSets   (const char* line);
    void        parseResolutionMetadata (const char* line,
                                            int& resolution,
                                            int& width,
                                            int& height,
                                            int& depth,
                                            long long& offset);
    int         parseHeaderSize         (const char* line);
    const char* getNextLine             (FILE* f);
    void        checkFilePtr            (FILE* ptr);

// test functionality
public: static int main(int argc, char** argv);

private: // members
    float*          mDataBuffer;        // ptr to all data resolutions
    std::vector<long long> mOffsets;         // byte offset to each resolution
    std::vector<float*>  mRawDataPtrs;       // ptr to data array for each res
    std::vector<int> mWidths;            // size (width) for each res
    std::vector<int> mHeights;           // size (height) for each res
    std::vector<int> mDepths;            // size (depth) for each res
    std::vector<int> mChunkSizes;        // num elements per chunk per res
    int             mHeaderSize;        // size (bytes) of file metadata
    std::string     mFilename;          // name of 'mrd' file.
    std::string     mGridFilename;      // name of grid file in header
    std::string     mDataFilename;      // data file name if different
    int             mNumResolutions;    // number of resolutions in the file
    int             mNumErrorDataSets;  // number of error data resolutions 
    float           mFileVersion;       // version of the 'mrd' file read
    std::string     mDataType;          // must be "float", else bad will be
    std::string     mDataRank;          // "scalar", "vector", "tensor"
    std::vector<int> mNumXchunks;        // number of chunks along 'width' axis
    std::vector<int> mNumYchunks;        // number of chunks along 'height' axis
    std::vector<int> mNumZchunks;        // number of chunks along 'depth' axis
    std::string     mVariableName;      // scalar or vector name
    float           mMinVal;            // min of entire MR dataset
    float           mMaxVal;            // max val of entire MR dataset
    bool            mHasMinVal;         // flag indicates if min val in file
    bool            mHasMaxVal;         // flag indicates if max val in file

private: // disable copy constructor and operator=
    MultiresFileReader(MultiresFileReader&);
    MultiresFileReader& operator=(MultiresFileReader&);
};

#endif // _MULTI_RES_FILE_READER_H_

