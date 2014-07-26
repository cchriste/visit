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
//                            avtMirandaFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_Miranda_FILE_FORMAT_H
#define AVT_Miranda_FILE_FORMAT_H

#include <avtMTMDFileFormat.h>

#include <string>
#include <vector>

class DBOptionsAttributes;


// ****************************************************************************
//  Class: avtMirandaFileFormat
//
//  Purpose:
//      Reads in Miranda files as a plugin to VisIt.
//
//  Programmer: dbremer -- generated by xml2avt
//  Creation:   Tue Jan 23 17:00:13 PST 2007
//
//  Modifications:
//    Dave Bremer, Wed Feb 20 15:25:12 PST 2008
//    Added support for version 1.2 of the .raw files, which specifies block
//    ordering using a tag of the form "fileorder: ZYX", rather than using a
//    separate grid file per process.
//
//    Dave Bremer, June 25, 2009
//    Added support for curvilinear meshes.  If curvilinear, the metadata file
//    will have the tag "curvilinear: yes".  Curvilinear data will be node
//    centered, while rectilinear data will remain zone centered.  The files
//    do not contain any duplicate points along the boundaries, so for each
//    mesh, var, vec, or mat read, neighboring domains be need to be opened 
//    and ghost nodes extracted.
// ****************************************************************************

class avtMirandaFileFormat : public avtMTMDFileFormat
{
  public:
                           avtMirandaFileFormat(const char *, DBOptionsAttributes *);
    virtual               ~avtMirandaFileFormat() {;};

    virtual void          *GetAuxiliaryData(const char *var, int timestate, 
                                            int domain, const char *type, 
                                            void *args, DestructorFunction &);
    
    virtual void           GetCycles(std::vector<int> &c)    {c = aCycles;}
    virtual void           GetTimes(std::vector<double> &t)  {t = aSimTimes;}

    virtual int            GetNTimesteps(void)  {return (int)aCycles.size();}

    virtual const char    *GetType(void)   { return "Miranda"; }
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, int, const char *);
    virtual vtkDataArray  *GetVar(int, int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, int, const char *);

  protected:
    std::string                 sFileVersion; 
    int                    dim;           // dimensionality.  2 or 3
    int                    flatDim;       // if dim==2, flatDim == 0,1,or 2
    

    double                 fOrigin[3];    // sample 0,0,0 is centered on fOrigin
    double                 fStride[3];    // spacing between samples
    
    int                    iGlobalDim[3];  // dimension of entire domain

    int                    iBlockSize[3]; // size of miranda compute domains per processor -- for version 2.0, this is smaller than the actual output to disk because it has zero shared nodes with adjacent blocks
    int                    iInteriorSize[3]; // format 2.0:  size of interior viz blocks written to disk for blocks 0 to iNumBlocks[i-2]
    int                    iBoundarySize[3];// format 2.0:  size of viz blocks written to disk for block iNumBlocks[i-1] (upper edge of domain)
    int                    iNumBlocks[3];
    
    std::vector<std::string> aVarNames;
    std::vector<int>         aVarNumComps;
    std::vector<std::vector<float> > aVarMinMax;    // Current Miranda dumps write min/max
                                          // of last timestep, so for now this 
                                          // is almost unused.

    std::vector<std::string> aMatNames;     // Material names. Size is 0 if there
                                          // is no material set. If present, the
                                          // material set comes after all vars.
    std::vector<int>       aCycles;
    std::vector<double>    aSimTimes;
    
    std::string            fileTemplate;
    std::string            gridTemplate;

    std::vector<int>       domainMap;     // has 3-int tuple for each block--
                                          // domain index to i,j,k block.

    int                    iFileOrder[3]; // 2,1,0; or 0,1,2; or another permutation
                                          // maps a domain index to i,j,k block
                                          // replaces gridTemplate and domainMap

    bool                   bCurvilinear;
    bool                   bZonal;

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *, int);
    virtual void           DomainToIJK(int domain, int &i, int &j, int &k);

    void                   GetBlockDims(int domain, int dims[3]);
    virtual vtkDataSet    *GetCurvilinearMesh2(int domain);
    virtual vtkDataSet    *GetCurvilinearMesh(int domain);
    virtual vtkDataSet    *GetRectilinearMesh(int domain);
     void                   InterleaveData(float *__restrict dst, float *__restrict src, int *dstDim, int nComp);
    virtual void           PackData(float *__restrict dst, const  float * const  *__restrict src, 
                                    const int *dstDim, int nComp, bool interleave);
    virtual void           ReadRawScalar(FILE *fd, int iComp, float *out, const char *filename, int domain);
    virtual void           FindNeighborDomains(int domain, int *neighbors, int *realdim);

    static  void           SkipToEndOfLine( ifstream &f, bool bCheckForBadTokens = true );
    static  double         GetFortranDouble( ifstream &f );
};


#endif
