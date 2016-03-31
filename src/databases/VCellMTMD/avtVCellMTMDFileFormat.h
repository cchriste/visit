/*****************************************************************************
*
* Copyright (c) 2000 - 2016, Lawrence Livermore National Security, LLC
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

/*
 * Virtual Cell (VCellMTMD) specific portions of this file are
 * Copyright (C) 1999-2011 University of Connecticut Health Center
 *
 * Licensed under the MIT License (the "License").
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at:
 *
 *  http://www.opensource.org/licenses/mit-license.php
 */

// ************************************************************************* //
//                            avtVCellMTMDFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_VCellMTMD_FILE_FORMAT_H
#define AVT_VCellMTMD_FILE_FORMAT_H

#include <avtMTMDFileFormat.h>

#include <avtMaterial.h>
#include <avtVarMetaData.h>
#include <vtkDoubleArray.h>

#include <vector>
#include <string>
#include <map>


static const std::string VOLMESH = std::string("volMesh");
static const std::string MEMBRMESH = std::string("membrMesh");
static const std::string POINTMESH = std::string("pointMesh");

static const std::string COMPARTMENTSUBDOMAIN_MATERIAL = std::string("compartmentSubDomains");
static const std::string MEMBRANESUBDOMAIN_MATERIAL = std::string("membraneSubDomains");



// ****************************************************************************
//  Class: avtVCellMTMDFileFormat
//
//  Purpose:
//      Reads in VCellMTMD files as a plugin to VisIt.
//
//  Programmer: frm -- generated by xml2avt
//  Creation:   Sat Oct 23 13:29:29 PST 2010
//
// ****************************************************************************

struct VCellLogEntry {
    int iteration;
    std::string simFileName;
    std::string zipFileName;
    double time;
};

struct MembraneRegionsMapVolumeRegion {
    int membregid;
    int volReg1;
    int volReg2;
};

struct MembraneElements{
    int membregid;
    int volIndex0;
    int volIndex1;
};


struct VCellMeshInfo {
    private:
        double* extents;
    public:
        int numElementsXYZ[3];
        float extentXYZ[3];
        float originXYZ[3];
        bool bInit;
        int * volumeRegionsMapSubvolume;
        unsigned short * volumeElementsMapVolumeRegion;
        MembraneRegionsMapVolumeRegion  * membraneRegionsMapVolumeRegion;
        std::map<std::string,int> compartmentSubdomainMap;
        std::map<std::string,std::map<std::string,int> > membraneSubdomainMap;
        int numMembraneElements;
        MembraneElements * membraneElements;
        bool isSmoldyn;


    VCellMeshInfo(){
        extents = NULL;
        bInit = false;
        numElementsXYZ[0] = 1;numElementsXYZ[1] = 1;numElementsXYZ[2] = 1;
        extentXYZ[0] = 1;extentXYZ[1] = 1;extentXYZ[2] = 1;
        originXYZ[0] = 0;originXYZ[1] = 0;originXYZ[2] = 0;
        compartmentSubdomainMap.clear();
        membraneSubdomainMap.clear();
    }
    ~VCellMeshInfo(){
        if(extents != NULL){delete [] extents;}
        if(volumeRegionsMapSubvolume != NULL){delete [] volumeRegionsMapSubvolume;}
        if(membraneRegionsMapVolumeRegion != NULL){delete [] membraneRegionsMapVolumeRegion;}
        if(volumeElementsMapVolumeRegion != NULL){delete [] volumeElementsMapVolumeRegion;}
        if(membraneElements != NULL){delete [] membraneElements;}
    }
    int getSpatialdimension(){
        return 1 + (numElementsXYZ[1]>1) + (numElementsXYZ[2]>1);
    }
    double* getExtents(){
        if(extents == NULL){
            extents = new double[6];

            extents[0] = originXYZ[0];
            extents[1] = extentXYZ[0];

            extents[2] = originXYZ[1];
            extents[3] = extentXYZ[1];

            extents[4] = originXYZ[2];
            extents[5] = extentXYZ[2];
        }
        return extents;
    }
};

struct CoordinateIndex {
    int x,y,z;
    CoordinateIndex(){
        x = 0;
        y = 0;
        z = 0;
    }
};

struct Coordinate {
    float x,y,z;
    Coordinate(){
        x = 0;
        y = 0;
        z = 0;
    }
    Coordinate(float px,float py,float pz){
        x = px;
        y = py;
        z = pz;
    }
};

typedef enum {
    VAR_UNKNOWN =            0,
    VAR_VOLUME =            1,
    VAR_MEMBRANE =            2,
    VAR_CONTOUR =            3,
    VAR_VOLUME_REGION =        4,
    VAR_MEMBRANE_REGION=    5,
    VAR_CONTOUR_REGION =    6
} VariableType;


typedef int int32;
typedef unsigned int uint32;

#define DATABLOCK_STRING_SIZE  124
#define MAGIC_STRING "VCell Data Dump"
#define VERSION_STRING  "2.0.1  "

struct FileHeader {
    char   magicString[16];
    char   versionString[8];
    int32  numBlocks;
    int32   firstBlockOffset;
    int32   sizeX;
    int32   sizeY;
    int32   sizeZ;
};
struct DataBlock {
    char   varName[DATABLOCK_STRING_SIZE];
    int32   varType;
    int32   size;
    int32   dataOffset;
    std::string domainName;
};

struct MembraneSubstVariables {
    std::string convertVolToMembrName;
    std::string membrFuncDomainName;
    std::vector<int> volumeIndexes;
};

struct RegionSubstVariables {
    std::string regionSubstName;
    int varType;
    std::string domainName;
};

class avtVCellMTMDFileFormat : public avtMTMDFileFormat
{
    private:
        std::vector<VCellLogEntry> vcellLogEntryList;
        std::vector<DataBlock> variableNames;
        std::vector<MembraneSubstVariables> membrSubstVariables;
        std::vector<RegionSubstVariables> regionSubstVariables;
        std::map<std::string,std::string> funcMapDomain;
        std::string baseFileName;
        std::string baseDirName;
        VCellMeshInfo vcellMeshInfo;
        bool bSimZip;
        int * membrQuadNodes;
        Coordinate * ucdMeshNodes;
        vtkDoubleArray *smoldynDataValues;

  public:
                       avtVCellMTMDFileFormat(const char *);
    virtual           ~avtVCellMTMDFileFormat();

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
     virtual void      *GetAuxiliaryData(const char *var, int timestep,int domain, const char *type, void *args,DestructorFunction &);
    

    //
    // If you know the times and cycle numbers, overload this function.
    // Otherwise, VisIt will make up some reasonable ones for you.
    //
    // virtual void        GetCycles(std::vector<int> &);
    virtual void        GetTimes(std::vector<double> &);
    //

    virtual int            GetNTimesteps(void);

    virtual const char    *GetType(void)   { return "VCellMTMD"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, int, const char *);
    virtual vtkDataArray  *GetVar(int, int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, int, const char *);


        //Added by frm-----------------------
    void        readLog(std::string& logFileName);
    void        readVCellMeshFile(void);
    void        readStateVariables();
    bool        zipUnzipWithRetry(const char* zipFileName,const char* simFileName, char* errmsg);
    void        readSimFileHeader(char *filename);
    void        calcUCDMembraneQuads();
    int            calcNodeCount(int axis);
    void        getCoordinateIndexFromVolumeIndex(CoordinateIndex& holder,int volIndex);
    void        calcUCDGridNodes();
    double        coordComponentFromSinglePlanePolicy(int argAxisFlag);
    void        getCoordinate(Coordinate& holder,CoordinateIndex coordIndex);
    void        readFunctions(avtDatabaseMetaData *md);
    void        readVariableValues(const char * simFileName,std::string varname,vtkDoubleArray * dataHolder);
    void        readDoubles(FILE *fp, double *data, int length);
    void        readVolumeSamples(ifstream & ifs);
    unsigned char fromHex(unsigned char* src);
    void        substituteMembraneFunctions(avtDatabaseMetaData *md);
    void        readSubDomains();
    std::string        findDomainName(std::string & varOrFuncName);
    void        createSubDomainMaterials(avtDatabaseMetaData *md);
    void        substituteVisitFunctionSyntax(avtDatabaseMetaData *md);
    void        convertRegionToTempVar(std::string & vcellFuncStr,std::string & domainName);
    void        addRegionSubstVar(std::string & regionSubstVarName,int varType,std::string & domainName);
    void        setMaterialRestricted(avtDatabaseMetaData *md,avtVarMetaData * varMetaData);
    void        stripToBaseName(std::string & filePath);
    //-----------------------------------

  protected:
    // DATA MEMBERS

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *, int);
};


#endif
