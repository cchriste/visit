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

//
// This code was contributed to the VisIt project by Corvid Technologies
// on February 10, 2010.
//

// ************************************************************************* //
//                            avtVelodyneFileFormat.h                        //
// ************************************************************************* //

#ifndef AVT_Velodyne_FILE_FORMAT_H
#define AVT_Velodyne_FILE_FORMAT_H

#include <vector>
#include <string>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <VelodyneReader.h>
#include <avtSTSDFileFormat.h>


// ****************************************************************************
//  Class: avtVelodyneFileFormat
//
//  Purpose:
//      Reads in Velodyne files as a plugin to VisIt.
//
//  Programmer: hpan -- generated by xml2avt
//  Creation:   Thu Aug 7 11:38:59 PDT 2008
//
// ****************************************************************************

class avtVelodyneFileFormat : public avtSTSDFileFormat
{
public:
  avtVelodyneFileFormat(const char *filename);
  virtual ~avtVelodyneFileFormat();

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
    virtual void  *GetAuxiliaryData(const char *var, const char *type,
                                    void *args, DestructorFunction &);

    //
    // These are used to declare what the current time and cycle are for the
    // file.  These should only be defined if the file format knows what the
    // time and/or cycle is.
    //
    virtual int       GetCycle(void);
    virtual double    GetTime(void);


    virtual const char    *GetType(void)   { return "Velodyne"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(const char *);
    virtual vtkDataArray  *GetVar(const char *);
    virtual vtkDataArray  *GetVectorVar(const char *);


 protected:
    int readNodeIndex();
    int readNodeCoord();
    int readdElements( int grp, int bufsz, int* elmt );
    int readNodeBaseVariableNames();
    int isNodeBasedVariable( const std::string& name );
    //int isVectorVariable( int grp, const char* varname );
    int isTensorVariable( int grp, const char* varname );
    void convert2dVectorTo3dVector( int num, float* val );
    void convert1dVectorTo3dVector( int num, float* val );
    // sym tensor (xx,yy,zz,xy,yz,zx)
    void convertSymTensorToFullTensor( int num, float* val );


  protected:
    // DATA MEMBERS
    VelodyneReader  *reader_;
    int             idx_mn_, idx_mx_;

    // node based variables, which are shared between Solid and Shell
    int            *map_;       // node index mapping, from read-in index to dataset order
    vtkPoints      *crd_;       // node coordindates
    int             nnvs_;      // number of node based variables
    std::vector<std::string>    nvname_;
    std::vector<vtkFloatArray*> nvdata_;
    vtkObjectBase  *pobj_;

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);

  protected:
    static const std::string node_name;
    static const std::string solid_name;
    static const std::string shell_name;
    static const std::string surface_name;
    static const std::string particle_name;
    static const std::string tiednode_name;
    static const std::string sph_name;
    static const std::string invalid_name;

    static std::string composeName( const std::string& m, const std::string& v, const char app='/' );
    static void decomposeName( const std::string& s, std::string& m, std::string& v );
    static int getTypeId( const char* name );
    static const std::string& getTypeName( int type );

};


#endif
