/*****************************************************************************
*
* Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
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
//                            avtIDXFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_IDX_FILE_FORMAT_H
#define AVT_IDX_FILE_FORMAT_H

#include <avtMTMDFileFormat.h>

#include <vector>
#include <DBOptionsAttributes.h>

#include <visus.h>
#include <visuscpp/db/dataset/visus_db_dataset.h>


// ****************************************************************************
//  Class: avtIDXFileFormat
//
//  Purpose:
//      Reads in IDX files as a plugin to VisIt.
//
//  Programmer: spetruzza, camc -- generated by xml2avt
//  Creation:   Tue Jan 7 12:20:07 MST 2014
//
// ****************************************************************************

class DummyNode;
struct avtView3D;
class avtIDXFileFormat : public avtMTMDFileFormat, public Visus::Object
{
  public:
                       avtIDXFileFormat(const char *, DBOptionsAttributes* attrs);
    virtual           ~avtIDXFileFormat();

    //
    // This is used to return unconvention data -- ranging from material
    // information to information about block connectivity.
    //
    // virtual void      *GetAuxiliaryData(const char *var, int timestep, 
    //                                     int domain, const char *type, void *args, 
    //                                     DestructorFunction &);
    //

    //
    // If you know the times and cycle numbers, overload this function.
    // Otherwise, VisIt will make up some reasonable ones for you.
    //
    // virtual void        GetCycles(std::vector<int> &);
    // virtual void        GetTimes(std::vector<double> &);
    //

    virtual int            GetNTimesteps(void);

    virtual const char    *GetType(void)   { return "IDX"; };

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *, int);

    virtual vtkDataSet    *GetMesh(int, int, const char *);
    virtual vtkDataArray  *GetVar(int, int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, int, const char *);

//    virtual void GetCycles(std::vector<int> &);
    virtual void GetTimes(std::vector<double> &);
    
    virtual void           FreeUpResources(void); 

  protected:

    std::string             filename;
    int                     nprocs;
    int                     rank;
    int                     dim;         //2d or 3d
    static int              num_instances;

    Visus::UniquePtr<Visus::Application> app;
    Visus::UniquePtr<Visus::Access>      access;
    Visus::SharedPtr<Visus::Dataset>     dataset;
    std::vector<Visus::Box>              boxes;
    std::vector<int*>                    boxes_bounds;
    bool multibox;
  
  private:
    
    bool reverse_endian;
    
    void calculateBoundsAndExtents();
    void loadBalance();
    
    inline int
    int16_Reverse_Endian(short val, unsigned char *output)
    {
        unsigned char *input  = ((unsigned char *)&val);
        
        output[0] = input[1];
        output[1] = input[0];
        
        return 2;
    }
    
    inline int
    int32_Reverse_Endian(int val, unsigned char *outbuf)
    {
        unsigned char *data = ((unsigned char *)&val) + 3;
        unsigned char *out = outbuf;
        
        *out++ = *data--;
        *out++ = *data--;
        *out++ = *data--;
        *out = *data;
        
        return 4;
    }
    
    inline int
    float32_Reverse_Endian(float val, unsigned char *outbuf)
    {
        unsigned char *data = ((unsigned char *)&val) + 3;
        unsigned char *out = outbuf;
        
        *out++ = *data--;
        *out++ = *data--;
        *out++ = *data--;
        *out = *data;
        
        return 4;
    }
    
    inline int
    double64_Reverse_Endian(double val, unsigned char *outbuf)
    {
        unsigned char *data = ((unsigned char *)&val) + 7;
        unsigned char *out = outbuf;
        
        *out++ = *data--;
        *out++ = *data--;
        *out++ = *data--;
        *out++ = *data--;
        *out++ = *data--;
        *out++ = *data--;
        *out++ = *data--;
        *out = *data;
        
        return 8;
    }


    VISUS_DECLARE_BINDABLE(avtIDXFileFormat);
};


#endif
