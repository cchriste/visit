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
//                           avtCurve2DFileFormat.h                          //
// ************************************************************************* //

#ifndef AVT_CURVE2D_FILE_FORMAT_H
#define AVT_CURVE2D_FILE_FORMAT_H

#include <avtSTSDFileFormat.h>

#include <vector>
#include <string>
#include <visitstream.h>


class     vtkRectilinearGrid;


// ****************************************************************************
//  Class: avtCurve2DFileFormat
//
//  Purpose:
//      A file format reader for curves.
//
//  Programmer: Hank Childs
//  Creation:   May 28, 2002
//
//  Modifications:
//
//    Hank Childs, Fri Aug  1 21:16:55 PDT 2003
//    Made the format be a STSD.
//
//    Kathleen Bonnell, Fri Oct 28 13:02:51 PDT 2005 
//    Added methods GetTime, GetCycle, and members curveTime, curveCycle.
//
//    Kathleen Bonnell, Mon Jul 31 10:15:00 PDT 2006 
//    Represent curve as 1D RectilinearGrid instead of PolyData. 
//
//    Kathleen Bonnell, Thu Aug  3 08:42:33 PDT 2006 
//    Added dataExtents. 
//
//    Mark C. Miller, Tue Oct 31 20:33:29 PST 2006
//    Added VALID_XVALUE token to support "zone-centered" curves
//
//    Kathleen Bonnell, Tue Jan 20 11:02:57 PST 2009
//    Added spatialExtents. 
//
//    Kathleen Biagas, Tue Jul 15 14:16:07 MST 2014
//    Change 'GetPoint' args from float to double.
//
// ****************************************************************************

typedef enum
{
    VALID_POINT       = 0,
    HEADER,          /* 1 */
    WHITESPACE,      /* 2 */
    INVALID_POINT,   /* 3 */
    VALID_XVALUE
} CurveToken;


class avtCurve2DFileFormat : public avtSTSDFileFormat
{
  public:
                          avtCurve2DFileFormat(const char *);
    virtual              ~avtCurve2DFileFormat();
    
    virtual const char   *GetType(void) { return "Curve File Format"; };

    virtual double        GetTime(void);
    virtual int           GetCycle(void);
    
    virtual vtkDataSet   *GetMesh(const char *);
    virtual vtkDataArray *GetVar(const char *);

    virtual void          PopulateDatabaseMetaData(avtDatabaseMetaData *);

  protected:
    std::string           filename;
    bool                  readFile;

    std::vector<vtkRectilinearGrid *> curves;
    std::vector<std::string>   curveNames;
    std::vector<double>        spatialExtents;
    std::vector<double>        dataExtents;
    double                     curveTime;
    int                        curveCycle;

    void                  ReadFile(void);
    CurveToken            GetPoint(ifstream &, double &, double &,
                                   std::string &);
};


#endif


