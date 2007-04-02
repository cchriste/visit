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
//                            avtSimV1FileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_SIMV1_FILE_FORMAT_H
#define AVT_SIMV1_FILE_FORMAT_H

#include <database_exports.h>

#include <avtSTMDFileFormat.h>
#include <VisItDataInterface_V1.h>
#include <avtSimulationInformation.h>
#include <avtMaterial.h>

#include <vector>
#include <set>
#include <string>

// ****************************************************************************
//  Class: avtSimV1FileFormat
//
//  Purpose:
//      Reads in a .sim1 file for VisIt in the MDServer.
//      Reads in simulation data as input to VisIt in the Engine.
//
//  Programmer: Jeremy Meredith
//  Creation:   March 10, 2005
//
//  Modifications:
//    Jeremy Meredith, Thu Apr 14 16:47:07 PDT 2005
//    Added Curve and Material support.
//
//    Jeremy Meredith, Wed May 11 11:02:34 PDT 2005
//    Added ghost zone support.  Added restricted load balancing support.
//
// ****************************************************************************

class avtSimV1FileFormat : public avtSTMDFileFormat
{
  public:
                       avtSimV1FileFormat(const char *);
    virtual           ~avtSimV1FileFormat() {;};

    virtual const char    *GetType(void)   { return "SimV1"; };
    virtual void           FreeUpResources(void); 
    virtual int            GetCycle() { return -1; }

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);
    virtual avtMaterial   *GetMaterial(int, const char *);
    virtual vtkDataSet    *GetCurve(const char *);

    virtual void          *GetAuxiliaryData(const char *var, int domain,
                                            const char *type, void *,
                                            DestructorFunction &df);

    virtual void           PopulateIOInformation(avtIOInformation& ioInfo);

  protected:
    avtSimulationInformation simInfo;
    VisIt_SimulationCallback cb;
    std::set<std::string>    curveMeshes;

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *);
};


#endif
