/*****************************************************************************
*
* Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
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
//                        avtVariableByZoneQuery.h                           //
// ************************************************************************* //

#ifndef AVT_VARIABLEBYZONE_QUERY_H
#define AVT_VARIABLEBYZONE_QUERY_H
#include <query_exports.h>

#include <avtPickByZoneQuery.h>

#include <PickAttributes.h>


class vtkDataSet;


// ****************************************************************************
//  Class: avtVariableByZoneQuery
//
//  Purpose:
//    A query that retrieves var information about a mesh given a 
//    particular domain and zone number.
//
//  Programmer: Kathleen Bonnell
//  Creation:   July 29, 2004
//
//  Modifications:
//    Kathleen Bonnell, Tue Nov  8 10:45:43 PST 2005
//    Added avtDataAttributes arg to Preparation.
//
//    Kathleen Bonnell, Tue Jul  8 15:40:45 PDT 2008
//    Added GetTimeCurveSpecs method.
//
//    Kathleen Bonnell, Tue Mar  1 16:07:47 PST 2011
//    Added SetNumVars method.
//
//    Kathleen Biagas, Mon Jun 20 10:31:43 PDT 2011
//    Added SetInputParams, removed SetNumVars, added domain, zone.
//
// ****************************************************************************

class QUERY_API avtVariableByZoneQuery : public avtPickByZoneQuery
{
  public:
                              avtVariableByZoneQuery();
    virtual                  ~avtVariableByZoneQuery();


    virtual const char       *GetType(void)   
                              { return "avtVariableByZoneQuery"; };
    virtual const char       *GetDescription(void)
                              { return "Retrieving var information on mesh."; }

    virtual void              SetInputParams(const MapNode &);

    virtual const MapNode    &GetTimeCurveSpecs(); 

  protected:
    virtual void                    Preparation(const avtDataAttributes &);
    virtual void                    PostExecute(void);

    int domain;
    int zone;
};


#endif
