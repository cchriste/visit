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
//                        avtLocateAndPickZoneQuery.h                        //
// ************************************************************************* //

#ifndef AVT_LOCATEANDPICKZONE_QUERY_H
#define AVT_LOCATEANDPICKZONE_QUERY_H
#include <query_exports.h>

#include <avtDatasetQuery.h>
#include <PickAttributes.h>

class avtLocateQuery;
class avtPickQuery;
class vtkDataSet;

// ****************************************************************************
//  Class: avtLocateAndPickZoneQuery
//
//  Purpose:
//
//  Programmer: Kathleen Bonnell 
//  Creation:   October 22, 2007 
//
//  Modifications:
//    Kathleen Bonnell, Tue Mar  1 13:02:19 PST 2011
//    Added SetNumVars.
//
//    Kathleen Biagas, Mon Jun 20 10:37:22 PDT 2011
//    Added SetInputParams, remove SetNumVars, added domain, zone.
//
// ****************************************************************************

class QUERY_API avtLocateAndPickZoneQuery : public avtDatasetQuery
{
  public:
                              avtLocateAndPickZoneQuery();
    virtual                  ~avtLocateAndPickZoneQuery();


    virtual const char       *GetType(void)
                                { return "avtLocateAndPickZoneQuery"; }
    virtual const char       *GetDescription(void)
                                { return "Locating and Picking Zone."; }
    virtual const char       *GetShortDescription(void)
                                { return "Pick Zone"; }

    virtual void              SetInputParams(const MapNode &); 

    virtual int               GetNFilters(void) { return 2; }

    virtual void              PerformQuery(QueryAttributes *);

    void                      SetPickAttsForTimeQuery(const PickAttributes *pa);

    virtual void              SetInvTransform(const avtMatrix *m);
    virtual void              SetNeedTransform(const bool v);

  protected:
    virtual void              Execute(vtkDataSet*, const int){;}

    PickAttributes            pickAtts;

  private:
    void                      SetPickAtts(const PickAttributes *);
    avtLocateQuery           *lcq;
    avtPickQuery             *zpq;
    int                       domain;
    int                       zone;
};


#endif


