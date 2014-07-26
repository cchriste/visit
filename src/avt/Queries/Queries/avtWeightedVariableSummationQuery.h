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
//                      avtWeightedVariableSummationQuery.h                  //
// ************************************************************************* //

#ifndef AVT_WEIGHTED_VARIABLE_SUMMATION_QUERY_H
#define AVT_WEIGHTED_VARIABLE_SUMMATION_QUERY_H

#include <query_exports.h>

#include <avtSummationQuery.h>

#include <string>

class     avtBinaryMultiplyExpression;
class     avtEdgeLength;
class     avtRevolvedVolume;
class     avtVMetricArea;
class     avtVMetricVolume;


// ****************************************************************************
//  Class: avtWeightedVariableSummationQuery
//
//  Purpose:
//      A query that will sum all of one variables values.
//
//  Programmer: Hank Childs
//  Creation:   February 3, 2004
//
//  Modifications:
//    Kathleen Bonnell, Wed Jul 28 08:50:51 PDT 2004
//    Added VerifyInput.
//
//    Kathleen Bonnell, Fri Feb  3 10:32:12 PST 2006 
//    Added revolvedVolume. 
//
//    Hank Childs, Thu May 11 13:28:50 PDT 2006
//    Added new virtual methods so that new queries can inherit from this.
//
//    Hank Childs, Wed Apr 28 05:25:52 PDT 2010
//    Add support for 1D cross sections.
//
// ****************************************************************************

class QUERY_API avtWeightedVariableSummationQuery : public avtSummationQuery
{
  public:
                         avtWeightedVariableSummationQuery();
    virtual             ~avtWeightedVariableSummationQuery();

    virtual const char  *GetType(void)  
                             { return "avtWeightedVariableSummationQuery"; };

  protected:
    avtEdgeLength               *length;
    avtVMetricArea              *area;
    avtVMetricVolume            *volume;
    avtRevolvedVolume           *revolvedVolume;
    avtBinaryMultiplyExpression *multiply;

    virtual avtDataObject_p    ApplyFilters(avtDataObject_p);
    virtual int                GetNFilters(void) { return 2; };
    virtual void               VerifyInput(void);

    virtual avtDataObject_p    CreateVariable(avtDataObject_p d) { return d; };
    virtual std::string        GetVarname(std::string &s) { return s; };
};


#endif


