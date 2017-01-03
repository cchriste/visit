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
//                             avtGyRadiusQuery.h                              //
// ************************************************************************* //

#ifndef AVT_GYRADIUS_QUERY_H
#define AVT_GYRADIUS_QUERY_H
#include <query_exports.h>
#include <avtDatasetQuery.h>
#include <avtExpressionEvaluatorFilter.h>

class vtkDataArray;
class vtkDataSet;
class avtCentroidQuery;

// ****************************************************************************
//  Class: avtGyRadiusQuery
//
//  Purpose:
//    A query that performs the following calculation: 
//    R = sqrt(sum( (Mi*(Ri - Rcm)^2) ) / sum(Mi)), where Mi is the mass of cell
//    i and (Ri - Rcm) is the distance of cell i from the center of mass.
//
//  Programmer: Kevin Griffin
//  Creation:   August 9, 2016
//
//  Modifications:
//
// ****************************************************************************

class QUERY_API avtGyRadiusQuery : virtual public avtDatasetQuery
{
public:
                            avtGyRadiusQuery();
    virtual                ~avtGyRadiusQuery();
    
    virtual const char     *GetType(void) { return "avtGyRadiusQuery"; };
    virtual const char     *GetDescription(void) { return "GyRadius."; };
    
    virtual void            SetInputParams(const MapNode &);
    static  void            GetDefaultInputParams(MapNode &);
    
protected:
    virtual void            PreExecute(void);
    virtual void            Execute(vtkDataSet *, const int);
    virtual void            PostExecute(void);
    virtual avtDataObject_p ApplyFilters(avtDataObject_p);
    
private:
    doubleVector    centroid;
    bool            overrideCentroid;
    double          totalSum;
    double          totalMass;
    bool            isZonal;
    std::string     varName;
    
    avtExpressionEvaluatorFilter *eef;
};

#endif
