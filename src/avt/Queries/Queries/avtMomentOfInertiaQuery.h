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
//                        avtMomentOfInertiaQuery.h                          //
// ************************************************************************* //

#ifndef AVT_MOMENT_OF_INERTIA_QUERY_H
#define AVT_MOMENT_OF_INERTIA_QUERY_H

#include <avtDatasetQuery.h>

#include <query_exports.h>

class vtkDataSet;
class vtkCell;

class avtAbsValExpression;
class avtBinaryMultiplyExpression;
class avtVMetricVolume;


// ****************************************************************************
//  Class: avtMomentOfInertiaQuery
//
//  Purpose:
//      This query calculates the moment of inertia tensor for 3D datasets.
//
//  Programmer: Hank Childs
//  Creation:   May 17, 2005
//
//  Modifications:
//
//    Hank Childs, Fri Aug 12 15:30:59 PDT 2005
//    Add absval.
//
// ****************************************************************************

class QUERY_API avtMomentOfInertiaQuery : public avtDatasetQuery
{
  public:
                                    avtMomentOfInertiaQuery();
    virtual                        ~avtMomentOfInertiaQuery();

    virtual const char             *GetType(void)
                                     {return "avtMomentOfInertiaQuery";};
    virtual const char             *GetDescription(void)
                                     {return "Calculating Moment of Inertia";};

  protected:
    double                          I[9];
    avtVMetricVolume               *volume;
    avtBinaryMultiplyExpression    *multiply;
    avtAbsValExpression            *absval;

    virtual void                    Execute(vtkDataSet *, const int);
    virtual void                    PreExecute(void);
    virtual void                    PostExecute(void);
    virtual void                    VerifyInput(void);
    virtual avtDataObject_p         ApplyFilters(avtDataObject_p);
};


#endif


