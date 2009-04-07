/*****************************************************************************
*
* Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-400142
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
//                         avtIVPVTKTimeVaryingField.h                       //
// ************************************************************************* //

#ifndef AVT_IVPVTKTimeVaryingFIELD_H
#define AVT_IVPVTKTimeVaryingFIELD_H

#include <avtIVPField.h>

#include <vtkVisItInterpolatedVelocityField.h>

#include <ivp_exports.h>


// ****************************************************************************
//  Class:  avtIVPVTKTimeVaryingField
//
//  Purpose:
//    A wrapper class to allow the use of vtkDataSets as IVP fields for 
//    streamline integration. Uses vtkVisItInterpolatedVelocityField on top of 
//    the supplied vtkDataSet. 
//
//  Programmer:  Dave Pugmire (on behalf of Hank Childs)
//  Creation:    Tue Feb 24 09:24:49 EST 2009
//
//  Modifications:
//
//   Dave Pugmire, Tue Mar 10 12:41:11 EDT 2009
//   Add GetValidTimeRange.
//
//   Hank Childs, Thu Apr  2 16:40:08 PDT 2009
//   Use vtkVisItInterpolatedVelocityField.
//
//   Hank Childs, Tue Apr  7 08:52:59 CDT 2009
//   Use a single vtkVisItInterpolatedVelocityField, which saves on
//   computation.
//
// ****************************************************************************

class IVP_API avtIVPVTKTimeVaryingField: public avtIVPField
{
  public:
                   avtIVPVTKTimeVaryingField(vtkVisItInterpolatedVelocityField *,
                                             double t1, double t2);
                   ~avtIVPVTKTimeVaryingField();

    // avtIVPField interface
    avtVec         operator()(const double& t, const avtVecRef& x) const;
    double         ComputeVorticity(const double& t, const avtVecRef& x) const;
    
    bool           IsInside( const double& t, const avtVecRef& x ) const;
    unsigned int   GetDimension() const;
    void           SetNormalized( bool v );
    virtual bool   GetValidTimeRange(double range[]) const;

  protected:

    double                               time1, time2;
    vtkVisItInterpolatedVelocityField   *iv1;
    bool           normalized;
    
};

#endif


