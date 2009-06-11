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
//                                avtIVPField.h                              //
// ************************************************************************* //

#ifndef AVT_IVPFIELD_H
#define AVT_IVPFIELD_H

#include <stdexcept>
#include <avtVec.h>
#include <ivp_exports.h>


// ****************************************************************************
//  Class: avtIVPField
//
//  Purpose:
//      avtIVPField is a base class for all manners of vector fields. 
//      Deriving from it should allow an adaptation of many different vector 
//      field types for the use of streamlines/IVP solutions by wrapping 
//      existing interpolation facilities. If the queried point is invalid, 
//      avtIVPField can throw an exception that will pass through avtIVPSolver.
//
//      The IVP right-hand side is made accessible to an IVP solver by means of //      the avtIVPField class, allowing the IVP solver to query points of the 
//      given vector field. 
//
//  Programmer: Christoph Garth
//  Creation:   February 25, 2008
//
//  Modifications:
//
//   Dave Pugmire, Tue Mar 10 12:41:11 EDT 2009
//   Add GetValidTimeRange.
//
//   Dave Pugmire, Mon Jun 8 2009, 11:44:01 EDT 2009
//   Added ComputeScalarVariable, HasGhostZones and GetExtents methods.
//
// ****************************************************************************

class IVP_API avtIVPField
{
  public:
    class Undefined: public std::exception
    {
      public:
        const char* what() const throw()
        {
            return "field undefined";
        }
    };

                         avtIVPField() {};

    virtual avtVec       operator()(const double& t, 
                                    const avtVecRef& x) const = 0;
    virtual double       ComputeVorticity(const double& t, 
                                          const avtVecRef& x ) const = 0;
    virtual double       ComputeScalarVariable(const double& t,
                                               const avtVecRef& x) const = 0;

    virtual bool         IsInside(const double& t, 
                                  const avtVecRef& x) const = 0;
    virtual unsigned int GetDimension() const = 0;
    virtual bool         GetValidTimeRange(double range[]) const = 0;
    virtual bool         HasGhostZones() const = 0;
    virtual void         GetExtents(double *extents) const = 0;
};

#endif


