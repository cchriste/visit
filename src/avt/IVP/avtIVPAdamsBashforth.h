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
//                               avtIVPAdamsBashforth.h                      //
// ************************************************************************* //

#ifndef AVT_IVPADAMSBASHFORTH_H
#define AVT_IVPADAMSBASHFORTH_H

#include <avtIVPSolver.h>
#include <ivp_exports.h>

#define ADAMS_BASHFORTH_NSTEPS 5

// ****************************************************************************
//  Class: avtIVPAdamsBashforth
//
//  Purpose:
//      An implementation of avtIVPSolver which models the 5th-order 
//      Adams-Bashforth multi-step method.
//
//  Programmer: Dave Pugmire
//  Creation:   August 5, 2008
//
//  Modifications:
//    Dave Pugmire, Fri Aug  8 16:05:34 EDT 2008
//    Improved version of A-B solver that builds function history from
//    initial HK4 steps.
//
//    Dave Pugmire, Tue Aug 19, 17:38:03 EDT 2008
//    Changed how distanced based termination is computed.
//
//    Dave Pugmire, Wed Aug 20, 12:54:44 EDT 2008
//    Add a tolerance and counter for handling stiffness detection.
//
//    Dave Pugmire, Mon Feb 23, 09:11:34 EST 2009
//    Reworked the termination code. Added a type enum and value. Made num steps
//    a termination criterion.
//
//    Dave Pugmire, Tue Feb 24 10:49:33 EST 2009
//    Replaced Euler step with RK4 step. Removed the Moulton corrector.
//
//    Dave Pugmire, Mon Mar  9 15:35:05 EDT 2009
//    Fix serialization for parallel integration.
//
//    Dave Pugmire, Tue Dec  1 11:50:18 EST 2009
//    Switch from avtVec to avtVector.
//
// ****************************************************************************

class IVP_API avtIVPAdamsBashforth: public avtIVPSolver
{
  public:
    avtIVPAdamsBashforth();
    ~avtIVPAdamsBashforth();

    // begin a new IVP solution
    virtual void     Reset( const double& t_start,
                            const avtVector &y_start,
                            const avtVector& v_start = avtVector(0,0,0) );

    // perform a single integration step
    // adaptive stepsize control retries until success or underflow
    virtual Result Step(avtIVPField* field, double t_max,
                        avtIVPStep* ivpstep = NULL);

    virtual void   OnExitDomain();

    virtual void   SetTolerances(const double& reltol, const double& abstol);

    virtual avtIVPAdamsBashforth* Clone() const
    {
        return new avtIVPAdamsBashforth( *this );
    }

  protected:
    // state serialization
    virtual void     AcceptStateVisitor(avtIVPStateHelper &aiss);
    
  private:
    int abStep, abNSteps;
    int numStep;
    int degenerate_iterations;
    double stiffness_eps;
    avtVector history[ADAMS_BASHFORTH_NSTEPS];
//    avtVector dhistory[ADAMS_BASHFORTH_NSTEPS];
};
#endif
