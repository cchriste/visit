/*****************************************************************************
*
* Copyright (c) 2000 - 2012, Lawrence Livermore National Security, LLC
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
//                             avtIntegralCurve.h                            //
// ************************************************************************* //

#ifndef AVT_INTEGRAL_CURVE_H
#define AVT_INTEGRAL_CURVE_H

#include <avtIVPSolver.h>
#include <MemStream.h>
#include <string>
#include <vector>
#include <avtVector.h>

class vtkObject;

typedef bool (*avtIntegralCurveCallback)(void);

// ****************************************************************************
// Class: BlockIDType
//
// Purpose:
//    Encapsulate the a domain/timestep.
//    
//
// Programmer: Dave Pugmire
// Creation:   Tue Mar 10 12:41:11 EDT 2009
//
// Modifications:
//
//   Dave Pugmire, Mon May 11 12:41:51 EDT 2009
//   Fix operator< so that that std::map works.
//
// ****************************************************************************

class IVP_API BlockIDType
{
  public:
    BlockIDType() :domain(-1), timeStep(0) {}
    BlockIDType(const int &d) :domain(d), timeStep(0) {}
    BlockIDType(const int &d, const int &t) :domain(d), timeStep(t) {}
    ~BlockIDType() {}

    void operator=(const BlockIDType &dt)
    {
        domain=dt.domain;
        timeStep=dt.timeStep;
    }

    bool operator==(const BlockIDType &dt) const
    {
        return (domain == dt.domain &&
                timeStep == dt.timeStep);
    }
    bool operator<(const BlockIDType &dt) const
    {
        return (domain < dt.domain) ||
               ((domain == dt.domain) && timeStep < dt.timeStep);
    }

    //Members
    int domain, timeStep;

    friend std::ostream& operator<<(std::ostream &out, const BlockIDType &d)
    {
        out<<"["<<d.domain<<", "<<d.timeStep<<"]";
        return out;
    }
};


// ****************************************************************************
//  Class: avtIntegralCurve
//
//  Purpose:
//      avtIntegralCurve is a straightforward implementation of integral curves,
//      based on avtIVPSolver.  Through this model instance, a user of 
//      avtIntegralCurve is able to select any IVP scheme to be used in the 
//      integration.
//
//  Programmer: Christoph Garth
//  Creation:   February 25, 2008
//
//  Modifications:
//
//    Dave Pugmire, Wed Aug 13 10:58:32 EDT 2008
//    Modify how data without ghost zones are handled.
//
//    Dave Pugmire, Tue Aug 19, 17:38:03 EDT 2008
//    Chagned how distanced based termination is computed.
//
//    Dave Pugmire, Wed Dec  3 08:33:42 EST 2008
//    Added maxSteps argument to Advance() to optionally control how many
//    integration steps are taken.
//
//    Dave Pugmire, Mon Feb 23, 09:11:34 EST 2009
//    Reworked the termination code. Added a type enum and value. Made num steps
//    a termination criterion. Code cleanup: We no longer need fwd/bwd solvers.
//    Removed the plane intersection code.
//
//   Dave Pugmire, Mon Jun 8 2009, 11:44:01 EDT 2009
//   Removed the wantVorticity, extents and ghostzone flags. Extents and ghost
//   zones are handled by the vtkDataSet itself. The wantVorticity was replaced
//   with a scalarValueType which can be 'or'-d together to specify what to
//   compute.
//
//   Dave Pugmire, Tue Aug 11 10:25:45 EDT 2009
//   Add new termination criterion: Number of intersections with an object.
//
//   Dave Pugmire, Tue Aug 18 08:47:40 EDT 2009
//   Don't record intersection points, just count them.
//
//   Dave Pugmire, Thu Sep 24 13:52:59 EDT 2009
//   Option to serialize steps.
//
//   Dave Pugmire, Tue Dec  1 11:50:18 EST 2009
//   Switch from avtVec to avtVector.
//
//   Dave Pugmire, Tue Dec 29 14:37:53 EST 2009
//   Generalize the compute scalar variable.
//
//   Dave Pugmire, Tue Feb 23 09:42:25 EST 2010
//   Use a vector instead of a list for the integration steps.
//
//   Dave Pugmire, Wed May 26 13:48:24 EDT 2010
//   New return enum.
//
//   Hank Childs, Thu Jun  3 10:44:46 PDT 2010
//   Remove TMin, PtStart, TStart, IsForward, IsBackward, and IsBothDir.
//   Rename TMax to GetCurrentTime, PtEnd to GetCurrentLocation.
//
//   Hank Childs, Fri Jun  4 15:45:39 CDT 2010
//   Combine this class with the contents of avtStreamlineWrapper.
//
//   Hank Childs, Fri Jun  4 21:30:18 CDT 2010
//   Separate out portions specific to Poincare and Streamline into
//   avtStateRecorderIntegralCurve.
//
//   Hank Childs, Tue Jun  8 09:30:45 CDT 2010
//   Put sequence tracking code into avtStateRecorderIntegralCurve.
//
//   Hank Childs, Mon Oct  4 15:03:43 PDT 2010
//   Remove termination code.  It now goes in derived types.
//
//   Dave Pugmire, Fri Nov  5 15:34:49 EDT 2010
//   Add counter to handle communication of ICs
//
//   Hank Childs, Sun Dec  5 11:43:46 PST 2010
//   Added data member for tracking when we encounter numerical problems.
//
//   Dave Pugmire, Fri Feb 18 14:52:18 EST 2011
//   Replaced minH with minHFactor for use when integrating upto a domain boundary.
//
//   Hank Childs, Tue Dec  6 16:23:47 PST 2011
//   Add virtual methods LessThan (for sorting) and 
//   PrepareForFinalCommunication.
//
//   David Camp, Wed Mar  7 10:43:07 PST 2012
//   Added a Serialize flag to the arguments. This is to support the restore
//   ICs code.
//
// ****************************************************************************

class IVP_API avtIntegralCurve
{
  public:

    enum Direction
    {
        DIRECTION_FORWARD  = 0,
        DIRECTION_BACKWARD = 1,
    };

    enum Status
    {
        STATUS_OK         = 0,
        STATUS_FINISHED   = 1,
        STATUS_TERMINATED = 2,
    };

    enum SerializeFlags
    {
        SERIALIZE_ALL     = -1,
        SERIALIZE_NO_OPT  = 0,
        SERIALIZE_STEPS   = 1,
        SERIALIZE_INC_SEQ = 2,
    };

    avtIntegralCurve( const avtIVPSolver* model, 
                      Direction dir, 
                      const double& t_start, 
                      const avtVector &p_start, 
                      const avtVector &v_start, 
                      long ID );

    avtIntegralCurve();
    virtual ~avtIntegralCurve();

    void Advance(avtIVPField* field);

    double    CurrentTime() const;
    void      CurrentLocation(avtVector &end);

    virtual void      Serialize(MemStream::Mode mode, MemStream &buff, 
                                avtIVPSolver *solver, SerializeFlags serializeFlags);

    virtual void      PrepareForSend(void) { ; };
    virtual void      ResetAfterSend(void) { ; };

    virtual bool      SameCurve(avtIntegralCurve *ic)
                               { return id == ic->id; };

    static bool       DomainCompare(const avtIntegralCurve *slA,
                                    const avtIntegralCurve *slB);

    bool              EncounteredNumericalProblems(void)
                             { return encounteredNumericalProblems; };

    void     SetPostStepCallback(avtIntegralCurveCallback func) {postStepCallbackFunction = func; }

    virtual avtIntegralCurve* MergeIntegralCurveSequence(std::vector<avtIntegralCurve *> &v) = 0;
    virtual void      PrepareForFinalCommunication(void) {;};

    // This is used for sorting, particularly for parallel communication
    virtual bool LessThan(const avtIntegralCurve *ic) const;

  protected:
    avtIntegralCurveCallback postStepCallbackFunction;

    avtIntegralCurve( const avtIntegralCurve& );
    avtIntegralCurve& operator=( const avtIntegralCurve& );
    
    virtual void AnalyzeStep( avtIVPStep& step,
                              avtIVPField* field ) = 0;
    virtual bool    UseFixedTerminationTime(void) { return false; };
    virtual double  FixedTerminationTime(void)    { return 0; };

  public:

    Status    status;
    Direction direction;

    // Helpers needed for figuring out which domain to use next
    std::vector<BlockIDType> seedPtDomainList;
    BlockIDType domain;
    long long sortKey;

    long id;
    int counter;

    bool   encounteredNumericalProblems;
    int    originatingRank;
  protected:

    avtIVPSolver*       ivp;

    static const double minHFactor;
};


// ostream operators for avtIntegralCurve's enum types
inline std::ostream& operator<<( std::ostream& out, 
                                 avtIntegralCurve::Status status )
{
    switch( status )
    {
    case avtIntegralCurve::STATUS_OK:
        return out << "OK";
    case avtIntegralCurve::STATUS_FINISHED:
        return out << "FINISHED";
    case avtIntegralCurve::STATUS_TERMINATED:
        return out << "TERMINATED";
    default:
        return out << "UNKNOWN";
    }
}

inline std::ostream& operator<<( std::ostream& out, 
                                 avtIntegralCurve::Direction dir )
{
    switch( dir )
    {
    case avtIntegralCurve::DIRECTION_FORWARD: 
        return out << "FORWARD";
    case avtIntegralCurve::DIRECTION_BACKWARD:
        return out << "BACKWARD";
    default:
        return out << "UNKNOWN";
    }
}


#endif // AVT_STREAMLINE_H


