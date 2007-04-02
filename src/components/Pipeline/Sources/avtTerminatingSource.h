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
//                           avtTerminatingSource.h                          //
// ************************************************************************* //

#ifndef AVT_TERMINATING_SOURCE_H
#define AVT_TERMINATING_SOURCE_H

#include <pipeline_exports.h>

#include <void_ref_ptr.h>

#include <avtQueryableSource.h>
#include <avtDataSpecification.h>
#include <avtPipelineSpecification.h>


class     avtInlinePipelineSource;
class     avtMetaData;


typedef avtDataSpecification_p (*LoadBalanceFunction)(void *,
                                                   avtPipelineSpecification_p);
typedef bool                   (*DynamicCheckFunction)(void *,
                                                   avtPipelineSpecification_p);
typedef void                   (*InitializeProgressCallback)(void *, int);


// ****************************************************************************
//  Class: avtTerminatingSource
//
//  Purpose:
//      This is the terminator of a pipeline's update/execute cycle.  As such,
//      it owns the actual pipeline object.
//
//  Programmer: Hank Childs
//  Creation:   May 23, 2001
//
//  Modifications:
//
//    Jeremy Meredith, Thu Sep 20 01:01:49 PDT 2001
//    Added dynamic checker callback function so it can make determinations
//    about the network it will execute.
//
//    Hank Childs, Thu Oct 25 16:44:24 PDT 2001
//    Add a virtual function for whether or not this source permits dynamic
//    load balancing.
//
//    Kathleen Bonnell, Fri Nov 15 09:07:36 PST 2002
//    Added virtual method Query.
//
//    Hank Childs, Mon Jul 28 16:33:34 PDT 2003
//    Derived from avtQueryableSource instead of avtDataObjectSource.
//
//    Kathleen Bonnell, Wed Nov 12 18:26:21 PST 2003 
//    Added virtual method FindElementForPoint. 
//
//    Kathleen Bonnell, Mon Dec 22 14:48:57 PST 2003 
//    Added virtual method GetDomainName. 
//
//    Kathleen Bonnell, Tue May 25 16:16:25 PDT 2004 
//    Add virtual method 'QueryZoneCenter'.
//
//    Jeremy Meredith, Wed Jun  9 09:13:39 PDT 2004
//    Added species aux data.
//
//    Kathleen Bonnell, Thu Jun 10 18:31:22 PDT 2004 
//    Renamed QueryZoneCenter to QueryCoords, added bool arg.
//
//    Kathleen Bonnell, Thu Dec 16 17:16:33 PST 2004 
//    Added another bool arg to QueryCoords. 
//
//    Kathleen Bonnell, Mon Jan  3 13:40:42 PST 2005 
//    Added GetSIL method. 
//
//    Kathleen Bonnell, Tue Jan 25 07:59:28 PST 2005 
//    Added const char *arg to QueryCoords. 
//
// ****************************************************************************

class PIPELINE_API avtTerminatingSource : virtual public avtQueryableSource
{
    friend class                   avtInlinePipelineSource;

  public:
                                   avtTerminatingSource();
    virtual                       ~avtTerminatingSource();

    virtual avtTerminatingSource  *GetTerminatingSource(void)  { return this;};

    avtMetaData                   *GetMetaData(void) { return metadata; };
    void                           GetMeshAuxiliaryData(const char *type,
                                       void *args, avtPipelineSpecification_p,
                                       VoidRefList &);
    void                           GetVariableAuxiliaryData(const char *type,
                                       void *args, avtPipelineSpecification_p,
                                       VoidRefList &);
    void                           GetMaterialAuxiliaryData(const char *type,
                                       void *args, avtPipelineSpecification_p,
                                       VoidRefList &);
    void                           GetSpeciesAuxiliaryData(const char *type,
                                       void *args, avtPipelineSpecification_p,
                                       VoidRefList &);

    virtual bool                   Update(avtPipelineSpecification_p);

    static void                    SetLoadBalancer(LoadBalanceFunction,void *);
    static void                    SetDynamicChecker(DynamicCheckFunction,
                                                     void *);
    static void                    RegisterInitializeProgressCallback(
                                           InitializeProgressCallback, void *);

    virtual avtDataSpecification_p GetFullDataSpecification(void);
    avtPipelineSpecification_p     GetGeneralPipelineSpecification(void);

    // Define this so derived types don't have to.
    virtual void                   Query(PickAttributes *){;};
    virtual avtSIL                *GetSIL(int timestep){return NULL;};
    virtual bool                   FindElementForPoint(const char*, const int, 
                                       const int, const char*, double[3], int &)
                                       { return false;};
    virtual void                   GetDomainName(const std::string &, const int, 
                                       const int, std::string &) {;};
    virtual bool                   QueryCoords(const std::string &, const int, 
                                       const int, const int, double[3],
                                       const bool, const bool = false, const char *mn=NULL)
                                       { return false;};

  protected:
    avtMetaData                   *metadata;
    static InitializeProgressCallback
                                   initializeProgressCallback;
    static void                   *initializeProgressCallbackArgs;

    virtual bool                   FetchData(avtDataSpecification_p) = 0;
    virtual void                   FetchMeshAuxiliaryData(const char *type,
                                       void *args, avtDataSpecification_p,
                                       VoidRefList &);
    virtual void                   FetchVariableAuxiliaryData(const char *type,
                                       void *args, avtDataSpecification_p,
                                       VoidRefList &);
    virtual void                   FetchMaterialAuxiliaryData(const char *type,
                                       void *args, avtDataSpecification_p,
                                       VoidRefList &);
    virtual void                   FetchSpeciesAuxiliaryData(const char *type,
                                       void *args, avtDataSpecification_p,
                                       VoidRefList &);

    avtDataSpecification_p         BalanceLoad(avtPipelineSpecification_p);

    virtual bool                   UseLoadBalancer(void) { return true; };
    virtual bool                   ArtificialPipeline(void) { return false; };
    virtual int                    NumStagesForFetch(avtDataSpecification_p);

 private:
    static LoadBalanceFunction     loadBalanceFunction;
    static void                   *loadBalanceFunctionArgs;
    static DynamicCheckFunction    dynamicCheckFunction;
    static void                   *dynamicCheckFunctionArgs;

    void                           InitPipeline(avtPipelineSpecification_p);
    virtual bool                   CanDoDynamicLoadBalancing(void);

};


#endif


