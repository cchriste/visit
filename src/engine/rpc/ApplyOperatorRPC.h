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

#ifndef APPLY_OPERATOR_RPC_H
#define APPLY_OPERATOR_RPC_H
#include <engine_rpc_exports.h>

#include <VisItRPC.h>
#include <string>

class ApplyOperatorRPC;

// ****************************************************************************
//  Class:  PrepareOperatorRPC
//
//  Purpose:
//    Signals the name of the operator about to be created so that
//    the ApplyOperatorRPC has space to store the correct attributes.
//
//  Programmer:  Jeremy Meredith
//  Creation:    March  2, 2001
//
// ****************************************************************************

class ENGINE_RPC_API PrepareOperatorRPC : public BlockingRPC
{
  public:
    PrepareOperatorRPC();
    ~PrepareOperatorRPC();

    const std::string TypeName() const { return "PrepareOperatorRPC";};

    void SetApplyOperatorRPC(ApplyOperatorRPC*);
    ApplyOperatorRPC *GetApplyOperatorRPC();

    void operator()(const std::string &n);
    void SelectAll();
    std::string GetID();
  private:
    ApplyOperatorRPC *applyOperatorRPC;
    std::string id;
};


// ****************************************************************************
//  Class:  ApplyOperatorRPC
//
//  Purpose:
//    Apply an operator.
//
//  Programmer:  Jeremy Meredith
//  Creation:    March  2, 2001
//
//  Modifications:
//    Brad Whitlock, Mon Mar 25 09:56:52 PDT 2002
//    Removed SetSocket.
//
// ****************************************************************************

class ENGINE_RPC_API ApplyOperatorRPC : public BlockingRPC
{
  public:
    ApplyOperatorRPC();
    virtual ~ApplyOperatorRPC();

    const std::string TypeName() const { return "ApplyOperatorRPC";};

    void operator()(const std::string&, const AttributeSubject*);

    virtual void SelectAll();

    std::string GetID();
    AttributeSubject *GetAtts();
    PrepareOperatorRPC &GetPrepareOperatorRPC();

    void SetAtts(AttributeSubject*);

    virtual void SetXfer(Xfer *x);

  private:
    AttributeSubject *atts;
    PrepareOperatorRPC prepareOperatorRPC;
};

#endif
