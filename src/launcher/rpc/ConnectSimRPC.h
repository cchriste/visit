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

#ifndef CONNECT_SIM_RPC_H
#define CONNECT_SIM_RPC_H
#include <vclrpc_exports.h>
#include <VisItRPC.h>
#include <vectortypes.h>

// ****************************************************************************
//  Class: ConnectSimRPC
//
//  Purpose:
//    This class encodes an RPC that tells the launcher to tell a simulation
//    to connect back to the viewer.
//
//  Notes:      
//
//  Programmer: Jeremy Meredith
//  Creation:   January 23, 2004
//
//  Modifications:
//    Jeremy Meredith, Mon May  9 14:41:35 PDT 2005
//    Added security key.
//
// ****************************************************************************

class LAUNCHER_RPC_API ConnectSimRPC : public BlockingRPC
{
  public:
    ConnectSimRPC();
    virtual ~ConnectSimRPC();

    // Invokation method
    void operator()(const stringVector &args, const std::string&, int,
                    const std::string&);

    // Property selection methods
    virtual void SelectAll();

    // Methods to access private data.
    const stringVector &GetLaunchArgs() const;
    const std::string  &GetSimHost() const;
    int                 GetSimPort() const;
    const std::string  &GetSimSecurityKey() const;

  private:
    stringVector launchArgs;
    std::string  simHost;
    int          simPort;
    std::string  simSecurityKey;
};


#endif
