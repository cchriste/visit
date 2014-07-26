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

#ifndef LAUNCHER_PROXY_H
#define LAUNCHER_PROXY_H
#include <visit-config.h>
#include <vclproxy_exports.h>
#include <RemoteProxyBase.h>
#include <LauncherState.h>
#include <LauncherMethods.h>
#include <map>

// ****************************************************************************
// Class: LauncherProxy
//
// Purpose:
//   This is a proxy class for the launcher program.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri May 2 16:11:43 PST 2003
//
// Modifications:
//    Jeremy Meredith, Tue Mar 30 10:08:18 PST 2004
//    I added support for simulations.
//
//    Jeremy Meredith, Wed May 11 09:04:52 PDT 2005
//    Added security key to simulation connection.
//
//    Jeremy Meredith, Thu May 24 10:21:48 EDT 2007
//    Added method to retrieve the SSH tunneling local-to-remote port map.
//
//    Thomas R. Treadway, Mon Oct  8 13:27:42 PDT 2007
//    Backing out SSH tunneling on Panther (MacOS X 10.3)
//
// ****************************************************************************

class LAUNCHER_PROXY_API LauncherProxy : public RemoteProxyBase
{
public:
    LauncherProxy();
    virtual ~LauncherProxy();

    virtual std::string GetComponentName() const;

    virtual void Create(const MachineProfile &profile, 
                        ConnectCallback *connectCallback = 0, 
                        void *connectCallbackData = 0,
                        bool createAsThoughLocal = false);

    std::map<int,int> GetPortTunnelMap();

    LauncherMethods *GetLauncherMethods() { return methods; }
    LauncherState   *GetLauncherState() { return state; }

protected:
    virtual void SetupComponentRPCs();
private:
    LauncherMethods *methods;
    LauncherState   *state;
};

#endif
