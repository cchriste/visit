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

#ifndef LAUNCHER_APPLICATION_H
#define LAUNCHER_APPLICATION_H
#include <visit-config.h>
#include <ParentProcess.h>
#include <Xfer.h>
#include <QuitRPC.h>
#include <KeepAliveRPC.h>
#include <LaunchRPC.h>
#include <ConnectSimRPC.h>
#include <RPCExecutor.h>
#include <LaunchService.h>
#include <LauncherState.h>
#include <map>

// ****************************************************************************
// Class: LauncherApplication
//
// Purpose:
//   This class contains the launcher application that is responsible for
//   launching VisIt components.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Fri May 2 16:55:28 PST 2003
//
// Modifications:
//   Brad Whitlock, Fri Mar 12 10:35:29 PDT 2004
//   I added KeepAliveRPC.
//
//   Jeremy Meredith, Tue Mar 30 10:06:33 PST 2004
//   I added support for connecting to running simulations.
//
//   Jeremy Meredith, Mon May  9 16:09:06 PDT 2005
//   I added security protocols for simulations.
//
//   Jeremy Meredith, Thu May 24 12:31:29 EDT 2007
//   Added support for SSH tunneling.  Specifically, added a bridge
//   so that the login node could act as a gateway between the SSH 
//   tunnel origination (on the login node) and the compute node
//   where a parallel engine might run.
//
//   Thomas R. Treadway, Mon Oct  8 13:27:42 PDT 2007
//   Backout SSH tunneling changes for Panther (MacOS X 10.3)
//
//   Brad Whitlock, Wed Nov 21 10:34:14 PST 2007
//   Added support for forwarding child process output to the client.
//
//   Brad Whitlock, Mon Apr 27 16:25:34 PST 2009
//   I changed the code so we can tunnel simulation data connections.
//
//   Brad Whitlock, Tue Nov 29 11:23:39 PST 2011
//   I moved some functionality into LaunchService.
//
// ****************************************************************************

class LauncherApplication
{
public:
    static LauncherApplication *Instance();
    virtual ~LauncherApplication();
    void Execute(int *argc, char **argv[]);
    void LaunchProcess(const stringVector &launchArgs);
    void ConnectSimulation(const stringVector &launchArgs,
                           const std::string &simHost, int simPort,
                           const std::string &simSecurityKey);

protected:
    LauncherApplication();
    void ProcessArguments(int *argc, char **argv[]);
    void Connect(int *argc, char **argv[]);
    void MainLoop();
    bool ProcessInput();
    void TurnOnAlarm();
    void TurnOffAlarm();

    static void AlarmHandler(int);
private:
    static LauncherApplication *instance;

    ParentProcess               parent;
    Xfer                        xfer;

    QuitRPC                     quitRPC;
    KeepAliveRPC                keepAliveRPC;
    LauncherState               launcherstate;
    RPCExecutor<QuitRPC>       *quitExecutor;
    RPCExecutor<KeepAliveRPC>  *keepAliveExecutor;
    RPCExecutor<LaunchRPC>     *launchExecutor;
    RPCExecutor<ConnectSimRPC> *connectSimExecutor;

    bool                        useSSHTunneling;
    bool                        keepGoing;
    int                         timeout;
    std::vector<Connection*>    childOutput;
    LaunchService               launch;
};

#endif
