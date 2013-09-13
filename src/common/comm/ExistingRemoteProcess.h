/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
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

#ifndef EXISTING_REMOTE_PROCESS_H
#define EXISTING_REMOTE_PROCESS_H
#include <comm_exports.h>
#include <RemoteProcess.h>
#include <ConnectCallback.h>

// ****************************************************************************
// Class: ExistingRemoteProcess
//
// Purpose:
//   Connects to a remote process that is already running.
//
// Notes:      The user-supplied ConnectCallback function is used to
//             initiate the connection with the remote process.
//
// Programmer: Brad Whitlock
// Creation:   Mon Nov 20 13:09:38 PST 2000
//
// Modifications:
//   Brad Whitlock, Mon May 5 17:24:01 PST 2003
//   I added the createAsThoughLocal argument to the Open method to force
//   the argument list used to create the process to be the local form
//   even though the process may be launched on a remote machine.
//
//   Jeremy Meredith, Thu Oct  9 14:04:08 PDT 2003
//   Added ability to manually specify a client host name or to have it
//   parsed from the SSH_CLIENT (or related) environment variables.  Added
//   ability to specify an SSH port.
//
//   Jeremy Meredith, Thu May 24 11:10:15 EDT 2007
//   Added SSH tunneling argument to Open.
//
//   Jeremy Meredith, Thu Feb 18 15:25:27 EST 2010
//   Split HostProfile int MachineProfile and LaunchProfile.
//
//   Eric Brugger, Mon May  2 16:44:20 PDT 2011
//   I added the ability to use a gateway machine when connecting to a
//   remote host.
//
//   Brad Whitlock, Tue Jun  5 17:24:03 PDT 2012
//   Pass in MachineProfile to Open.
//
// ****************************************************************************

class COMM_API ExistingRemoteProcess : public RemoteProcess
{
public:
    ExistingRemoteProcess(const std::string &rProgram);
    virtual ~ExistingRemoteProcess();
    virtual bool Open(const MachineProfile &profile,
                      int numRead, int numWrite,
                      bool createAsThoughLocal = false);
    void SetConnectCallback(ConnectCallback *cb);
    void SetConnectCallbackData(void *data);
private:
    ConnectCallback *connectCallback;
    void            *connectCallbackData;
};

#endif
