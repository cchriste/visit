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

#ifndef REMOTE_PROCESS_H
#define REMOTE_PROCESS_H
#include <visit-config.h>
#include <comm_exports.h>
#if defined(_WIN32)
#include <winsock2.h>
#else
#include <netdb.h>
#include <netinet/in.h>
#endif
#include <vectortypes.h>
#include <MachineProfile.h>
#include <map>

class Connection;

// ****************************************************************************
// Class: RemoteProcess
//
// Purpose:
//   Starts a process on a remote machine using ssh and connects some
//   sockets to it.
//
// Notes:      
//
// Programmer: Brad Whitlock
// Creation:   Thu Jul 13 15:47:04 PST 2000
//
// Modifications:
//    Jeremy Meredith, Tue Aug  8 13:49:42 PDT 2000
//    Changed it to allow more than one read/write socket.   
//
//    Brad Whitlock, Fri Aug 25 10:49:10 PDT 2000
//    Changed it so it can launch local processes instead of always
//    having to go through ssh.
//
//    Brad Whitlock, Mon Nov 20 16:12:38 PST 2000
//    Changed private permissions to protected.
//
//    Jeremy Meredith, Fri Apr 27 11:41:05 PDT 2001
//    Added authentication callback function.
//
//    Brad Whitlock, Mon Sep 24 11:26:13 PDT 2001
//    Added GetLocalHostName() method.
//
//    Brad Whitlock, Thu Feb 21 10:02:10 PDT 2002
//    Added GetLocalUserName() method.
//
//    Brad Whitlock, Wed Apr 3 12:01:38 PDT 2002
//    Added WaitForTermination() method.
//
//    Brad Whitlock, Mon Aug 26 12:47:38 PDT 2002
//    I added methods to return the SSH program being used as well as any
//    arguments that it might use.
//
//    Brad Whitlock, Thu Aug 29 17:46:48 PST 2002
//    I added another argument to the authentication callback so we can
//    pass the remote login name.
//
//    Brad Whitlock, Thu Sep 26 16:48:55 PST 2002
//    I added a progress callback and I split AcceptSocket into two methods.
//
//    Brad Whitlock, Mon Dec 16 14:28:02 PST 2002
//    I added securityKey and a method to close the listen socket.
//
//    Brad Whitlock, Mon May 5 13:24:47 PST 2003
//    I refactored the code a little to simplify handling of command lines.
//
//    Jeremy Meredith, Thu Jul  3 14:49:23 PDT 2003
//    Added ability to disable ptys.
//
//    Jeremy Meredith, Thu Oct  9 13:56:56 PDT 2003
//    Added ability to manually specify a client host name or to have it
//    parsed from the SSH_CLIENT (or related) environment variables.  Added
//    ability to specify an SSH port.
//
//    Jeremy Meredith, Thu May 24 11:08:03 EDT 2007
//    Added support for SSH port tunneling.  Also made CreateCommandLine
//    non-const; we need to store the port forwarding map at that time.
//
//    Thomas R. Treadway, Mon Oct  8 13:27:42 PDT 2007
//    Backing out SSH tunneling on Panther (MacOS X 10.3)
//
//    Jeremy Meredith, Thu Feb 18 15:25:27 EST 2010
//    Split HostProfile int MachineProfile and LaunchProfile.
//
//    Eric Brugger, Mon May  2 16:35:56 PDT 2011
//    I added the ability to use a gateway machine when connecting to a
//    remote host.
//
//    Eric Brugger, Mon Sep 26 16:40:20 PDT 2011
//    I modified the remote launching to pass the remote user name to the
//    ssh command to the gateway machine instead of to the ssh command to
//    the remote machine.
//
//    Brad Whitlock, Fri Jan 13 14:55:55 PST 2012
//    I modifed the launching routines to handle SSH separately from the 
//    regular command line so we can better support gateway machines with
//    SSH tunneling.
//
//    Brad Whitlock, Tue Jun  5 17:23:25 PDT 2012
//    Change Open so it takes a MachineProfile instead of a bunch of separate
//    options.
//
//    Brad Whitlock, Wed Jun 13 11:12:25 PDT 2012
//    I added KillProcess.
//
// ****************************************************************************

class COMM_API RemoteProcess
{
public:
    RemoteProcess(const std::string &rProgram);
    virtual ~RemoteProcess();
    virtual bool Open(const MachineProfile &profile,
                      int numRead, int numWrite,
                      bool createAsThoughLocal = false);
    void WaitForTermination();
    void AddArgument(const std::string &arg);

    const std::string &GetLocalHostName() const;
    const std::string &GetLocalUserName() const;
    Connection *GetReadConnection(int i=0) const;
    Connection *GetWriteConnection(int i=0) const;
    int  GetProcessId() const;
    void SetProgressCallback(bool (*)(void *, int), void *);
    std::map<int,int> GetPortTunnelMap() { return portTunnelMap; }

    static void SetAuthenticationCallback(void (*)(const char *, const char *, int));
    static void SetChangeUserNameCallback(bool (*)(const std::string &,std::string&));
    static void DisablePTY();
    static void SetCustomConnectionCallback(Connection* (*)(int, void *), void* cbData);
protected:
    bool StartMakingConnection(const std::string &remoteHost, int numRead, int numWrite);
    void FinishMakingConnection(int numRead, int numWrite);
    const char *SecureShell() const;
    const char *SecureShellArgs() const;
    bool CallProgressCallback(int stage);
    bool HostIsLocal(const std::string &rHost) const;
    void CreatePortNumbers(int *local, int *remote, int *gateway, int nPorts) const;
    void CreateSSHCommandLine(stringVector &args, const MachineProfile &profile);
    void CreateCommandLine(stringVector &args, const MachineProfile &profile, 
                           int numRead, int numWrite);
    virtual void Launch(const stringVector &);
    void LaunchLocal(const stringVector &);
    void LaunchRemote(const std::string &host, const std::string &remoteUserName, 
                      const stringVector &);
    void KillProcess();
protected:
    int                      listenPortNum;
    std::string              localHost, localUserName;
    std::string              securityKey;
private:
    bool GetSocketAndPort();
    int  AcceptSocket();
    int  SingleThreadedAcceptSocket();
    int  MultiThreadedAcceptSocket();
    void CloseListenSocket();
    void ExchangeTypeRepresentations();
    char **CreateSplitCommandLine(const stringVector &args, int &argc) const;
    void DestroySplitCommandLine(char **args, int argc) const;
    char *StrDup(const std::string &) const;
private:
    DESCRIPTOR               listenSocketNum;
    struct sockaddr_in       sin;
    std::string              remoteProgram;
    std::vector<std::string> argList;
    int                      remoteProgramPid;
    Connection             **readConnections, **writeConnections;
    int                      nReadConnections, nWriteConnections;
    bool                   (*progressCallback)(void *, int);
    void                    *progressCallbackData;
    std::map<int,int>        portTunnelMap;

    static void            (*getAuthentication)(const char *, const char *, int);
    static bool            (*changeUsername)(const std::string &, std::string&);
    static bool              disablePTY;
    static Connection*     (*customConnectionCallback)(int,void*);
    static void             *customConnectionCallbackData;
};

#endif
