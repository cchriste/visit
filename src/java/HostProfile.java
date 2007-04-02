// ***************************************************************************
//
// Copyright (c) 2000 - 2007, The Regents of the University of California
// Produced at the Lawrence Livermore National Laboratory
// All rights reserved.
//
// This file is part of VisIt. For details, see http://www.llnl.gov/visit/. The
// full copyright notice is contained in the file COPYRIGHT located at the root
// of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
//
// Redistribution  and  use  in  source  and  binary  forms,  with  or  without
// modification, are permitted provided that the following conditions are met:
//
//  - Redistributions of  source code must  retain the above  copyright notice,
//    this list of conditions and the disclaimer below.
//  - Redistributions in binary form must reproduce the above copyright notice,
//    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
//    documentation and/or materials provided with the distribution.
//  - Neither the name of the UC/LLNL nor  the names of its contributors may be
//    used to  endorse or  promote products derived from  this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED.  IN  NO  EVENT  SHALL  THE  REGENTS  OF  THE  UNIVERSITY OF
// CALIFORNIA, THE U.S.  DEPARTMENT  OF  ENERGY OR CONTRIBUTORS BE  LIABLE  FOR
// ANY  DIRECT,  INDIRECT,  INCIDENTAL,  SPECIAL,  EXEMPLARY,  OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
// SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
// CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
// LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
// OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ***************************************************************************

package llnl.visit;

import java.util.Vector;

// ****************************************************************************
// Class: HostProfile
//
// Purpose:
//    This class contains information needed to launch a remote,VisIt engine.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Mar 14 17:54:55 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class HostProfile extends AttributeSubject
{
    // Enum values
    public final static int CLIENTHOSTDETERMINATION_MACHINENAME = 0;
    public final static int CLIENTHOSTDETERMINATION_MANUALLYSPECIFIED = 1;
    public final static int CLIENTHOSTDETERMINATION_PARSEDFROMSSHCLIENT = 2;


    public HostProfile()
    {
        super(38);

        profileName = new String("notset");
        host = new String("localhost");
        userName = new String("notset");
        timeout = 480;
        numProcessors = 1;
        numNodesSet = false;
        numNodes = 1;
        partitionSet = false;
        partition = new String("");
        bankSet = false;
        bank = new String("");
        timeLimitSet = false;
        timeLimit = new String("");
        launchMethodSet = false;
        launchMethod = new String("");
        forceStatic = true;
        forceDynamic = false;
        active = false;
        arguments = new Vector();
        parallel = false;
        launchArgsSet = false;
        launchArgs = new String("");
        sublaunchArgsSet = false;
        sublaunchArgs = new String("");
        hostAliases = new String("");
        shareOneBatchJob = false;
        sshPortSpecified = false;
        sshPort = 22;
        clientHostDetermination = CLIENTHOSTDETERMINATION_MACHINENAME;
        manualClientHostName = new String("");
        machinefileSet = false;
        machinefile = new String("");
        visitSetsUpEnv = false;
        canDoHWAccel = false;
        havePreCommand = false;
        hwAccelPreCommand = new String("");
        havePostCommand = false;
        hwAccelPostCommand = new String("");
    }

    public HostProfile(HostProfile obj)
    {
        super(38);

        int i;

        profileName = new String(obj.profileName);
        host = new String(obj.host);
        userName = new String(obj.userName);
        timeout = obj.timeout;
        numProcessors = obj.numProcessors;
        numNodesSet = obj.numNodesSet;
        numNodes = obj.numNodes;
        partitionSet = obj.partitionSet;
        partition = new String(obj.partition);
        bankSet = obj.bankSet;
        bank = new String(obj.bank);
        timeLimitSet = obj.timeLimitSet;
        timeLimit = new String(obj.timeLimit);
        launchMethodSet = obj.launchMethodSet;
        launchMethod = new String(obj.launchMethod);
        forceStatic = obj.forceStatic;
        forceDynamic = obj.forceDynamic;
        active = obj.active;
        arguments = new Vector(obj.arguments.size());
        for(i = 0; i < obj.arguments.size(); ++i)
            arguments.addElement(new String((String)obj.arguments.elementAt(i)));

        parallel = obj.parallel;
        launchArgsSet = obj.launchArgsSet;
        launchArgs = new String(obj.launchArgs);
        sublaunchArgsSet = obj.sublaunchArgsSet;
        sublaunchArgs = new String(obj.sublaunchArgs);
        hostAliases = new String(obj.hostAliases);
        shareOneBatchJob = obj.shareOneBatchJob;
        sshPortSpecified = obj.sshPortSpecified;
        sshPort = obj.sshPort;
        clientHostDetermination = obj.clientHostDetermination;
        manualClientHostName = new String(obj.manualClientHostName);
        machinefileSet = obj.machinefileSet;
        machinefile = new String(obj.machinefile);
        visitSetsUpEnv = obj.visitSetsUpEnv;
        canDoHWAccel = obj.canDoHWAccel;
        havePreCommand = obj.havePreCommand;
        hwAccelPreCommand = new String(obj.hwAccelPreCommand);
        havePostCommand = obj.havePostCommand;
        hwAccelPostCommand = new String(obj.hwAccelPostCommand);

        SelectAll();
    }

    public boolean equals(HostProfile obj)
    {
        int i;

        // Create the return value
        return ((profileName == obj.profileName) &&
                (host == obj.host) &&
                (userName == obj.userName) &&
                (timeout == obj.timeout) &&
                (numProcessors == obj.numProcessors) &&
                (numNodesSet == obj.numNodesSet) &&
                (numNodes == obj.numNodes) &&
                (partitionSet == obj.partitionSet) &&
                (partition == obj.partition) &&
                (bankSet == obj.bankSet) &&
                (bank == obj.bank) &&
                (timeLimitSet == obj.timeLimitSet) &&
                (timeLimit == obj.timeLimit) &&
                (launchMethodSet == obj.launchMethodSet) &&
                (launchMethod == obj.launchMethod) &&
                (forceStatic == obj.forceStatic) &&
                (forceDynamic == obj.forceDynamic) &&
                (active == obj.active) &&
                (arguments == obj.arguments) &&
                (parallel == obj.parallel) &&
                (launchArgsSet == obj.launchArgsSet) &&
                (launchArgs == obj.launchArgs) &&
                (sublaunchArgsSet == obj.sublaunchArgsSet) &&
                (sublaunchArgs == obj.sublaunchArgs) &&
                (hostAliases == obj.hostAliases) &&
                (shareOneBatchJob == obj.shareOneBatchJob) &&
                (sshPortSpecified == obj.sshPortSpecified) &&
                (sshPort == obj.sshPort) &&
                (clientHostDetermination == obj.clientHostDetermination) &&
                (manualClientHostName == obj.manualClientHostName) &&
                (machinefileSet == obj.machinefileSet) &&
                (machinefile == obj.machinefile) &&
                (visitSetsUpEnv == obj.visitSetsUpEnv) &&
                (canDoHWAccel == obj.canDoHWAccel) &&
                (havePreCommand == obj.havePreCommand) &&
                (hwAccelPreCommand == obj.hwAccelPreCommand) &&
                (havePostCommand == obj.havePostCommand) &&
                (hwAccelPostCommand == obj.hwAccelPostCommand));
    }

    // Property setting methods
    public void SetProfileName(String profileName_)
    {
        profileName = profileName_;
        Select(0);
    }

    public void SetHost(String host_)
    {
        host = host_;
        Select(1);
    }

    public void SetUserName(String userName_)
    {
        userName = userName_;
        Select(2);
    }

    public void SetTimeout(int timeout_)
    {
        timeout = timeout_;
        Select(3);
    }

    public void SetNumProcessors(int numProcessors_)
    {
        numProcessors = numProcessors_;
        Select(4);
    }

    public void SetNumNodesSet(boolean numNodesSet_)
    {
        numNodesSet = numNodesSet_;
        Select(5);
    }

    public void SetNumNodes(int numNodes_)
    {
        numNodes = numNodes_;
        Select(6);
    }

    public void SetPartitionSet(boolean partitionSet_)
    {
        partitionSet = partitionSet_;
        Select(7);
    }

    public void SetPartition(String partition_)
    {
        partition = partition_;
        Select(8);
    }

    public void SetBankSet(boolean bankSet_)
    {
        bankSet = bankSet_;
        Select(9);
    }

    public void SetBank(String bank_)
    {
        bank = bank_;
        Select(10);
    }

    public void SetTimeLimitSet(boolean timeLimitSet_)
    {
        timeLimitSet = timeLimitSet_;
        Select(11);
    }

    public void SetTimeLimit(String timeLimit_)
    {
        timeLimit = timeLimit_;
        Select(12);
    }

    public void SetLaunchMethodSet(boolean launchMethodSet_)
    {
        launchMethodSet = launchMethodSet_;
        Select(13);
    }

    public void SetLaunchMethod(String launchMethod_)
    {
        launchMethod = launchMethod_;
        Select(14);
    }

    public void SetForceStatic(boolean forceStatic_)
    {
        forceStatic = forceStatic_;
        Select(15);
    }

    public void SetForceDynamic(boolean forceDynamic_)
    {
        forceDynamic = forceDynamic_;
        Select(16);
    }

    public void SetActive(boolean active_)
    {
        active = active_;
        Select(17);
    }

    public void SetArguments(Vector arguments_)
    {
        arguments = arguments_;
        Select(18);
    }

    public void SetParallel(boolean parallel_)
    {
        parallel = parallel_;
        Select(19);
    }

    public void SetLaunchArgsSet(boolean launchArgsSet_)
    {
        launchArgsSet = launchArgsSet_;
        Select(20);
    }

    public void SetLaunchArgs(String launchArgs_)
    {
        launchArgs = launchArgs_;
        Select(21);
    }

    public void SetSublaunchArgsSet(boolean sublaunchArgsSet_)
    {
        sublaunchArgsSet = sublaunchArgsSet_;
        Select(22);
    }

    public void SetSublaunchArgs(String sublaunchArgs_)
    {
        sublaunchArgs = sublaunchArgs_;
        Select(23);
    }

    public void SetHostAliases(String hostAliases_)
    {
        hostAliases = hostAliases_;
        Select(24);
    }

    public void SetShareOneBatchJob(boolean shareOneBatchJob_)
    {
        shareOneBatchJob = shareOneBatchJob_;
        Select(25);
    }

    public void SetSshPortSpecified(boolean sshPortSpecified_)
    {
        sshPortSpecified = sshPortSpecified_;
        Select(26);
    }

    public void SetSshPort(int sshPort_)
    {
        sshPort = sshPort_;
        Select(27);
    }

    public void SetClientHostDetermination(int clientHostDetermination_)
    {
        clientHostDetermination = clientHostDetermination_;
        Select(28);
    }

    public void SetManualClientHostName(String manualClientHostName_)
    {
        manualClientHostName = manualClientHostName_;
        Select(29);
    }

    public void SetMachinefileSet(boolean machinefileSet_)
    {
        machinefileSet = machinefileSet_;
        Select(30);
    }

    public void SetMachinefile(String machinefile_)
    {
        machinefile = machinefile_;
        Select(31);
    }

    public void SetVisitSetsUpEnv(boolean visitSetsUpEnv_)
    {
        visitSetsUpEnv = visitSetsUpEnv_;
        Select(32);
    }

    public void SetCanDoHWAccel(boolean canDoHWAccel_)
    {
        canDoHWAccel = canDoHWAccel_;
        Select(33);
    }

    public void SetHavePreCommand(boolean havePreCommand_)
    {
        havePreCommand = havePreCommand_;
        Select(34);
    }

    public void SetHwAccelPreCommand(String hwAccelPreCommand_)
    {
        hwAccelPreCommand = hwAccelPreCommand_;
        Select(35);
    }

    public void SetHavePostCommand(boolean havePostCommand_)
    {
        havePostCommand = havePostCommand_;
        Select(36);
    }

    public void SetHwAccelPostCommand(String hwAccelPostCommand_)
    {
        hwAccelPostCommand = hwAccelPostCommand_;
        Select(37);
    }

    // Property getting methods
    public String  GetProfileName() { return profileName; }
    public String  GetHost() { return host; }
    public String  GetUserName() { return userName; }
    public int     GetTimeout() { return timeout; }
    public int     GetNumProcessors() { return numProcessors; }
    public boolean GetNumNodesSet() { return numNodesSet; }
    public int     GetNumNodes() { return numNodes; }
    public boolean GetPartitionSet() { return partitionSet; }
    public String  GetPartition() { return partition; }
    public boolean GetBankSet() { return bankSet; }
    public String  GetBank() { return bank; }
    public boolean GetTimeLimitSet() { return timeLimitSet; }
    public String  GetTimeLimit() { return timeLimit; }
    public boolean GetLaunchMethodSet() { return launchMethodSet; }
    public String  GetLaunchMethod() { return launchMethod; }
    public boolean GetForceStatic() { return forceStatic; }
    public boolean GetForceDynamic() { return forceDynamic; }
    public boolean GetActive() { return active; }
    public Vector  GetArguments() { return arguments; }
    public boolean GetParallel() { return parallel; }
    public boolean GetLaunchArgsSet() { return launchArgsSet; }
    public String  GetLaunchArgs() { return launchArgs; }
    public boolean GetSublaunchArgsSet() { return sublaunchArgsSet; }
    public String  GetSublaunchArgs() { return sublaunchArgs; }
    public String  GetHostAliases() { return hostAliases; }
    public boolean GetShareOneBatchJob() { return shareOneBatchJob; }
    public boolean GetSshPortSpecified() { return sshPortSpecified; }
    public int     GetSshPort() { return sshPort; }
    public int     GetClientHostDetermination() { return clientHostDetermination; }
    public String  GetManualClientHostName() { return manualClientHostName; }
    public boolean GetMachinefileSet() { return machinefileSet; }
    public String  GetMachinefile() { return machinefile; }
    public boolean GetVisitSetsUpEnv() { return visitSetsUpEnv; }
    public boolean GetCanDoHWAccel() { return canDoHWAccel; }
    public boolean GetHavePreCommand() { return havePreCommand; }
    public String  GetHwAccelPreCommand() { return hwAccelPreCommand; }
    public boolean GetHavePostCommand() { return havePostCommand; }
    public String  GetHwAccelPostCommand() { return hwAccelPostCommand; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(profileName);
        if(WriteSelect(1, buf))
            buf.WriteString(host);
        if(WriteSelect(2, buf))
            buf.WriteString(userName);
        if(WriteSelect(3, buf))
            buf.WriteInt(timeout);
        if(WriteSelect(4, buf))
            buf.WriteInt(numProcessors);
        if(WriteSelect(5, buf))
            buf.WriteBool(numNodesSet);
        if(WriteSelect(6, buf))
            buf.WriteInt(numNodes);
        if(WriteSelect(7, buf))
            buf.WriteBool(partitionSet);
        if(WriteSelect(8, buf))
            buf.WriteString(partition);
        if(WriteSelect(9, buf))
            buf.WriteBool(bankSet);
        if(WriteSelect(10, buf))
            buf.WriteString(bank);
        if(WriteSelect(11, buf))
            buf.WriteBool(timeLimitSet);
        if(WriteSelect(12, buf))
            buf.WriteString(timeLimit);
        if(WriteSelect(13, buf))
            buf.WriteBool(launchMethodSet);
        if(WriteSelect(14, buf))
            buf.WriteString(launchMethod);
        if(WriteSelect(15, buf))
            buf.WriteBool(forceStatic);
        if(WriteSelect(16, buf))
            buf.WriteBool(forceDynamic);
        if(WriteSelect(17, buf))
            buf.WriteBool(active);
        if(WriteSelect(18, buf))
            buf.WriteStringVector(arguments);
        if(WriteSelect(19, buf))
            buf.WriteBool(parallel);
        if(WriteSelect(20, buf))
            buf.WriteBool(launchArgsSet);
        if(WriteSelect(21, buf))
            buf.WriteString(launchArgs);
        if(WriteSelect(22, buf))
            buf.WriteBool(sublaunchArgsSet);
        if(WriteSelect(23, buf))
            buf.WriteString(sublaunchArgs);
        if(WriteSelect(24, buf))
            buf.WriteString(hostAliases);
        if(WriteSelect(25, buf))
            buf.WriteBool(shareOneBatchJob);
        if(WriteSelect(26, buf))
            buf.WriteBool(sshPortSpecified);
        if(WriteSelect(27, buf))
            buf.WriteInt(sshPort);
        if(WriteSelect(28, buf))
            buf.WriteInt(clientHostDetermination);
        if(WriteSelect(29, buf))
            buf.WriteString(manualClientHostName);
        if(WriteSelect(30, buf))
            buf.WriteBool(machinefileSet);
        if(WriteSelect(31, buf))
            buf.WriteString(machinefile);
        if(WriteSelect(32, buf))
            buf.WriteBool(visitSetsUpEnv);
        if(WriteSelect(33, buf))
            buf.WriteBool(canDoHWAccel);
        if(WriteSelect(34, buf))
            buf.WriteBool(havePreCommand);
        if(WriteSelect(35, buf))
            buf.WriteString(hwAccelPreCommand);
        if(WriteSelect(36, buf))
            buf.WriteBool(havePostCommand);
        if(WriteSelect(37, buf))
            buf.WriteString(hwAccelPostCommand);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetProfileName(buf.ReadString());
                break;
            case 1:
                SetHost(buf.ReadString());
                break;
            case 2:
                SetUserName(buf.ReadString());
                break;
            case 3:
                SetTimeout(buf.ReadInt());
                break;
            case 4:
                SetNumProcessors(buf.ReadInt());
                break;
            case 5:
                SetNumNodesSet(buf.ReadBool());
                break;
            case 6:
                SetNumNodes(buf.ReadInt());
                break;
            case 7:
                SetPartitionSet(buf.ReadBool());
                break;
            case 8:
                SetPartition(buf.ReadString());
                break;
            case 9:
                SetBankSet(buf.ReadBool());
                break;
            case 10:
                SetBank(buf.ReadString());
                break;
            case 11:
                SetTimeLimitSet(buf.ReadBool());
                break;
            case 12:
                SetTimeLimit(buf.ReadString());
                break;
            case 13:
                SetLaunchMethodSet(buf.ReadBool());
                break;
            case 14:
                SetLaunchMethod(buf.ReadString());
                break;
            case 15:
                SetForceStatic(buf.ReadBool());
                break;
            case 16:
                SetForceDynamic(buf.ReadBool());
                break;
            case 17:
                SetActive(buf.ReadBool());
                break;
            case 18:
                SetArguments(buf.ReadStringVector());
                break;
            case 19:
                SetParallel(buf.ReadBool());
                break;
            case 20:
                SetLaunchArgsSet(buf.ReadBool());
                break;
            case 21:
                SetLaunchArgs(buf.ReadString());
                break;
            case 22:
                SetSublaunchArgsSet(buf.ReadBool());
                break;
            case 23:
                SetSublaunchArgs(buf.ReadString());
                break;
            case 24:
                SetHostAliases(buf.ReadString());
                break;
            case 25:
                SetShareOneBatchJob(buf.ReadBool());
                break;
            case 26:
                SetSshPortSpecified(buf.ReadBool());
                break;
            case 27:
                SetSshPort(buf.ReadInt());
                break;
            case 28:
                SetClientHostDetermination(buf.ReadInt());
                break;
            case 29:
                SetManualClientHostName(buf.ReadString());
                break;
            case 30:
                SetMachinefileSet(buf.ReadBool());
                break;
            case 31:
                SetMachinefile(buf.ReadString());
                break;
            case 32:
                SetVisitSetsUpEnv(buf.ReadBool());
                break;
            case 33:
                SetCanDoHWAccel(buf.ReadBool());
                break;
            case 34:
                SetHavePreCommand(buf.ReadBool());
                break;
            case 35:
                SetHwAccelPreCommand(buf.ReadString());
                break;
            case 36:
                SetHavePostCommand(buf.ReadBool());
                break;
            case 37:
                SetHwAccelPostCommand(buf.ReadString());
                break;
            }
        }
    }


    // Attributes
    private String  profileName;
    private String  host;
    private String  userName;
    private int     timeout;
    private int     numProcessors;
    private boolean numNodesSet;
    private int     numNodes;
    private boolean partitionSet;
    private String  partition;
    private boolean bankSet;
    private String  bank;
    private boolean timeLimitSet;
    private String  timeLimit;
    private boolean launchMethodSet;
    private String  launchMethod;
    private boolean forceStatic;
    private boolean forceDynamic;
    private boolean active;
    private Vector  arguments; // vector of String objects
    private boolean parallel;
    private boolean launchArgsSet;
    private String  launchArgs;
    private boolean sublaunchArgsSet;
    private String  sublaunchArgs;
    private String  hostAliases;
    private boolean shareOneBatchJob;
    private boolean sshPortSpecified;
    private int     sshPort;
    private int     clientHostDetermination;
    private String  manualClientHostName;
    private boolean machinefileSet;
    private String  machinefile;
    private boolean visitSetsUpEnv;
    private boolean canDoHWAccel;
    private boolean havePreCommand;
    private String  hwAccelPreCommand;
    private boolean havePostCommand;
    private String  hwAccelPostCommand;
}

