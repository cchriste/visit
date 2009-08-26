// ***************************************************************************
//
// Copyright (c) 2000 - 2009, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400124
// All rights reserved.
//
// This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
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
//    documentation and/or other materials provided with the distribution.
//  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
//    be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
// LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
// DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
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
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class HostProfile extends AttributeSubject
{
    private static int numAdditionalAttributes = 44;

    // Enum values
    public final static int CLIENTHOSTDETERMINATION_MACHINENAME = 0;
    public final static int CLIENTHOSTDETERMINATION_MANUALLYSPECIFIED = 1;
    public final static int CLIENTHOSTDETERMINATION_PARSEDFROMSSHCLIENT = 2;


    public HostProfile()
    {
        super(numAdditionalAttributes);

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
        sublaunchPreCmdSet = false;
        sublaunchPreCmd = new String("");
        sublaunchPostCmdSet = false;
        sublaunchPostCmd = new String("");
        hostAliases = new String("");
        hostNickname = new String("");
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
        tunnelSSH = false;
    }

    public HostProfile(int nMoreFields)
    {
        super(numAdditionalAttributes + nMoreFields);

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
        sublaunchPreCmdSet = false;
        sublaunchPreCmd = new String("");
        sublaunchPostCmdSet = false;
        sublaunchPostCmd = new String("");
        hostAliases = new String("");
        hostNickname = new String("");
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
        tunnelSSH = false;
    }

    public HostProfile(HostProfile obj)
    {
        super(numAdditionalAttributes);

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
        sublaunchPreCmdSet = obj.sublaunchPreCmdSet;
        sublaunchPreCmd = new String(obj.sublaunchPreCmd);
        sublaunchPostCmdSet = obj.sublaunchPostCmdSet;
        sublaunchPostCmd = new String(obj.sublaunchPostCmd);
        hostAliases = new String(obj.hostAliases);
        hostNickname = new String(obj.hostNickname);
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
        tunnelSSH = obj.tunnelSSH;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return numAdditionalAttributes;
    }

    public boolean equals(HostProfile obj)
    {
        int i;

        // Compare the elements in the arguments vector.
        boolean arguments_equal = (obj.arguments.size() == arguments.size());
        for(i = 0; (i < arguments.size()) && arguments_equal; ++i)
        {
            // Make references to String from Object.
            String arguments1 = (String)arguments.elementAt(i);
            String arguments2 = (String)obj.arguments.elementAt(i);
            arguments_equal = arguments1.equals(arguments2);
        }
        // Create the return value
        return ((profileName.equals(obj.profileName)) &&
                (host.equals(obj.host)) &&
                (userName.equals(obj.userName)) &&
                (timeout == obj.timeout) &&
                (numProcessors == obj.numProcessors) &&
                (numNodesSet == obj.numNodesSet) &&
                (numNodes == obj.numNodes) &&
                (partitionSet == obj.partitionSet) &&
                (partition.equals(obj.partition)) &&
                (bankSet == obj.bankSet) &&
                (bank.equals(obj.bank)) &&
                (timeLimitSet == obj.timeLimitSet) &&
                (timeLimit.equals(obj.timeLimit)) &&
                (launchMethodSet == obj.launchMethodSet) &&
                (launchMethod.equals(obj.launchMethod)) &&
                (forceStatic == obj.forceStatic) &&
                (forceDynamic == obj.forceDynamic) &&
                (active == obj.active) &&
                arguments_equal &&
                (parallel == obj.parallel) &&
                (launchArgsSet == obj.launchArgsSet) &&
                (launchArgs.equals(obj.launchArgs)) &&
                (sublaunchArgsSet == obj.sublaunchArgsSet) &&
                (sublaunchArgs.equals(obj.sublaunchArgs)) &&
                (sublaunchPreCmdSet == obj.sublaunchPreCmdSet) &&
                (sublaunchPreCmd.equals(obj.sublaunchPreCmd)) &&
                (sublaunchPostCmdSet == obj.sublaunchPostCmdSet) &&
                (sublaunchPostCmd.equals(obj.sublaunchPostCmd)) &&
                (hostAliases.equals(obj.hostAliases)) &&
                (hostNickname.equals(obj.hostNickname)) &&
                (shareOneBatchJob == obj.shareOneBatchJob) &&
                (sshPortSpecified == obj.sshPortSpecified) &&
                (sshPort == obj.sshPort) &&
                (clientHostDetermination == obj.clientHostDetermination) &&
                (manualClientHostName.equals(obj.manualClientHostName)) &&
                (machinefileSet == obj.machinefileSet) &&
                (machinefile.equals(obj.machinefile)) &&
                (visitSetsUpEnv == obj.visitSetsUpEnv) &&
                (canDoHWAccel == obj.canDoHWAccel) &&
                (havePreCommand == obj.havePreCommand) &&
                (hwAccelPreCommand.equals(obj.hwAccelPreCommand)) &&
                (havePostCommand == obj.havePostCommand) &&
                (hwAccelPostCommand.equals(obj.hwAccelPostCommand)) &&
                (tunnelSSH == obj.tunnelSSH));
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

    public void SetSublaunchPreCmdSet(boolean sublaunchPreCmdSet_)
    {
        sublaunchPreCmdSet = sublaunchPreCmdSet_;
        Select(24);
    }

    public void SetSublaunchPreCmd(String sublaunchPreCmd_)
    {
        sublaunchPreCmd = sublaunchPreCmd_;
        Select(25);
    }

    public void SetSublaunchPostCmdSet(boolean sublaunchPostCmdSet_)
    {
        sublaunchPostCmdSet = sublaunchPostCmdSet_;
        Select(26);
    }

    public void SetSublaunchPostCmd(String sublaunchPostCmd_)
    {
        sublaunchPostCmd = sublaunchPostCmd_;
        Select(27);
    }

    public void SetHostAliases(String hostAliases_)
    {
        hostAliases = hostAliases_;
        Select(28);
    }

    public void SetHostNickname(String hostNickname_)
    {
        hostNickname = hostNickname_;
        Select(29);
    }

    public void SetShareOneBatchJob(boolean shareOneBatchJob_)
    {
        shareOneBatchJob = shareOneBatchJob_;
        Select(30);
    }

    public void SetSshPortSpecified(boolean sshPortSpecified_)
    {
        sshPortSpecified = sshPortSpecified_;
        Select(31);
    }

    public void SetSshPort(int sshPort_)
    {
        sshPort = sshPort_;
        Select(32);
    }

    public void SetClientHostDetermination(int clientHostDetermination_)
    {
        clientHostDetermination = clientHostDetermination_;
        Select(33);
    }

    public void SetManualClientHostName(String manualClientHostName_)
    {
        manualClientHostName = manualClientHostName_;
        Select(34);
    }

    public void SetMachinefileSet(boolean machinefileSet_)
    {
        machinefileSet = machinefileSet_;
        Select(35);
    }

    public void SetMachinefile(String machinefile_)
    {
        machinefile = machinefile_;
        Select(36);
    }

    public void SetVisitSetsUpEnv(boolean visitSetsUpEnv_)
    {
        visitSetsUpEnv = visitSetsUpEnv_;
        Select(37);
    }

    public void SetCanDoHWAccel(boolean canDoHWAccel_)
    {
        canDoHWAccel = canDoHWAccel_;
        Select(38);
    }

    public void SetHavePreCommand(boolean havePreCommand_)
    {
        havePreCommand = havePreCommand_;
        Select(39);
    }

    public void SetHwAccelPreCommand(String hwAccelPreCommand_)
    {
        hwAccelPreCommand = hwAccelPreCommand_;
        Select(40);
    }

    public void SetHavePostCommand(boolean havePostCommand_)
    {
        havePostCommand = havePostCommand_;
        Select(41);
    }

    public void SetHwAccelPostCommand(String hwAccelPostCommand_)
    {
        hwAccelPostCommand = hwAccelPostCommand_;
        Select(42);
    }

    public void SetTunnelSSH(boolean tunnelSSH_)
    {
        tunnelSSH = tunnelSSH_;
        Select(43);
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
    public boolean GetSublaunchPreCmdSet() { return sublaunchPreCmdSet; }
    public String  GetSublaunchPreCmd() { return sublaunchPreCmd; }
    public boolean GetSublaunchPostCmdSet() { return sublaunchPostCmdSet; }
    public String  GetSublaunchPostCmd() { return sublaunchPostCmd; }
    public String  GetHostAliases() { return hostAliases; }
    public String  GetHostNickname() { return hostNickname; }
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
    public boolean GetTunnelSSH() { return tunnelSSH; }

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
            buf.WriteBool(sublaunchPreCmdSet);
        if(WriteSelect(25, buf))
            buf.WriteString(sublaunchPreCmd);
        if(WriteSelect(26, buf))
            buf.WriteBool(sublaunchPostCmdSet);
        if(WriteSelect(27, buf))
            buf.WriteString(sublaunchPostCmd);
        if(WriteSelect(28, buf))
            buf.WriteString(hostAliases);
        if(WriteSelect(29, buf))
            buf.WriteString(hostNickname);
        if(WriteSelect(30, buf))
            buf.WriteBool(shareOneBatchJob);
        if(WriteSelect(31, buf))
            buf.WriteBool(sshPortSpecified);
        if(WriteSelect(32, buf))
            buf.WriteInt(sshPort);
        if(WriteSelect(33, buf))
            buf.WriteInt(clientHostDetermination);
        if(WriteSelect(34, buf))
            buf.WriteString(manualClientHostName);
        if(WriteSelect(35, buf))
            buf.WriteBool(machinefileSet);
        if(WriteSelect(36, buf))
            buf.WriteString(machinefile);
        if(WriteSelect(37, buf))
            buf.WriteBool(visitSetsUpEnv);
        if(WriteSelect(38, buf))
            buf.WriteBool(canDoHWAccel);
        if(WriteSelect(39, buf))
            buf.WriteBool(havePreCommand);
        if(WriteSelect(40, buf))
            buf.WriteString(hwAccelPreCommand);
        if(WriteSelect(41, buf))
            buf.WriteBool(havePostCommand);
        if(WriteSelect(42, buf))
            buf.WriteString(hwAccelPostCommand);
        if(WriteSelect(43, buf))
            buf.WriteBool(tunnelSSH);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
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
            SetSublaunchPreCmdSet(buf.ReadBool());
            break;
        case 25:
            SetSublaunchPreCmd(buf.ReadString());
            break;
        case 26:
            SetSublaunchPostCmdSet(buf.ReadBool());
            break;
        case 27:
            SetSublaunchPostCmd(buf.ReadString());
            break;
        case 28:
            SetHostAliases(buf.ReadString());
            break;
        case 29:
            SetHostNickname(buf.ReadString());
            break;
        case 30:
            SetShareOneBatchJob(buf.ReadBool());
            break;
        case 31:
            SetSshPortSpecified(buf.ReadBool());
            break;
        case 32:
            SetSshPort(buf.ReadInt());
            break;
        case 33:
            SetClientHostDetermination(buf.ReadInt());
            break;
        case 34:
            SetManualClientHostName(buf.ReadString());
            break;
        case 35:
            SetMachinefileSet(buf.ReadBool());
            break;
        case 36:
            SetMachinefile(buf.ReadString());
            break;
        case 37:
            SetVisitSetsUpEnv(buf.ReadBool());
            break;
        case 38:
            SetCanDoHWAccel(buf.ReadBool());
            break;
        case 39:
            SetHavePreCommand(buf.ReadBool());
            break;
        case 40:
            SetHwAccelPreCommand(buf.ReadString());
            break;
        case 41:
            SetHavePostCommand(buf.ReadBool());
            break;
        case 42:
            SetHwAccelPostCommand(buf.ReadString());
            break;
        case 43:
            SetTunnelSSH(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("profileName", profileName, indent) + "\n";
        str = str + stringToString("host", host, indent) + "\n";
        str = str + stringToString("userName", userName, indent) + "\n";
        str = str + intToString("timeout", timeout, indent) + "\n";
        str = str + intToString("numProcessors", numProcessors, indent) + "\n";
        str = str + boolToString("numNodesSet", numNodesSet, indent) + "\n";
        str = str + intToString("numNodes", numNodes, indent) + "\n";
        str = str + boolToString("partitionSet", partitionSet, indent) + "\n";
        str = str + stringToString("partition", partition, indent) + "\n";
        str = str + boolToString("bankSet", bankSet, indent) + "\n";
        str = str + stringToString("bank", bank, indent) + "\n";
        str = str + boolToString("timeLimitSet", timeLimitSet, indent) + "\n";
        str = str + stringToString("timeLimit", timeLimit, indent) + "\n";
        str = str + boolToString("launchMethodSet", launchMethodSet, indent) + "\n";
        str = str + stringToString("launchMethod", launchMethod, indent) + "\n";
        str = str + boolToString("forceStatic", forceStatic, indent) + "\n";
        str = str + boolToString("forceDynamic", forceDynamic, indent) + "\n";
        str = str + boolToString("active", active, indent) + "\n";
        str = str + stringVectorToString("arguments", arguments, indent) + "\n";
        str = str + boolToString("parallel", parallel, indent) + "\n";
        str = str + boolToString("launchArgsSet", launchArgsSet, indent) + "\n";
        str = str + stringToString("launchArgs", launchArgs, indent) + "\n";
        str = str + boolToString("sublaunchArgsSet", sublaunchArgsSet, indent) + "\n";
        str = str + stringToString("sublaunchArgs", sublaunchArgs, indent) + "\n";
        str = str + boolToString("sublaunchPreCmdSet", sublaunchPreCmdSet, indent) + "\n";
        str = str + stringToString("sublaunchPreCmd", sublaunchPreCmd, indent) + "\n";
        str = str + boolToString("sublaunchPostCmdSet", sublaunchPostCmdSet, indent) + "\n";
        str = str + stringToString("sublaunchPostCmd", sublaunchPostCmd, indent) + "\n";
        str = str + stringToString("hostAliases", hostAliases, indent) + "\n";
        str = str + stringToString("hostNickname", hostNickname, indent) + "\n";
        str = str + boolToString("shareOneBatchJob", shareOneBatchJob, indent) + "\n";
        str = str + boolToString("sshPortSpecified", sshPortSpecified, indent) + "\n";
        str = str + intToString("sshPort", sshPort, indent) + "\n";
        str = str + indent + "clientHostDetermination = ";
        if(clientHostDetermination == CLIENTHOSTDETERMINATION_MACHINENAME)
            str = str + "CLIENTHOSTDETERMINATION_MACHINENAME";
        if(clientHostDetermination == CLIENTHOSTDETERMINATION_MANUALLYSPECIFIED)
            str = str + "CLIENTHOSTDETERMINATION_MANUALLYSPECIFIED";
        if(clientHostDetermination == CLIENTHOSTDETERMINATION_PARSEDFROMSSHCLIENT)
            str = str + "CLIENTHOSTDETERMINATION_PARSEDFROMSSHCLIENT";
        str = str + "\n";
        str = str + stringToString("manualClientHostName", manualClientHostName, indent) + "\n";
        str = str + boolToString("machinefileSet", machinefileSet, indent) + "\n";
        str = str + stringToString("machinefile", machinefile, indent) + "\n";
        str = str + boolToString("visitSetsUpEnv", visitSetsUpEnv, indent) + "\n";
        str = str + boolToString("canDoHWAccel", canDoHWAccel, indent) + "\n";
        str = str + boolToString("havePreCommand", havePreCommand, indent) + "\n";
        str = str + stringToString("hwAccelPreCommand", hwAccelPreCommand, indent) + "\n";
        str = str + boolToString("havePostCommand", havePostCommand, indent) + "\n";
        str = str + stringToString("hwAccelPostCommand", hwAccelPostCommand, indent) + "\n";
        str = str + boolToString("tunnelSSH", tunnelSSH, indent) + "\n";
        return str;
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
    private boolean sublaunchPreCmdSet;
    private String  sublaunchPreCmd;
    private boolean sublaunchPostCmdSet;
    private String  sublaunchPostCmd;
    private String  hostAliases;
    private String  hostNickname;
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
    private boolean tunnelSSH;
}

