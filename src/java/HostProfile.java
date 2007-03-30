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
// Creation:   Fri May 16 10:54:57 PDT 2003
//
// Modifications:
//   
// ****************************************************************************

public class HostProfile extends AttributeSubject
{
    public HostProfile()
    {
        super(24);

        profileName = new String("notset");
        host = new String("localhost");
        userName = new String("notset");
        timeout = 240;
        numProcessors = 1;
        numNodesSet = false;
        numNodes = 1;
        partitionSet = false;
        partition = new String("(null)");
        bankSet = false;
        bank = new String("(null)");
        timeLimitSet = false;
        timeLimit = new String("(null)");
        launchMethodSet = false;
        launchMethod = new String("(null)");
        forceStatic = true;
        forceDynamic = false;
        active = false;
        arguments = new Vector();
        parallel = false;
        launchArgsSet = false;
        launchArgs = new String("(null)");
        hostAliases = new String("(null)");
        shareOneBatchJob = false;
    }

    public HostProfile(HostProfile obj)
    {
        super(24);

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
        hostAliases = new String(obj.hostAliases);
        shareOneBatchJob = obj.shareOneBatchJob;

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
                (hostAliases == obj.hostAliases) &&
                (shareOneBatchJob == obj.shareOneBatchJob));
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

    public void SetHostAliases(String hostAliases_)
    {
        hostAliases = hostAliases_;
        Select(22);
    }

    public void SetShareOneBatchJob(boolean shareOneBatchJob_)
    {
        shareOneBatchJob = shareOneBatchJob_;
        Select(23);
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
    public String  GetHostAliases() { return hostAliases; }
    public boolean GetShareOneBatchJob() { return shareOneBatchJob; }

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
            buf.WriteString(hostAliases);
        if(WriteSelect(23, buf))
            buf.WriteBool(shareOneBatchJob);
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
                SetHostAliases(buf.ReadString());
                break;
            case 23:
                SetShareOneBatchJob(buf.ReadBool());
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
    private String  hostAliases;
    private boolean shareOneBatchJob;
}

