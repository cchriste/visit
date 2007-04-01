package llnl.visit;

import java.util.Vector;
import java.lang.Integer;

// ****************************************************************************
// Class: ViewerRPC
//
// Purpose:
//    This class contains the attributes for controlling the viewer.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Mon Mar 15 14:48:00 PST 2004
//
// Modifications:
//   
// ****************************************************************************

public class ViewerRPC extends AttributeSubject
{
    // Constants
    public final static int VIEWERRPCTYPE_CLOSERPC = 0;
    public final static int VIEWERRPCTYPE_ADDWINDOWRPC = 1;
    public final static int VIEWERRPCTYPE_DELETEWINDOWRPC = 2;
    public final static int VIEWERRPCTYPE_SETWINDOWLAYOUTRPC = 3;
    public final static int VIEWERRPCTYPE_SETACTIVEWINDOWRPC = 4;
    public final static int VIEWERRPCTYPE_CLEARWINDOWRPC = 5;
    public final static int VIEWERRPCTYPE_CLEARALLWINDOWSRPC = 6;
    public final static int VIEWERRPCTYPE_OPENDATABASERPC = 7;
    public final static int VIEWERRPCTYPE_REOPENDATABASERPC = 8;
    public final static int VIEWERRPCTYPE_REPLACEDATABASERPC = 9;
    public final static int VIEWERRPCTYPE_OVERLAYDATABASERPC = 10;
    public final static int VIEWERRPCTYPE_OPENCOMPUTEENGINERPC = 11;
    public final static int VIEWERRPCTYPE_CLOSECOMPUTEENGINERPC = 12;
    public final static int VIEWERRPCTYPE_ANIMATIONSETNFRAMESRPC = 13;
    public final static int VIEWERRPCTYPE_ANIMATIONPLAYRPC = 14;
    public final static int VIEWERRPCTYPE_ANIMATIONREVERSEPLAYRPC = 15;
    public final static int VIEWERRPCTYPE_ANIMATIONSTOPRPC = 16;
    public final static int VIEWERRPCTYPE_ANIMATIONNEXTFRAMERPC = 17;
    public final static int VIEWERRPCTYPE_ANIMATIONPREVIOUSFRAMERPC = 18;
    public final static int VIEWERRPCTYPE_ANIMATIONSETFRAMERPC = 19;
    public final static int VIEWERRPCTYPE_ADDPLOTRPC = 20;
    public final static int VIEWERRPCTYPE_SETPLOTFRAMERANGERPC = 21;
    public final static int VIEWERRPCTYPE_DELETEPLOTKEYFRAMERPC = 22;
    public final static int VIEWERRPCTYPE_MOVEPLOTKEYFRAMERPC = 23;
    public final static int VIEWERRPCTYPE_DELETEACTIVEPLOTSRPC = 24;
    public final static int VIEWERRPCTYPE_HIDEACTIVEPLOTSRPC = 25;
    public final static int VIEWERRPCTYPE_DRAWPLOTSRPC = 26;
    public final static int VIEWERRPCTYPE_DISABLEREDRAWRPC = 27;
    public final static int VIEWERRPCTYPE_REDRAWRPC = 28;
    public final static int VIEWERRPCTYPE_SETACTIVEPLOTSRPC = 29;
    public final static int VIEWERRPCTYPE_CHANGEACTIVEPLOTSVARRPC = 30;
    public final static int VIEWERRPCTYPE_ADDOPERATORRPC = 31;
    public final static int VIEWERRPCTYPE_PROMOTEOPERATORRPC = 32;
    public final static int VIEWERRPCTYPE_DEMOTEOPERATORRPC = 33;
    public final static int VIEWERRPCTYPE_REMOVEOPERATORRPC = 34;
    public final static int VIEWERRPCTYPE_REMOVELASTOPERATORRPC = 35;
    public final static int VIEWERRPCTYPE_REMOVEALLOPERATORSRPC = 36;
    public final static int VIEWERRPCTYPE_SAVEWINDOWRPC = 37;
    public final static int VIEWERRPCTYPE_SETDEFAULTPLOTOPTIONSRPC = 38;
    public final static int VIEWERRPCTYPE_SETPLOTOPTIONSRPC = 39;
    public final static int VIEWERRPCTYPE_SETDEFAULTOPERATOROPTIONSRPC = 40;
    public final static int VIEWERRPCTYPE_SETOPERATOROPTIONSRPC = 41;
    public final static int VIEWERRPCTYPE_WRITECONFIGFILERPC = 42;
    public final static int VIEWERRPCTYPE_CONNECTTOMETADATASERVERRPC = 43;
    public final static int VIEWERRPCTYPE_ICONIFYALLWINDOWSRPC = 44;
    public final static int VIEWERRPCTYPE_DEICONIFYALLWINDOWSRPC = 45;
    public final static int VIEWERRPCTYPE_SHOWALLWINDOWSRPC = 46;
    public final static int VIEWERRPCTYPE_HIDEALLWINDOWSRPC = 47;
    public final static int VIEWERRPCTYPE_UPDATECOLORTABLERPC = 48;
    public final static int VIEWERRPCTYPE_SETANNOTATIONATTRIBUTESRPC = 49;
    public final static int VIEWERRPCTYPE_SETDEFAULTANNOTATIONATTRIBUTESRPC = 50;
    public final static int VIEWERRPCTYPE_RESETANNOTATIONATTRIBUTESRPC = 51;
    public final static int VIEWERRPCTYPE_SETKEYFRAMEATTRIBUTESRPC = 52;
    public final static int VIEWERRPCTYPE_SETPLOTSILRESTRICTIONRPC = 53;
    public final static int VIEWERRPCTYPE_SETVIEWCURVERPC = 54;
    public final static int VIEWERRPCTYPE_SETVIEW2DRPC = 55;
    public final static int VIEWERRPCTYPE_SETVIEW3DRPC = 56;
    public final static int VIEWERRPCTYPE_RESETPLOTOPTIONSRPC = 57;
    public final static int VIEWERRPCTYPE_RESETOPERATOROPTIONSRPC = 58;
    public final static int VIEWERRPCTYPE_SETAPPEARANCERPC = 59;
    public final static int VIEWERRPCTYPE_PROCESSEXPRESSIONSRPC = 60;
    public final static int VIEWERRPCTYPE_SETLIGHTLISTRPC = 61;
    public final static int VIEWERRPCTYPE_SETDEFAULTLIGHTLISTRPC = 62;
    public final static int VIEWERRPCTYPE_RESETLIGHTLISTRPC = 63;
    public final static int VIEWERRPCTYPE_SETANIMATIONATTRIBUTESRPC = 64;
    public final static int VIEWERRPCTYPE_SETWINDOWAREARPC = 65;
    public final static int VIEWERRPCTYPE_PRINTWINDOWRPC = 66;
    public final static int VIEWERRPCTYPE_RESETVIEWRPC = 67;
    public final static int VIEWERRPCTYPE_RECENTERVIEWRPC = 68;
    public final static int VIEWERRPCTYPE_TOGGLEMAINTAINVIEWMODERPC = 69;
    public final static int VIEWERRPCTYPE_TOGGLEBOUNDINGBOXMODERPC = 70;
    public final static int VIEWERRPCTYPE_TOGGLECAMERAVIEWMODERPC = 71;
    public final static int VIEWERRPCTYPE_TOGGLEPERSPECTIVEVIEWRPC = 72;
    public final static int VIEWERRPCTYPE_TOGGLESPINMODERPC = 73;
    public final static int VIEWERRPCTYPE_TOGGLELOCKTIMERPC = 74;
    public final static int VIEWERRPCTYPE_TOGGLELOCKTOOLSRPC = 75;
    public final static int VIEWERRPCTYPE_TOGGLELOCKVIEWMODERPC = 76;
    public final static int VIEWERRPCTYPE_TOGGLEFULLFRAMERPC = 77;
    public final static int VIEWERRPCTYPE_UNDOVIEWRPC = 78;
    public final static int VIEWERRPCTYPE_INVERTBACKGROUNDRPC = 79;
    public final static int VIEWERRPCTYPE_CLEARPICKPOINTSRPC = 80;
    public final static int VIEWERRPCTYPE_SETWINDOWMODERPC = 81;
    public final static int VIEWERRPCTYPE_ENABLETOOLRPC = 82;
    public final static int VIEWERRPCTYPE_COPYVIEWTOWINDOWRPC = 83;
    public final static int VIEWERRPCTYPE_COPYLIGHTINGTOWINDOWRPC = 84;
    public final static int VIEWERRPCTYPE_COPYANNOTATIONSTOWINDOWRPC = 85;
    public final static int VIEWERRPCTYPE_COPYPLOTSTOWINDOWRPC = 86;
    public final static int VIEWERRPCTYPE_CLEARCACHERPC = 87;
    public final static int VIEWERRPCTYPE_CLEARCACHEFORALLENGINESRPC = 88;
    public final static int VIEWERRPCTYPE_SETVIEWEXTENTSTYPERPC = 89;
    public final static int VIEWERRPCTYPE_CLEARREFLINESRPC = 90;
    public final static int VIEWERRPCTYPE_SETRENDERINGATTRIBUTESRPC = 91;
    public final static int VIEWERRPCTYPE_DATABASEQUERYRPC = 92;
    public final static int VIEWERRPCTYPE_POINTQUERYRPC = 93;
    public final static int VIEWERRPCTYPE_LINEQUERYRPC = 94;
    public final static int VIEWERRPCTYPE_CLONEWINDOWRPC = 95;
    public final static int VIEWERRPCTYPE_SETMATERIALATTRIBUTESRPC = 96;
    public final static int VIEWERRPCTYPE_SETDEFAULTMATERIALATTRIBUTESRPC = 97;
    public final static int VIEWERRPCTYPE_RESETMATERIALATTRIBUTESRPC = 98;
    public final static int VIEWERRPCTYPE_SETPLOTDATABASESTATERPC = 99;
    public final static int VIEWERRPCTYPE_DELETEPLOTDATABASEKEYFRAMERPC = 100;
    public final static int VIEWERRPCTYPE_MOVEPLOTDATABASEKEYFRAMERPC = 101;
    public final static int VIEWERRPCTYPE_CLEARVIEWKEYFRAMESRPC = 102;
    public final static int VIEWERRPCTYPE_DELETEVIEWKEYFRAMERPC = 103;
    public final static int VIEWERRPCTYPE_MOVEVIEWKEYFRAMERPC = 104;
    public final static int VIEWERRPCTYPE_SETVIEWKEYFRAMERPC = 105;
    public final static int VIEWERRPCTYPE_OPENMDSERVERRPC = 106;
    public final static int VIEWERRPCTYPE_ENABLETOOLBARRPC = 107;
    public final static int VIEWERRPCTYPE_HIDETOOLBARSRPC = 108;
    public final static int VIEWERRPCTYPE_HIDETOOLBARSFORALLWINDOWSRPC = 109;
    public final static int VIEWERRPCTYPE_SHOWTOOLBARSRPC = 110;
    public final static int VIEWERRPCTYPE_SHOWTOOLBARSFORALLWINDOWSRPC = 111;
    public final static int VIEWERRPCTYPE_SETTOOLBARICONSIZERPC = 112;
    public final static int VIEWERRPCTYPE_SAVEVIEWRPC = 113;
    public final static int VIEWERRPCTYPE_SETGLOBALLINEOUTATTRIBUTESRPC = 114;
    public final static int VIEWERRPCTYPE_SETPICKATTRIBUTESRPC = 115;
    public final static int VIEWERRPCTYPE_EXPORTCOLORTABLERPC = 116;
    public final static int VIEWERRPCTYPE_EXPORTENTIRESTATERPC = 117;
    public final static int VIEWERRPCTYPE_IMPORTENTIRESTATERPC = 118;
    public final static int VIEWERRPCTYPE_RESETPICKATTRIBUTESRPC = 119;
    public final static int VIEWERRPCTYPE_ADDANNOTATIONOBJECTRPC = 120;
    public final static int VIEWERRPCTYPE_HIDEACTIVEANNOTATIONOBJECTSRPC = 121;
    public final static int VIEWERRPCTYPE_DELETEACTIVEANNOTATIONOBJECTSRPC = 122;
    public final static int VIEWERRPCTYPE_RAISEACTIVEANNOTATIONOBJECTSRPC = 123;
    public final static int VIEWERRPCTYPE_LOWERACTIVEANNOTATIONOBJECTSRPC = 124;
    public final static int VIEWERRPCTYPE_SETANNOTATIONOBJECTOPTIONSRPC = 125;
    public final static int VIEWERRPCTYPE_SETDEFAULTANNOTATIONOBJECTLISTRPC = 126;
    public final static int VIEWERRPCTYPE_RESETANNOTATIONOBJECTLISTRPC = 127;
    public final static int VIEWERRPCTYPE_RESETPICKLETTERRPC = 128;
    public final static int VIEWERRPCTYPE_SETDEFAULTPICKATTRIBUTESRPC = 129;
    public final static int VIEWERRPCTYPE_CHOOSECENTEROFROTATIONRPC = 130;
    public final static int VIEWERRPCTYPE_SETCENTEROFROTATIONRPC = 131;
    public final static int VIEWERRPCTYPE_MAXRPC = 132;


    public ViewerRPC()
    {
        super(28);

        RPCType = VIEWERRPCTYPE_CLOSERPC;
        windowLayout = 1;
        windowId = 0;
        windowMode = 0;
        windowArea = new String("(null)");
        database = new String("(null)");
        programHost = new String("(null)");
        programOptions = new Vector();
        nFrames = 0;
        frameNumber = 0;
        frameRange = new int[2];
        frameRange[0] = 0;
        frameRange[1] = 0;
        frame = 0;
        plotType = 0;
        operatorType = 0;
        variable = new String("(null)");
        activePlotIds = new Vector();
        activeOperatorIds = new Vector();
        expandedPlotIds = new Vector();
        colorTableName = new String("(null)");
        queryName = new String("(null)");
        queryPoint1 = new double[3];
        queryPoint1[0] = 0;
        queryPoint1[1] = 0;
        queryPoint1[2] = 0;
        queryPoint2 = new double[3];
        queryPoint2[0] = 0;
        queryPoint2[1] = 0;
        queryPoint2[2] = 0;
        queryVariables = new Vector();
        toolId = 0;
        boolFlag = false;
        intArg1 = 0;
        intArg2 = 0;
        intArg3 = 0;
    }

    public ViewerRPC(ViewerRPC obj)
    {
        super(28);

        int i;

        RPCType = obj.RPCType;
        windowLayout = obj.windowLayout;
        windowId = obj.windowId;
        windowMode = obj.windowMode;
        windowArea = new String(obj.windowArea);
        database = new String(obj.database);
        programHost = new String(obj.programHost);
        programOptions = new Vector(obj.programOptions.size());
        for(i = 0; i < obj.programOptions.size(); ++i)
            programOptions.addElement(new String((String)obj.programOptions.elementAt(i)));

        nFrames = obj.nFrames;
        frameNumber = obj.frameNumber;
        frameRange = new int[2];
        frameRange[0] = obj.frameRange[0];
        frameRange[1] = obj.frameRange[1];

        frame = obj.frame;
        plotType = obj.plotType;
        operatorType = obj.operatorType;
        variable = new String(obj.variable);
        activePlotIds = new Vector();
        for(i = 0; i < obj.activePlotIds.size(); ++i)
        {
            Integer iv = (Integer)obj.activePlotIds.elementAt(i);
            activePlotIds.addElement(new Integer(iv.intValue()));
        }
        activeOperatorIds = new Vector();
        for(i = 0; i < obj.activeOperatorIds.size(); ++i)
        {
            Integer iv = (Integer)obj.activeOperatorIds.elementAt(i);
            activeOperatorIds.addElement(new Integer(iv.intValue()));
        }
        expandedPlotIds = new Vector();
        for(i = 0; i < obj.expandedPlotIds.size(); ++i)
        {
            Integer iv = (Integer)obj.expandedPlotIds.elementAt(i);
            expandedPlotIds.addElement(new Integer(iv.intValue()));
        }
        colorTableName = new String(obj.colorTableName);
        queryName = new String(obj.queryName);
        queryPoint1 = new double[3];
        queryPoint1[0] = obj.queryPoint1[0];
        queryPoint1[1] = obj.queryPoint1[1];
        queryPoint1[2] = obj.queryPoint1[2];

        queryPoint2 = new double[3];
        queryPoint2[0] = obj.queryPoint2[0];
        queryPoint2[1] = obj.queryPoint2[1];
        queryPoint2[2] = obj.queryPoint2[2];

        queryVariables = new Vector(obj.queryVariables.size());
        for(i = 0; i < obj.queryVariables.size(); ++i)
            queryVariables.addElement(new String((String)obj.queryVariables.elementAt(i)));

        toolId = obj.toolId;
        boolFlag = obj.boolFlag;
        intArg1 = obj.intArg1;
        intArg2 = obj.intArg2;
        intArg3 = obj.intArg3;

        SelectAll();
    }

    public boolean equals(ViewerRPC obj)
    {
        int i;

        // Compare the frameRange arrays.
        boolean frameRange_equal = true;
        for(i = 0; i < 2 && frameRange_equal; ++i)
            frameRange_equal = (frameRange[i] == obj.frameRange[i]);

        // Compare the queryPoint1 arrays.
        boolean queryPoint1_equal = true;
        for(i = 0; i < 3 && queryPoint1_equal; ++i)
            queryPoint1_equal = (queryPoint1[i] == obj.queryPoint1[i]);

        // Compare the queryPoint2 arrays.
        boolean queryPoint2_equal = true;
        for(i = 0; i < 3 && queryPoint2_equal; ++i)
            queryPoint2_equal = (queryPoint2[i] == obj.queryPoint2[i]);

        // Create the return value
        return ((RPCType == obj.RPCType) &&
                (windowLayout == obj.windowLayout) &&
                (windowId == obj.windowId) &&
                (windowMode == obj.windowMode) &&
                (windowArea == obj.windowArea) &&
                (database == obj.database) &&
                (programHost == obj.programHost) &&
                (programOptions == obj.programOptions) &&
                (nFrames == obj.nFrames) &&
                (frameNumber == obj.frameNumber) &&
                frameRange_equal &&
                (frame == obj.frame) &&
                (plotType == obj.plotType) &&
                (operatorType == obj.operatorType) &&
                (variable == obj.variable) &&
                (activePlotIds == obj.activePlotIds) &&
                (activeOperatorIds == obj.activeOperatorIds) &&
                (expandedPlotIds == obj.expandedPlotIds) &&
                (colorTableName == obj.colorTableName) &&
                (queryName == obj.queryName) &&
                queryPoint1_equal &&
                queryPoint2_equal &&
                (queryVariables == obj.queryVariables) &&
                (toolId == obj.toolId) &&
                (boolFlag == obj.boolFlag) &&
                (intArg1 == obj.intArg1) &&
                (intArg2 == obj.intArg2) &&
                (intArg3 == obj.intArg3));
    }

    // Property setting methods
    public void SetRPCType(int RPCType_)
    {
        RPCType = RPCType_;
        Select(0);
    }

    public void SetWindowLayout(int windowLayout_)
    {
        windowLayout = windowLayout_;
        Select(1);
    }

    public void SetWindowId(int windowId_)
    {
        windowId = windowId_;
        Select(2);
    }

    public void SetWindowMode(int windowMode_)
    {
        windowMode = windowMode_;
        Select(3);
    }

    public void SetWindowArea(String windowArea_)
    {
        windowArea = windowArea_;
        Select(4);
    }

    public void SetDatabase(String database_)
    {
        database = database_;
        Select(5);
    }

    public void SetProgramHost(String programHost_)
    {
        programHost = programHost_;
        Select(6);
    }

    public void SetProgramOptions(Vector programOptions_)
    {
        programOptions = programOptions_;
        Select(7);
    }

    public void SetNFrames(int nFrames_)
    {
        nFrames = nFrames_;
        Select(8);
    }

    public void SetFrameNumber(int frameNumber_)
    {
        frameNumber = frameNumber_;
        Select(9);
    }

    public void SetFrameRange(int[] frameRange_)
    {
        frameRange[0] = frameRange_[0];
        frameRange[1] = frameRange_[1];
        Select(10);
    }

    public void SetFrameRange(int e0, int e1)
    {
        frameRange[0] = e0;
        frameRange[1] = e1;
        Select(10);
    }

    public void SetFrame(int frame_)
    {
        frame = frame_;
        Select(11);
    }

    public void SetPlotType(int plotType_)
    {
        plotType = plotType_;
        Select(12);
    }

    public void SetOperatorType(int operatorType_)
    {
        operatorType = operatorType_;
        Select(13);
    }

    public void SetVariable(String variable_)
    {
        variable = variable_;
        Select(14);
    }

    public void SetActivePlotIds(Vector activePlotIds_)
    {
        activePlotIds = activePlotIds_;
        Select(15);
    }

    public void SetActiveOperatorIds(Vector activeOperatorIds_)
    {
        activeOperatorIds = activeOperatorIds_;
        Select(16);
    }

    public void SetExpandedPlotIds(Vector expandedPlotIds_)
    {
        expandedPlotIds = expandedPlotIds_;
        Select(17);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(18);
    }

    public void SetQueryName(String queryName_)
    {
        queryName = queryName_;
        Select(19);
    }

    public void SetQueryPoint1(double[] queryPoint1_)
    {
        queryPoint1[0] = queryPoint1_[0];
        queryPoint1[1] = queryPoint1_[1];
        queryPoint1[2] = queryPoint1_[2];
        Select(20);
    }

    public void SetQueryPoint1(double e0, double e1, double e2)
    {
        queryPoint1[0] = e0;
        queryPoint1[1] = e1;
        queryPoint1[2] = e2;
        Select(20);
    }

    public void SetQueryPoint2(double[] queryPoint2_)
    {
        queryPoint2[0] = queryPoint2_[0];
        queryPoint2[1] = queryPoint2_[1];
        queryPoint2[2] = queryPoint2_[2];
        Select(21);
    }

    public void SetQueryPoint2(double e0, double e1, double e2)
    {
        queryPoint2[0] = e0;
        queryPoint2[1] = e1;
        queryPoint2[2] = e2;
        Select(21);
    }

    public void SetQueryVariables(Vector queryVariables_)
    {
        queryVariables = queryVariables_;
        Select(22);
    }

    public void SetToolId(int toolId_)
    {
        toolId = toolId_;
        Select(23);
    }

    public void SetBoolFlag(boolean boolFlag_)
    {
        boolFlag = boolFlag_;
        Select(24);
    }

    public void SetIntArg1(int intArg1_)
    {
        intArg1 = intArg1_;
        Select(25);
    }

    public void SetIntArg2(int intArg2_)
    {
        intArg2 = intArg2_;
        Select(26);
    }

    public void SetIntArg3(int intArg3_)
    {
        intArg3 = intArg3_;
        Select(27);
    }

    // Property getting methods
    public int      GetRPCType() { return RPCType; }
    public int      GetWindowLayout() { return windowLayout; }
    public int      GetWindowId() { return windowId; }
    public int      GetWindowMode() { return windowMode; }
    public String   GetWindowArea() { return windowArea; }
    public String   GetDatabase() { return database; }
    public String   GetProgramHost() { return programHost; }
    public Vector   GetProgramOptions() { return programOptions; }
    public int      GetNFrames() { return nFrames; }
    public int      GetFrameNumber() { return frameNumber; }
    public int[]    GetFrameRange() { return frameRange; }
    public int      GetFrame() { return frame; }
    public int      GetPlotType() { return plotType; }
    public int      GetOperatorType() { return operatorType; }
    public String   GetVariable() { return variable; }
    public Vector   GetActivePlotIds() { return activePlotIds; }
    public Vector   GetActiveOperatorIds() { return activeOperatorIds; }
    public Vector   GetExpandedPlotIds() { return expandedPlotIds; }
    public String   GetColorTableName() { return colorTableName; }
    public String   GetQueryName() { return queryName; }
    public double[] GetQueryPoint1() { return queryPoint1; }
    public double[] GetQueryPoint2() { return queryPoint2; }
    public Vector   GetQueryVariables() { return queryVariables; }
    public int      GetToolId() { return toolId; }
    public boolean  GetBoolFlag() { return boolFlag; }
    public int      GetIntArg1() { return intArg1; }
    public int      GetIntArg2() { return intArg2; }
    public int      GetIntArg3() { return intArg3; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(RPCType);
        if(WriteSelect(1, buf))
            buf.WriteInt(windowLayout);
        if(WriteSelect(2, buf))
            buf.WriteInt(windowId);
        if(WriteSelect(3, buf))
            buf.WriteInt(windowMode);
        if(WriteSelect(4, buf))
            buf.WriteString(windowArea);
        if(WriteSelect(5, buf))
            buf.WriteString(database);
        if(WriteSelect(6, buf))
            buf.WriteString(programHost);
        if(WriteSelect(7, buf))
            buf.WriteStringVector(programOptions);
        if(WriteSelect(8, buf))
            buf.WriteInt(nFrames);
        if(WriteSelect(9, buf))
            buf.WriteInt(frameNumber);
        if(WriteSelect(10, buf))
            buf.WriteIntArray(frameRange);
        if(WriteSelect(11, buf))
            buf.WriteInt(frame);
        if(WriteSelect(12, buf))
            buf.WriteInt(plotType);
        if(WriteSelect(13, buf))
            buf.WriteInt(operatorType);
        if(WriteSelect(14, buf))
            buf.WriteString(variable);
        if(WriteSelect(15, buf))
            buf.WriteIntVector(activePlotIds);
        if(WriteSelect(16, buf))
            buf.WriteIntVector(activeOperatorIds);
        if(WriteSelect(17, buf))
            buf.WriteIntVector(expandedPlotIds);
        if(WriteSelect(18, buf))
            buf.WriteString(colorTableName);
        if(WriteSelect(19, buf))
            buf.WriteString(queryName);
        if(WriteSelect(20, buf))
            buf.WriteDoubleArray(queryPoint1);
        if(WriteSelect(21, buf))
            buf.WriteDoubleArray(queryPoint2);
        if(WriteSelect(22, buf))
            buf.WriteStringVector(queryVariables);
        if(WriteSelect(23, buf))
            buf.WriteInt(toolId);
        if(WriteSelect(24, buf))
            buf.WriteBool(boolFlag);
        if(WriteSelect(25, buf))
            buf.WriteInt(intArg1);
        if(WriteSelect(26, buf))
            buf.WriteInt(intArg2);
        if(WriteSelect(27, buf))
            buf.WriteInt(intArg3);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetRPCType(buf.ReadInt());
                break;
            case 1:
                SetWindowLayout(buf.ReadInt());
                break;
            case 2:
                SetWindowId(buf.ReadInt());
                break;
            case 3:
                SetWindowMode(buf.ReadInt());
                break;
            case 4:
                SetWindowArea(buf.ReadString());
                break;
            case 5:
                SetDatabase(buf.ReadString());
                break;
            case 6:
                SetProgramHost(buf.ReadString());
                break;
            case 7:
                SetProgramOptions(buf.ReadStringVector());
                break;
            case 8:
                SetNFrames(buf.ReadInt());
                break;
            case 9:
                SetFrameNumber(buf.ReadInt());
                break;
            case 10:
                SetFrameRange(buf.ReadIntArray());
                break;
            case 11:
                SetFrame(buf.ReadInt());
                break;
            case 12:
                SetPlotType(buf.ReadInt());
                break;
            case 13:
                SetOperatorType(buf.ReadInt());
                break;
            case 14:
                SetVariable(buf.ReadString());
                break;
            case 15:
                SetActivePlotIds(buf.ReadIntVector());
                break;
            case 16:
                SetActiveOperatorIds(buf.ReadIntVector());
                break;
            case 17:
                SetExpandedPlotIds(buf.ReadIntVector());
                break;
            case 18:
                SetColorTableName(buf.ReadString());
                break;
            case 19:
                SetQueryName(buf.ReadString());
                break;
            case 20:
                SetQueryPoint1(buf.ReadDoubleArray());
                break;
            case 21:
                SetQueryPoint2(buf.ReadDoubleArray());
                break;
            case 22:
                SetQueryVariables(buf.ReadStringVector());
                break;
            case 23:
                SetToolId(buf.ReadInt());
                break;
            case 24:
                SetBoolFlag(buf.ReadBool());
                break;
            case 25:
                SetIntArg1(buf.ReadInt());
                break;
            case 26:
                SetIntArg2(buf.ReadInt());
                break;
            case 27:
                SetIntArg3(buf.ReadInt());
                break;
            }
        }
    }


    // Attributes
    private int      RPCType;
    private int      windowLayout;
    private int      windowId;
    private int      windowMode;
    private String   windowArea;
    private String   database;
    private String   programHost;
    private Vector   programOptions; // vector of String objects
    private int      nFrames;
    private int      frameNumber;
    private int[]    frameRange;
    private int      frame;
    private int      plotType;
    private int      operatorType;
    private String   variable;
    private Vector   activePlotIds; // vector of Integer objects
    private Vector   activeOperatorIds; // vector of Integer objects
    private Vector   expandedPlotIds; // vector of Integer objects
    private String   colorTableName;
    private String   queryName;
    private double[] queryPoint1;
    private double[] queryPoint2;
    private Vector   queryVariables; // vector of String objects
    private int      toolId;
    private boolean  boolFlag;
    private int      intArg1;
    private int      intArg2;
    private int      intArg3;
}

