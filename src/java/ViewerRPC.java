// ***************************************************************************
//
// Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400142
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
import java.lang.Integer;
import java.lang.Double;

// ****************************************************************************
// Class: ViewerRPC
//
// Purpose:
//    This class contains the attributes for controlling the viewer.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Fri Dec 14 14:30:55 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class ViewerRPC extends AttributeSubject
{
    // Enum values
    public final static int VIEWERRPCTYPE_CLOSERPC = 0;
    public final static int VIEWERRPCTYPE_DETACHRPC = 1;
    public final static int VIEWERRPCTYPE_ADDWINDOWRPC = 2;
    public final static int VIEWERRPCTYPE_DELETEWINDOWRPC = 3;
    public final static int VIEWERRPCTYPE_SETWINDOWLAYOUTRPC = 4;
    public final static int VIEWERRPCTYPE_SETACTIVEWINDOWRPC = 5;
    public final static int VIEWERRPCTYPE_CLEARWINDOWRPC = 6;
    public final static int VIEWERRPCTYPE_CLEARALLWINDOWSRPC = 7;
    public final static int VIEWERRPCTYPE_OPENDATABASERPC = 8;
    public final static int VIEWERRPCTYPE_CLOSEDATABASERPC = 9;
    public final static int VIEWERRPCTYPE_ACTIVATEDATABASERPC = 10;
    public final static int VIEWERRPCTYPE_CHECKFORNEWSTATESRPC = 11;
    public final static int VIEWERRPCTYPE_CREATEDATABASECORRELATIONRPC = 12;
    public final static int VIEWERRPCTYPE_ALTERDATABASECORRELATIONRPC = 13;
    public final static int VIEWERRPCTYPE_DELETEDATABASECORRELATIONRPC = 14;
    public final static int VIEWERRPCTYPE_REOPENDATABASERPC = 15;
    public final static int VIEWERRPCTYPE_REPLACEDATABASERPC = 16;
    public final static int VIEWERRPCTYPE_OVERLAYDATABASERPC = 17;
    public final static int VIEWERRPCTYPE_OPENCOMPUTEENGINERPC = 18;
    public final static int VIEWERRPCTYPE_CLOSECOMPUTEENGINERPC = 19;
    public final static int VIEWERRPCTYPE_ANIMATIONSETNFRAMESRPC = 20;
    public final static int VIEWERRPCTYPE_ANIMATIONPLAYRPC = 21;
    public final static int VIEWERRPCTYPE_ANIMATIONREVERSEPLAYRPC = 22;
    public final static int VIEWERRPCTYPE_ANIMATIONSTOPRPC = 23;
    public final static int VIEWERRPCTYPE_TIMESLIDERNEXTSTATERPC = 24;
    public final static int VIEWERRPCTYPE_TIMESLIDERPREVIOUSSTATERPC = 25;
    public final static int VIEWERRPCTYPE_SETTIMESLIDERSTATERPC = 26;
    public final static int VIEWERRPCTYPE_SETACTIVETIMESLIDERRPC = 27;
    public final static int VIEWERRPCTYPE_ADDPLOTRPC = 28;
    public final static int VIEWERRPCTYPE_SETPLOTFRAMERANGERPC = 29;
    public final static int VIEWERRPCTYPE_DELETEPLOTKEYFRAMERPC = 30;
    public final static int VIEWERRPCTYPE_MOVEPLOTKEYFRAMERPC = 31;
    public final static int VIEWERRPCTYPE_DELETEACTIVEPLOTSRPC = 32;
    public final static int VIEWERRPCTYPE_HIDEACTIVEPLOTSRPC = 33;
    public final static int VIEWERRPCTYPE_DRAWPLOTSRPC = 34;
    public final static int VIEWERRPCTYPE_DISABLEREDRAWRPC = 35;
    public final static int VIEWERRPCTYPE_REDRAWRPC = 36;
    public final static int VIEWERRPCTYPE_SETACTIVEPLOTSRPC = 37;
    public final static int VIEWERRPCTYPE_CHANGEACTIVEPLOTSVARRPC = 38;
    public final static int VIEWERRPCTYPE_ADDOPERATORRPC = 39;
    public final static int VIEWERRPCTYPE_ADDINITIALIZEDOPERATORRPC = 40;
    public final static int VIEWERRPCTYPE_PROMOTEOPERATORRPC = 41;
    public final static int VIEWERRPCTYPE_DEMOTEOPERATORRPC = 42;
    public final static int VIEWERRPCTYPE_REMOVEOPERATORRPC = 43;
    public final static int VIEWERRPCTYPE_REMOVELASTOPERATORRPC = 44;
    public final static int VIEWERRPCTYPE_REMOVEALLOPERATORSRPC = 45;
    public final static int VIEWERRPCTYPE_SAVEWINDOWRPC = 46;
    public final static int VIEWERRPCTYPE_SETDEFAULTPLOTOPTIONSRPC = 47;
    public final static int VIEWERRPCTYPE_SETPLOTOPTIONSRPC = 48;
    public final static int VIEWERRPCTYPE_SETDEFAULTOPERATOROPTIONSRPC = 49;
    public final static int VIEWERRPCTYPE_SETOPERATOROPTIONSRPC = 50;
    public final static int VIEWERRPCTYPE_WRITECONFIGFILERPC = 51;
    public final static int VIEWERRPCTYPE_CONNECTTOMETADATASERVERRPC = 52;
    public final static int VIEWERRPCTYPE_ICONIFYALLWINDOWSRPC = 53;
    public final static int VIEWERRPCTYPE_DEICONIFYALLWINDOWSRPC = 54;
    public final static int VIEWERRPCTYPE_SHOWALLWINDOWSRPC = 55;
    public final static int VIEWERRPCTYPE_HIDEALLWINDOWSRPC = 56;
    public final static int VIEWERRPCTYPE_UPDATECOLORTABLERPC = 57;
    public final static int VIEWERRPCTYPE_SETANNOTATIONATTRIBUTESRPC = 58;
    public final static int VIEWERRPCTYPE_SETDEFAULTANNOTATIONATTRIBUTESRPC = 59;
    public final static int VIEWERRPCTYPE_RESETANNOTATIONATTRIBUTESRPC = 60;
    public final static int VIEWERRPCTYPE_SETKEYFRAMEATTRIBUTESRPC = 61;
    public final static int VIEWERRPCTYPE_SETPLOTSILRESTRICTIONRPC = 62;
    public final static int VIEWERRPCTYPE_SETVIEWCURVERPC = 63;
    public final static int VIEWERRPCTYPE_SETVIEW2DRPC = 64;
    public final static int VIEWERRPCTYPE_SETVIEW3DRPC = 65;
    public final static int VIEWERRPCTYPE_RESETPLOTOPTIONSRPC = 66;
    public final static int VIEWERRPCTYPE_RESETOPERATOROPTIONSRPC = 67;
    public final static int VIEWERRPCTYPE_SETAPPEARANCERPC = 68;
    public final static int VIEWERRPCTYPE_PROCESSEXPRESSIONSRPC = 69;
    public final static int VIEWERRPCTYPE_SETLIGHTLISTRPC = 70;
    public final static int VIEWERRPCTYPE_SETDEFAULTLIGHTLISTRPC = 71;
    public final static int VIEWERRPCTYPE_RESETLIGHTLISTRPC = 72;
    public final static int VIEWERRPCTYPE_SETANIMATIONATTRIBUTESRPC = 73;
    public final static int VIEWERRPCTYPE_SETWINDOWAREARPC = 74;
    public final static int VIEWERRPCTYPE_PRINTWINDOWRPC = 75;
    public final static int VIEWERRPCTYPE_RESETVIEWRPC = 76;
    public final static int VIEWERRPCTYPE_RECENTERVIEWRPC = 77;
    public final static int VIEWERRPCTYPE_TOGGLEMAINTAINVIEWMODERPC = 78;
    public final static int VIEWERRPCTYPE_TOGGLEMAINTAINDATAMODERPC = 79;
    public final static int VIEWERRPCTYPE_TOGGLEBOUNDINGBOXMODERPC = 80;
    public final static int VIEWERRPCTYPE_TOGGLECAMERAVIEWMODERPC = 81;
    public final static int VIEWERRPCTYPE_TOGGLEPERSPECTIVEVIEWRPC = 82;
    public final static int VIEWERRPCTYPE_TOGGLESPINMODERPC = 83;
    public final static int VIEWERRPCTYPE_TOGGLELOCKTIMERPC = 84;
    public final static int VIEWERRPCTYPE_TOGGLELOCKTOOLSRPC = 85;
    public final static int VIEWERRPCTYPE_TOGGLELOCKVIEWMODERPC = 86;
    public final static int VIEWERRPCTYPE_TOGGLEFULLFRAMERPC = 87;
    public final static int VIEWERRPCTYPE_UNDOVIEWRPC = 88;
    public final static int VIEWERRPCTYPE_REDOVIEWRPC = 89;
    public final static int VIEWERRPCTYPE_INVERTBACKGROUNDRPC = 90;
    public final static int VIEWERRPCTYPE_CLEARPICKPOINTSRPC = 91;
    public final static int VIEWERRPCTYPE_SETWINDOWMODERPC = 92;
    public final static int VIEWERRPCTYPE_ENABLETOOLRPC = 93;
    public final static int VIEWERRPCTYPE_COPYVIEWTOWINDOWRPC = 94;
    public final static int VIEWERRPCTYPE_COPYLIGHTINGTOWINDOWRPC = 95;
    public final static int VIEWERRPCTYPE_COPYANNOTATIONSTOWINDOWRPC = 96;
    public final static int VIEWERRPCTYPE_COPYPLOTSTOWINDOWRPC = 97;
    public final static int VIEWERRPCTYPE_CLEARCACHERPC = 98;
    public final static int VIEWERRPCTYPE_CLEARCACHEFORALLENGINESRPC = 99;
    public final static int VIEWERRPCTYPE_SETVIEWEXTENTSTYPERPC = 100;
    public final static int VIEWERRPCTYPE_CLEARREFLINESRPC = 101;
    public final static int VIEWERRPCTYPE_SETRENDERINGATTRIBUTESRPC = 102;
    public final static int VIEWERRPCTYPE_DATABASEQUERYRPC = 103;
    public final static int VIEWERRPCTYPE_POINTQUERYRPC = 104;
    public final static int VIEWERRPCTYPE_LINEQUERYRPC = 105;
    public final static int VIEWERRPCTYPE_CLONEWINDOWRPC = 106;
    public final static int VIEWERRPCTYPE_SETMATERIALATTRIBUTESRPC = 107;
    public final static int VIEWERRPCTYPE_SETDEFAULTMATERIALATTRIBUTESRPC = 108;
    public final static int VIEWERRPCTYPE_RESETMATERIALATTRIBUTESRPC = 109;
    public final static int VIEWERRPCTYPE_SETPLOTDATABASESTATERPC = 110;
    public final static int VIEWERRPCTYPE_DELETEPLOTDATABASEKEYFRAMERPC = 111;
    public final static int VIEWERRPCTYPE_MOVEPLOTDATABASEKEYFRAMERPC = 112;
    public final static int VIEWERRPCTYPE_CLEARVIEWKEYFRAMESRPC = 113;
    public final static int VIEWERRPCTYPE_DELETEVIEWKEYFRAMERPC = 114;
    public final static int VIEWERRPCTYPE_MOVEVIEWKEYFRAMERPC = 115;
    public final static int VIEWERRPCTYPE_SETVIEWKEYFRAMERPC = 116;
    public final static int VIEWERRPCTYPE_OPENMDSERVERRPC = 117;
    public final static int VIEWERRPCTYPE_ENABLETOOLBARRPC = 118;
    public final static int VIEWERRPCTYPE_HIDETOOLBARSRPC = 119;
    public final static int VIEWERRPCTYPE_HIDETOOLBARSFORALLWINDOWSRPC = 120;
    public final static int VIEWERRPCTYPE_SHOWTOOLBARSRPC = 121;
    public final static int VIEWERRPCTYPE_SHOWTOOLBARSFORALLWINDOWSRPC = 122;
    public final static int VIEWERRPCTYPE_SETTOOLBARICONSIZERPC = 123;
    public final static int VIEWERRPCTYPE_SAVEVIEWRPC = 124;
    public final static int VIEWERRPCTYPE_SETGLOBALLINEOUTATTRIBUTESRPC = 125;
    public final static int VIEWERRPCTYPE_SETPICKATTRIBUTESRPC = 126;
    public final static int VIEWERRPCTYPE_EXPORTCOLORTABLERPC = 127;
    public final static int VIEWERRPCTYPE_EXPORTENTIRESTATERPC = 128;
    public final static int VIEWERRPCTYPE_IMPORTENTIRESTATERPC = 129;
    public final static int VIEWERRPCTYPE_IMPORTENTIRESTATEWITHDIFFERENTSOURCESRPC = 130;
    public final static int VIEWERRPCTYPE_RESETPICKATTRIBUTESRPC = 131;
    public final static int VIEWERRPCTYPE_ADDANNOTATIONOBJECTRPC = 132;
    public final static int VIEWERRPCTYPE_HIDEACTIVEANNOTATIONOBJECTSRPC = 133;
    public final static int VIEWERRPCTYPE_DELETEACTIVEANNOTATIONOBJECTSRPC = 134;
    public final static int VIEWERRPCTYPE_RAISEACTIVEANNOTATIONOBJECTSRPC = 135;
    public final static int VIEWERRPCTYPE_LOWERACTIVEANNOTATIONOBJECTSRPC = 136;
    public final static int VIEWERRPCTYPE_SETANNOTATIONOBJECTOPTIONSRPC = 137;
    public final static int VIEWERRPCTYPE_SETDEFAULTANNOTATIONOBJECTLISTRPC = 138;
    public final static int VIEWERRPCTYPE_RESETANNOTATIONOBJECTLISTRPC = 139;
    public final static int VIEWERRPCTYPE_RESETPICKLETTERRPC = 140;
    public final static int VIEWERRPCTYPE_SETDEFAULTPICKATTRIBUTESRPC = 141;
    public final static int VIEWERRPCTYPE_CHOOSECENTEROFROTATIONRPC = 142;
    public final static int VIEWERRPCTYPE_SETCENTEROFROTATIONRPC = 143;
    public final static int VIEWERRPCTYPE_SETQUERYOVERTIMEATTRIBUTESRPC = 144;
    public final static int VIEWERRPCTYPE_SETDEFAULTQUERYOVERTIMEATTRIBUTESRPC = 145;
    public final static int VIEWERRPCTYPE_RESETQUERYOVERTIMEATTRIBUTESRPC = 146;
    public final static int VIEWERRPCTYPE_RESETLINEOUTCOLORRPC = 147;
    public final static int VIEWERRPCTYPE_SETINTERACTORATTRIBUTESRPC = 148;
    public final static int VIEWERRPCTYPE_SETDEFAULTINTERACTORATTRIBUTESRPC = 149;
    public final static int VIEWERRPCTYPE_RESETINTERACTORATTRIBUTESRPC = 150;
    public final static int VIEWERRPCTYPE_GETPROCINFORPC = 151;
    public final static int VIEWERRPCTYPE_SENDSIMULATIONCOMMANDRPC = 152;
    public final static int VIEWERRPCTYPE_UPDATEDBPLUGININFORPC = 153;
    public final static int VIEWERRPCTYPE_EXPORTDBRPC = 154;
    public final static int VIEWERRPCTYPE_SETTRYHARDERCYCLESTIMESRPC = 155;
    public final static int VIEWERRPCTYPE_OPENCLIENTRPC = 156;
    public final static int VIEWERRPCTYPE_SUPPRESSQUERYOUTPUTRPC = 157;
    public final static int VIEWERRPCTYPE_SETQUERYFLOATFORMATRPC = 158;
    public final static int VIEWERRPCTYPE_SETMESHMANAGEMENTATTRIBUTESRPC = 159;
    public final static int VIEWERRPCTYPE_SETDEFAULTMESHMANAGEMENTATTRIBUTESRPC = 160;
    public final static int VIEWERRPCTYPE_RESETMESHMANAGEMENTATTRIBUTESRPC = 161;
    public final static int VIEWERRPCTYPE_RESIZEWINDOWRPC = 162;
    public final static int VIEWERRPCTYPE_MOVEWINDOWRPC = 163;
    public final static int VIEWERRPCTYPE_MOVEANDRESIZEWINDOWRPC = 164;
    public final static int VIEWERRPCTYPE_SETSTATELOGGINGRPC = 165;
    public final static int VIEWERRPCTYPE_CONSTRUCTDDFRPC = 166;
    public final static int VIEWERRPCTYPE_UPDATEPLOTINFOATTSRPC = 167;
    public final static int VIEWERRPCTYPE_REQUESTMETADATARPC = 168;
    public final static int VIEWERRPCTYPE_SETTREATALLDBSASTIMEVARYINGRPC = 169;
    public final static int VIEWERRPCTYPE_SETCREATEMESHQUALITYEXPRESSIONSRPC = 170;
    public final static int VIEWERRPCTYPE_SETCREATETIMEDERIVATIVEEXPRESSIONSRPC = 171;
    public final static int VIEWERRPCTYPE_SETCREATEVECTORMAGNITUDEEXPRESSIONSRPC = 172;
    public final static int VIEWERRPCTYPE_COPYACTIVEPLOTSRPC = 173;
    public final static int VIEWERRPCTYPE_MAXRPC = 174;
    public final static int VIEWERRPCTYPE_DISCONNECTPLOTFROMTIMESLIDERRPC = 175;


    public ViewerRPC()
    {
        super(33);

        RPCType = VIEWERRPCTYPE_CLOSERPC;
        windowLayout = 1;
        windowId = 0;
        windowMode = 0;
        windowArea = new String("");
        database = new String("");
        programHost = new String("");
        programSim = new String("");
        programOptions = new Vector();
        nFrames = 0;
        stateNumber = 0;
        frameRange = new int[2];
        frameRange[0] = 0;
        frameRange[1] = 0;
        frame = 0;
        plotType = 0;
        operatorType = 0;
        variable = new String("");
        activePlotIds = new Vector();
        activeOperatorIds = new Vector();
        expandedPlotIds = new Vector();
        colorTableName = new String("");
        queryName = new String("");
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
        stringArg1 = new String("");
        stringArg2 = new String("");
        doubleArg1 = new Vector();
        doubleArg1.addElement(new Double(0));
        doubleArg2 = new Vector();
        doubleArg2.addElement(new Double(0));
    }

    public ViewerRPC(ViewerRPC obj)
    {
        super(33);

        int i;

        RPCType = obj.RPCType;
        windowLayout = obj.windowLayout;
        windowId = obj.windowId;
        windowMode = obj.windowMode;
        windowArea = new String(obj.windowArea);
        database = new String(obj.database);
        programHost = new String(obj.programHost);
        programSim = new String(obj.programSim);
        programOptions = new Vector(obj.programOptions.size());
        for(i = 0; i < obj.programOptions.size(); ++i)
            programOptions.addElement(new String((String)obj.programOptions.elementAt(i)));

        nFrames = obj.nFrames;
        stateNumber = obj.stateNumber;
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
        stringArg1 = new String(obj.stringArg1);
        stringArg2 = new String(obj.stringArg2);
        doubleArg1 = new Vector(obj.doubleArg1.size());
        for(i = 0; i < obj.doubleArg1.size(); ++i)
        {
            Double dv = (Double)obj.doubleArg1.elementAt(i);
            doubleArg1.addElement(new Double(dv.doubleValue()));
        }

        doubleArg2 = new Vector(obj.doubleArg2.size());
        for(i = 0; i < obj.doubleArg2.size(); ++i)
        {
            Double dv = (Double)obj.doubleArg2.elementAt(i);
            doubleArg2.addElement(new Double(dv.doubleValue()));
        }


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
                (programSim == obj.programSim) &&
                (programOptions == obj.programOptions) &&
                (nFrames == obj.nFrames) &&
                (stateNumber == obj.stateNumber) &&
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
                (intArg3 == obj.intArg3) &&
                (stringArg1 == obj.stringArg1) &&
                (stringArg2 == obj.stringArg2) &&
                (doubleArg1 == obj.doubleArg1) &&
                (doubleArg2 == obj.doubleArg2));
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

    public void SetProgramSim(String programSim_)
    {
        programSim = programSim_;
        Select(7);
    }

    public void SetProgramOptions(Vector programOptions_)
    {
        programOptions = programOptions_;
        Select(8);
    }

    public void SetNFrames(int nFrames_)
    {
        nFrames = nFrames_;
        Select(9);
    }

    public void SetStateNumber(int stateNumber_)
    {
        stateNumber = stateNumber_;
        Select(10);
    }

    public void SetFrameRange(int[] frameRange_)
    {
        frameRange[0] = frameRange_[0];
        frameRange[1] = frameRange_[1];
        Select(11);
    }

    public void SetFrameRange(int e0, int e1)
    {
        frameRange[0] = e0;
        frameRange[1] = e1;
        Select(11);
    }

    public void SetFrame(int frame_)
    {
        frame = frame_;
        Select(12);
    }

    public void SetPlotType(int plotType_)
    {
        plotType = plotType_;
        Select(13);
    }

    public void SetOperatorType(int operatorType_)
    {
        operatorType = operatorType_;
        Select(14);
    }

    public void SetVariable(String variable_)
    {
        variable = variable_;
        Select(15);
    }

    public void SetActivePlotIds(Vector activePlotIds_)
    {
        activePlotIds = activePlotIds_;
        Select(16);
    }

    public void SetActiveOperatorIds(Vector activeOperatorIds_)
    {
        activeOperatorIds = activeOperatorIds_;
        Select(17);
    }

    public void SetExpandedPlotIds(Vector expandedPlotIds_)
    {
        expandedPlotIds = expandedPlotIds_;
        Select(18);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(19);
    }

    public void SetQueryName(String queryName_)
    {
        queryName = queryName_;
        Select(20);
    }

    public void SetQueryPoint1(double[] queryPoint1_)
    {
        queryPoint1[0] = queryPoint1_[0];
        queryPoint1[1] = queryPoint1_[1];
        queryPoint1[2] = queryPoint1_[2];
        Select(21);
    }

    public void SetQueryPoint1(double e0, double e1, double e2)
    {
        queryPoint1[0] = e0;
        queryPoint1[1] = e1;
        queryPoint1[2] = e2;
        Select(21);
    }

    public void SetQueryPoint2(double[] queryPoint2_)
    {
        queryPoint2[0] = queryPoint2_[0];
        queryPoint2[1] = queryPoint2_[1];
        queryPoint2[2] = queryPoint2_[2];
        Select(22);
    }

    public void SetQueryPoint2(double e0, double e1, double e2)
    {
        queryPoint2[0] = e0;
        queryPoint2[1] = e1;
        queryPoint2[2] = e2;
        Select(22);
    }

    public void SetQueryVariables(Vector queryVariables_)
    {
        queryVariables = queryVariables_;
        Select(23);
    }

    public void SetToolId(int toolId_)
    {
        toolId = toolId_;
        Select(24);
    }

    public void SetBoolFlag(boolean boolFlag_)
    {
        boolFlag = boolFlag_;
        Select(25);
    }

    public void SetIntArg1(int intArg1_)
    {
        intArg1 = intArg1_;
        Select(26);
    }

    public void SetIntArg2(int intArg2_)
    {
        intArg2 = intArg2_;
        Select(27);
    }

    public void SetIntArg3(int intArg3_)
    {
        intArg3 = intArg3_;
        Select(28);
    }

    public void SetStringArg1(String stringArg1_)
    {
        stringArg1 = stringArg1_;
        Select(29);
    }

    public void SetStringArg2(String stringArg2_)
    {
        stringArg2 = stringArg2_;
        Select(30);
    }

    public void SetDoubleArg1(Vector doubleArg1_)
    {
        doubleArg1 = doubleArg1_;
        Select(31);
    }

    public void SetDoubleArg2(Vector doubleArg2_)
    {
        doubleArg2 = doubleArg2_;
        Select(32);
    }

    // Property getting methods
    public int      GetRPCType() { return RPCType; }
    public int      GetWindowLayout() { return windowLayout; }
    public int      GetWindowId() { return windowId; }
    public int      GetWindowMode() { return windowMode; }
    public String   GetWindowArea() { return windowArea; }
    public String   GetDatabase() { return database; }
    public String   GetProgramHost() { return programHost; }
    public String   GetProgramSim() { return programSim; }
    public Vector   GetProgramOptions() { return programOptions; }
    public int      GetNFrames() { return nFrames; }
    public int      GetStateNumber() { return stateNumber; }
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
    public String   GetStringArg1() { return stringArg1; }
    public String   GetStringArg2() { return stringArg2; }
    public Vector   GetDoubleArg1() { return doubleArg1; }
    public Vector   GetDoubleArg2() { return doubleArg2; }

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
            buf.WriteString(programSim);
        if(WriteSelect(8, buf))
            buf.WriteStringVector(programOptions);
        if(WriteSelect(9, buf))
            buf.WriteInt(nFrames);
        if(WriteSelect(10, buf))
            buf.WriteInt(stateNumber);
        if(WriteSelect(11, buf))
            buf.WriteIntArray(frameRange);
        if(WriteSelect(12, buf))
            buf.WriteInt(frame);
        if(WriteSelect(13, buf))
            buf.WriteInt(plotType);
        if(WriteSelect(14, buf))
            buf.WriteInt(operatorType);
        if(WriteSelect(15, buf))
            buf.WriteString(variable);
        if(WriteSelect(16, buf))
            buf.WriteIntVector(activePlotIds);
        if(WriteSelect(17, buf))
            buf.WriteIntVector(activeOperatorIds);
        if(WriteSelect(18, buf))
            buf.WriteIntVector(expandedPlotIds);
        if(WriteSelect(19, buf))
            buf.WriteString(colorTableName);
        if(WriteSelect(20, buf))
            buf.WriteString(queryName);
        if(WriteSelect(21, buf))
            buf.WriteDoubleArray(queryPoint1);
        if(WriteSelect(22, buf))
            buf.WriteDoubleArray(queryPoint2);
        if(WriteSelect(23, buf))
            buf.WriteStringVector(queryVariables);
        if(WriteSelect(24, buf))
            buf.WriteInt(toolId);
        if(WriteSelect(25, buf))
            buf.WriteBool(boolFlag);
        if(WriteSelect(26, buf))
            buf.WriteInt(intArg1);
        if(WriteSelect(27, buf))
            buf.WriteInt(intArg2);
        if(WriteSelect(28, buf))
            buf.WriteInt(intArg3);
        if(WriteSelect(29, buf))
            buf.WriteString(stringArg1);
        if(WriteSelect(30, buf))
            buf.WriteString(stringArg2);
        if(WriteSelect(31, buf))
            buf.WriteDoubleVector(doubleArg1);
        if(WriteSelect(32, buf))
            buf.WriteDoubleVector(doubleArg2);
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
                SetProgramSim(buf.ReadString());
                break;
            case 8:
                SetProgramOptions(buf.ReadStringVector());
                break;
            case 9:
                SetNFrames(buf.ReadInt());
                break;
            case 10:
                SetStateNumber(buf.ReadInt());
                break;
            case 11:
                SetFrameRange(buf.ReadIntArray());
                break;
            case 12:
                SetFrame(buf.ReadInt());
                break;
            case 13:
                SetPlotType(buf.ReadInt());
                break;
            case 14:
                SetOperatorType(buf.ReadInt());
                break;
            case 15:
                SetVariable(buf.ReadString());
                break;
            case 16:
                SetActivePlotIds(buf.ReadIntVector());
                break;
            case 17:
                SetActiveOperatorIds(buf.ReadIntVector());
                break;
            case 18:
                SetExpandedPlotIds(buf.ReadIntVector());
                break;
            case 19:
                SetColorTableName(buf.ReadString());
                break;
            case 20:
                SetQueryName(buf.ReadString());
                break;
            case 21:
                SetQueryPoint1(buf.ReadDoubleArray());
                break;
            case 22:
                SetQueryPoint2(buf.ReadDoubleArray());
                break;
            case 23:
                SetQueryVariables(buf.ReadStringVector());
                break;
            case 24:
                SetToolId(buf.ReadInt());
                break;
            case 25:
                SetBoolFlag(buf.ReadBool());
                break;
            case 26:
                SetIntArg1(buf.ReadInt());
                break;
            case 27:
                SetIntArg2(buf.ReadInt());
                break;
            case 28:
                SetIntArg3(buf.ReadInt());
                break;
            case 29:
                SetStringArg1(buf.ReadString());
                break;
            case 30:
                SetStringArg2(buf.ReadString());
                break;
            case 31:
                SetDoubleArg1(buf.ReadDoubleVector());
                break;
            case 32:
                SetDoubleArg2(buf.ReadDoubleVector());
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
    private String   programSim;
    private Vector   programOptions; // vector of String objects
    private int      nFrames;
    private int      stateNumber;
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
    private String   stringArg1;
    private String   stringArg2;
    private Vector   doubleArg1; // vector of Double objects
    private Vector   doubleArg2; // vector of Double objects
}

