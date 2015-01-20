// ***************************************************************************
//
// Copyright (c) 2000 - 2015, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-442911
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

import java.lang.Integer;
import java.util.Vector;

// ****************************************************************************
// Class: Plot
//
// Purpose:
//    This class is a plot element in a plot list.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class Plot extends AttributeSubject
{
    private static int Plot_numAdditionalAtts = 21;

    // Enum values
    public final static int STATETYPE_NEWLYCREATED = 0;
    public final static int STATETYPE_PENDING = 1;
    public final static int STATETYPE_COMPLETED = 2;
    public final static int STATETYPE_ERROR = 3;


    public Plot()
    {
        super(Plot_numAdditionalAtts);

        stateType = STATETYPE_NEWLYCREATED;
        plotType = 0;
        plotName = new String("");
        activeFlag = false;
        hiddenFlag = false;
        expandedFlag = false;
        plotVar = new String("notset");
        databaseName = new String("notset");
        operators = new Vector();
        operatorNames = new Vector();
        activeOperator = -1;
        id = -1;
        embeddedPlotId = -1;
        beginFrame = -999;
        endFrame = 999;
        keyframes = new Vector();
        databaseKeyframes = new Vector();
        isFromSimulation = false;
        followsTime = true;
        description = new String("");
        selection = new String("");
    }

    public Plot(int nMoreFields)
    {
        super(Plot_numAdditionalAtts + nMoreFields);

        stateType = STATETYPE_NEWLYCREATED;
        plotType = 0;
        plotName = new String("");
        activeFlag = false;
        hiddenFlag = false;
        expandedFlag = false;
        plotVar = new String("notset");
        databaseName = new String("notset");
        operators = new Vector();
        operatorNames = new Vector();
        activeOperator = -1;
        id = -1;
        embeddedPlotId = -1;
        beginFrame = -999;
        endFrame = 999;
        keyframes = new Vector();
        databaseKeyframes = new Vector();
        isFromSimulation = false;
        followsTime = true;
        description = new String("");
        selection = new String("");
    }

    public Plot(Plot obj)
    {
        super(Plot_numAdditionalAtts);

        int i;

        stateType = obj.stateType;
        plotType = obj.plotType;
        plotName = new String(obj.plotName);
        activeFlag = obj.activeFlag;
        hiddenFlag = obj.hiddenFlag;
        expandedFlag = obj.expandedFlag;
        plotVar = new String(obj.plotVar);
        databaseName = new String(obj.databaseName);
        operators = new Vector();
        for(i = 0; i < obj.operators.size(); ++i)
        {
            Integer iv = (Integer)obj.operators.elementAt(i);
            operators.addElement(new Integer(iv.intValue()));
        }
        operatorNames = new Vector(obj.operatorNames.size());
        for(i = 0; i < obj.operatorNames.size(); ++i)
            operatorNames.addElement(new String((String)obj.operatorNames.elementAt(i)));

        activeOperator = obj.activeOperator;
        id = obj.id;
        embeddedPlotId = obj.embeddedPlotId;
        beginFrame = obj.beginFrame;
        endFrame = obj.endFrame;
        keyframes = new Vector();
        for(i = 0; i < obj.keyframes.size(); ++i)
        {
            Integer iv = (Integer)obj.keyframes.elementAt(i);
            keyframes.addElement(new Integer(iv.intValue()));
        }
        databaseKeyframes = new Vector();
        for(i = 0; i < obj.databaseKeyframes.size(); ++i)
        {
            Integer iv = (Integer)obj.databaseKeyframes.elementAt(i);
            databaseKeyframes.addElement(new Integer(iv.intValue()));
        }
        isFromSimulation = obj.isFromSimulation;
        followsTime = obj.followsTime;
        description = new String(obj.description);
        selection = new String(obj.selection);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return Plot_numAdditionalAtts;
    }

    public boolean equals(Plot obj)
    {
        int i;

        // Compare the elements in the operators vector.
        boolean operators_equal = (obj.operators.size() == operators.size());
        for(i = 0; (i < operators.size()) && operators_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer operators1 = (Integer)operators.elementAt(i);
            Integer operators2 = (Integer)obj.operators.elementAt(i);
            operators_equal = operators1.equals(operators2);
        }
        // Compare the elements in the operatorNames vector.
        boolean operatorNames_equal = (obj.operatorNames.size() == operatorNames.size());
        for(i = 0; (i < operatorNames.size()) && operatorNames_equal; ++i)
        {
            // Make references to String from Object.
            String operatorNames1 = (String)operatorNames.elementAt(i);
            String operatorNames2 = (String)obj.operatorNames.elementAt(i);
            operatorNames_equal = operatorNames1.equals(operatorNames2);
        }
        // Compare the elements in the keyframes vector.
        boolean keyframes_equal = (obj.keyframes.size() == keyframes.size());
        for(i = 0; (i < keyframes.size()) && keyframes_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer keyframes1 = (Integer)keyframes.elementAt(i);
            Integer keyframes2 = (Integer)obj.keyframes.elementAt(i);
            keyframes_equal = keyframes1.equals(keyframes2);
        }
        // Compare the elements in the databaseKeyframes vector.
        boolean databaseKeyframes_equal = (obj.databaseKeyframes.size() == databaseKeyframes.size());
        for(i = 0; (i < databaseKeyframes.size()) && databaseKeyframes_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer databaseKeyframes1 = (Integer)databaseKeyframes.elementAt(i);
            Integer databaseKeyframes2 = (Integer)obj.databaseKeyframes.elementAt(i);
            databaseKeyframes_equal = databaseKeyframes1.equals(databaseKeyframes2);
        }
        // Create the return value
        return ((stateType == obj.stateType) &&
                (plotType == obj.plotType) &&
                (plotName.equals(obj.plotName)) &&
                (activeFlag == obj.activeFlag) &&
                (hiddenFlag == obj.hiddenFlag) &&
                (expandedFlag == obj.expandedFlag) &&
                (plotVar.equals(obj.plotVar)) &&
                (databaseName.equals(obj.databaseName)) &&
                operators_equal &&
                operatorNames_equal &&
                (activeOperator == obj.activeOperator) &&
                (id == obj.id) &&
                (embeddedPlotId == obj.embeddedPlotId) &&
                (beginFrame == obj.beginFrame) &&
                (endFrame == obj.endFrame) &&
                keyframes_equal &&
                databaseKeyframes_equal &&
                (isFromSimulation == obj.isFromSimulation) &&
                (followsTime == obj.followsTime) &&
                (description.equals(obj.description)) &&
                (selection.equals(obj.selection)));
    }

    // Property setting methods
    public void SetStateType(int stateType_)
    {
        stateType = stateType_;
        Select(0);
    }

    public void SetPlotType(int plotType_)
    {
        plotType = plotType_;
        Select(1);
    }

    public void SetPlotName(String plotName_)
    {
        plotName = plotName_;
        Select(2);
    }

    public void SetActiveFlag(boolean activeFlag_)
    {
        activeFlag = activeFlag_;
        Select(3);
    }

    public void SetHiddenFlag(boolean hiddenFlag_)
    {
        hiddenFlag = hiddenFlag_;
        Select(4);
    }

    public void SetExpandedFlag(boolean expandedFlag_)
    {
        expandedFlag = expandedFlag_;
        Select(5);
    }

    public void SetPlotVar(String plotVar_)
    {
        plotVar = plotVar_;
        Select(6);
    }

    public void SetDatabaseName(String databaseName_)
    {
        databaseName = databaseName_;
        Select(7);
    }

    public void SetOperators(Vector operators_)
    {
        operators = operators_;
        Select(8);
    }

    public void SetOperatorNames(Vector operatorNames_)
    {
        operatorNames = operatorNames_;
        Select(9);
    }

    public void SetActiveOperator(int activeOperator_)
    {
        activeOperator = activeOperator_;
        Select(10);
    }

    public void SetId(int id_)
    {
        id = id_;
        Select(11);
    }

    public void SetEmbeddedPlotId(int embeddedPlotId_)
    {
        embeddedPlotId = embeddedPlotId_;
        Select(12);
    }

    public void SetBeginFrame(int beginFrame_)
    {
        beginFrame = beginFrame_;
        Select(13);
    }

    public void SetEndFrame(int endFrame_)
    {
        endFrame = endFrame_;
        Select(14);
    }

    public void SetKeyframes(Vector keyframes_)
    {
        keyframes = keyframes_;
        Select(15);
    }

    public void SetDatabaseKeyframes(Vector databaseKeyframes_)
    {
        databaseKeyframes = databaseKeyframes_;
        Select(16);
    }

    public void SetIsFromSimulation(boolean isFromSimulation_)
    {
        isFromSimulation = isFromSimulation_;
        Select(17);
    }

    public void SetFollowsTime(boolean followsTime_)
    {
        followsTime = followsTime_;
        Select(18);
    }

    public void SetDescription(String description_)
    {
        description = description_;
        Select(19);
    }

    public void SetSelection(String selection_)
    {
        selection = selection_;
        Select(20);
    }

    // Property getting methods
    public int     GetStateType() { return stateType; }
    public int     GetPlotType() { return plotType; }
    public String  GetPlotName() { return plotName; }
    public boolean GetActiveFlag() { return activeFlag; }
    public boolean GetHiddenFlag() { return hiddenFlag; }
    public boolean GetExpandedFlag() { return expandedFlag; }
    public String  GetPlotVar() { return plotVar; }
    public String  GetDatabaseName() { return databaseName; }
    public Vector  GetOperators() { return operators; }
    public Vector  GetOperatorNames() { return operatorNames; }
    public int     GetActiveOperator() { return activeOperator; }
    public int     GetId() { return id; }
    public int     GetEmbeddedPlotId() { return embeddedPlotId; }
    public int     GetBeginFrame() { return beginFrame; }
    public int     GetEndFrame() { return endFrame; }
    public Vector  GetKeyframes() { return keyframes; }
    public Vector  GetDatabaseKeyframes() { return databaseKeyframes; }
    public boolean GetIsFromSimulation() { return isFromSimulation; }
    public boolean GetFollowsTime() { return followsTime; }
    public String  GetDescription() { return description; }
    public String  GetSelection() { return selection; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(stateType);
        if(WriteSelect(1, buf))
            buf.WriteInt(plotType);
        if(WriteSelect(2, buf))
            buf.WriteString(plotName);
        if(WriteSelect(3, buf))
            buf.WriteBool(activeFlag);
        if(WriteSelect(4, buf))
            buf.WriteBool(hiddenFlag);
        if(WriteSelect(5, buf))
            buf.WriteBool(expandedFlag);
        if(WriteSelect(6, buf))
            buf.WriteString(plotVar);
        if(WriteSelect(7, buf))
            buf.WriteString(databaseName);
        if(WriteSelect(8, buf))
            buf.WriteIntVector(operators);
        if(WriteSelect(9, buf))
            buf.WriteStringVector(operatorNames);
        if(WriteSelect(10, buf))
            buf.WriteInt(activeOperator);
        if(WriteSelect(11, buf))
            buf.WriteInt(id);
        if(WriteSelect(12, buf))
            buf.WriteInt(embeddedPlotId);
        if(WriteSelect(13, buf))
            buf.WriteInt(beginFrame);
        if(WriteSelect(14, buf))
            buf.WriteInt(endFrame);
        if(WriteSelect(15, buf))
            buf.WriteIntVector(keyframes);
        if(WriteSelect(16, buf))
            buf.WriteIntVector(databaseKeyframes);
        if(WriteSelect(17, buf))
            buf.WriteBool(isFromSimulation);
        if(WriteSelect(18, buf))
            buf.WriteBool(followsTime);
        if(WriteSelect(19, buf))
            buf.WriteString(description);
        if(WriteSelect(20, buf))
            buf.WriteString(selection);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetStateType(buf.ReadInt());
            break;
        case 1:
            SetPlotType(buf.ReadInt());
            break;
        case 2:
            SetPlotName(buf.ReadString());
            break;
        case 3:
            SetActiveFlag(buf.ReadBool());
            break;
        case 4:
            SetHiddenFlag(buf.ReadBool());
            break;
        case 5:
            SetExpandedFlag(buf.ReadBool());
            break;
        case 6:
            SetPlotVar(buf.ReadString());
            break;
        case 7:
            SetDatabaseName(buf.ReadString());
            break;
        case 8:
            SetOperators(buf.ReadIntVector());
            break;
        case 9:
            SetOperatorNames(buf.ReadStringVector());
            break;
        case 10:
            SetActiveOperator(buf.ReadInt());
            break;
        case 11:
            SetId(buf.ReadInt());
            break;
        case 12:
            SetEmbeddedPlotId(buf.ReadInt());
            break;
        case 13:
            SetBeginFrame(buf.ReadInt());
            break;
        case 14:
            SetEndFrame(buf.ReadInt());
            break;
        case 15:
            SetKeyframes(buf.ReadIntVector());
            break;
        case 16:
            SetDatabaseKeyframes(buf.ReadIntVector());
            break;
        case 17:
            SetIsFromSimulation(buf.ReadBool());
            break;
        case 18:
            SetFollowsTime(buf.ReadBool());
            break;
        case 19:
            SetDescription(buf.ReadString());
            break;
        case 20:
            SetSelection(buf.ReadString());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "stateType = ";
        if(stateType == STATETYPE_NEWLYCREATED)
            str = str + "STATETYPE_NEWLYCREATED";
        if(stateType == STATETYPE_PENDING)
            str = str + "STATETYPE_PENDING";
        if(stateType == STATETYPE_COMPLETED)
            str = str + "STATETYPE_COMPLETED";
        if(stateType == STATETYPE_ERROR)
            str = str + "STATETYPE_ERROR";
        str = str + "\n";
        str = str + intToString("plotType", plotType, indent) + "\n";
        str = str + stringToString("plotName", plotName, indent) + "\n";
        str = str + boolToString("activeFlag", activeFlag, indent) + "\n";
        str = str + boolToString("hiddenFlag", hiddenFlag, indent) + "\n";
        str = str + boolToString("expandedFlag", expandedFlag, indent) + "\n";
        str = str + stringToString("plotVar", plotVar, indent) + "\n";
        str = str + stringToString("databaseName", databaseName, indent) + "\n";
        str = str + intVectorToString("operators", operators, indent) + "\n";
        str = str + stringVectorToString("operatorNames", operatorNames, indent) + "\n";
        str = str + intToString("activeOperator", activeOperator, indent) + "\n";
        str = str + intToString("id", id, indent) + "\n";
        str = str + intToString("embeddedPlotId", embeddedPlotId, indent) + "\n";
        str = str + intToString("beginFrame", beginFrame, indent) + "\n";
        str = str + intToString("endFrame", endFrame, indent) + "\n";
        str = str + intVectorToString("keyframes", keyframes, indent) + "\n";
        str = str + intVectorToString("databaseKeyframes", databaseKeyframes, indent) + "\n";
        str = str + boolToString("isFromSimulation", isFromSimulation, indent) + "\n";
        str = str + boolToString("followsTime", followsTime, indent) + "\n";
        str = str + stringToString("description", description, indent) + "\n";
        str = str + stringToString("selection", selection, indent) + "\n";
        return str;
    }


    // Attributes
    private int     stateType;
    private int     plotType;
    private String  plotName;
    private boolean activeFlag;
    private boolean hiddenFlag;
    private boolean expandedFlag;
    private String  plotVar;
    private String  databaseName;
    private Vector  operators; // vector of Integer objects
    private Vector  operatorNames; // vector of String objects
    private int     activeOperator;
    private int     id;
    private int     embeddedPlotId;
    private int     beginFrame;
    private int     endFrame;
    private Vector  keyframes; // vector of Integer objects
    private Vector  databaseKeyframes; // vector of Integer objects
    private boolean isFromSimulation;
    private boolean followsTime;
    private String  description;
    private String  selection;
}

