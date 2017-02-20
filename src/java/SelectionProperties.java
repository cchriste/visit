// ***************************************************************************
//
// Copyright (c) 2000 - 2017, Lawrence Livermore National Security, LLC
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

import java.util.Vector;
import java.lang.Double;

// ****************************************************************************
// Class: SelectionProperties
//
// Purpose:
//    Contains attributes for a selection
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class SelectionProperties extends AttributeSubject
{
    private static int SelectionProperties_numAdditionalAtts = 19;

    // Enum values
    public final static int SELECTIONTYPE_BASICSELECTION = 0;
    public final static int SELECTIONTYPE_CUMULATIVEQUERYSELECTION = 1;

    public final static int COMBINATIONTYPE_COMBINEAND = 0;
    public final static int COMBINATIONTYPE_COMBINEOR = 1;

    public final static int HISTOGRAMTYPE_HISTOGRAMTIME = 0;
    public final static int HISTOGRAMTYPE_HISTOGRAMMATCHES = 1;
    public final static int HISTOGRAMTYPE_HISTOGRAMID = 2;
    public final static int HISTOGRAMTYPE_HISTOGRAMVARIABLE = 3;

    public final static int IDVARIABLETYPE_USEZONEIDFORID = 0;
    public final static int IDVARIABLETYPE_USEGLOBALZONEIDFORID = 1;
    public final static int IDVARIABLETYPE_USELOCATIONSFORID = 2;
    public final static int IDVARIABLETYPE_USEVARIABLEFORID = 3;


    public SelectionProperties()
    {
        super(SelectionProperties_numAdditionalAtts);

        name = new String("");
        source = new String("");
        host = new String("");
        selectionType = SELECTIONTYPE_BASICSELECTION;
        idVariableType = IDVARIABLETYPE_USEZONEIDFORID;
        idVariable = new String("");
        variables = new Vector();
        variableMins = new Vector();
        variableMaxs = new Vector();
        minTimeState = 0;
        maxTimeState = -1;
        timeStateStride = 1;
        combineRule = COMBINATIONTYPE_COMBINEOR;
        histogramType = HISTOGRAMTYPE_HISTOGRAMTIME;
        histogramNumBins = 10;
        histogramAutoScaleNumBins = false;
        histogramStartBin = 0;
        histogramEndBin = 9;
        histogramVariable = new String("");
    }

    public SelectionProperties(int nMoreFields)
    {
        super(SelectionProperties_numAdditionalAtts + nMoreFields);

        name = new String("");
        source = new String("");
        host = new String("");
        selectionType = SELECTIONTYPE_BASICSELECTION;
        idVariableType = IDVARIABLETYPE_USEZONEIDFORID;
        idVariable = new String("");
        variables = new Vector();
        variableMins = new Vector();
        variableMaxs = new Vector();
        minTimeState = 0;
        maxTimeState = -1;
        timeStateStride = 1;
        combineRule = COMBINATIONTYPE_COMBINEOR;
        histogramType = HISTOGRAMTYPE_HISTOGRAMTIME;
        histogramNumBins = 10;
        histogramAutoScaleNumBins = false;
        histogramStartBin = 0;
        histogramEndBin = 9;
        histogramVariable = new String("");
    }

    public SelectionProperties(SelectionProperties obj)
    {
        super(SelectionProperties_numAdditionalAtts);

        int i;

        name = new String(obj.name);
        source = new String(obj.source);
        host = new String(obj.host);
        selectionType = obj.selectionType;
        idVariableType = obj.idVariableType;
        idVariable = new String(obj.idVariable);
        variables = new Vector(obj.variables.size());
        for(i = 0; i < obj.variables.size(); ++i)
            variables.addElement(new String((String)obj.variables.elementAt(i)));

        variableMins = new Vector(obj.variableMins.size());
        for(i = 0; i < obj.variableMins.size(); ++i)
        {
            Double dv = (Double)obj.variableMins.elementAt(i);
            variableMins.addElement(new Double(dv.doubleValue()));
        }

        variableMaxs = new Vector(obj.variableMaxs.size());
        for(i = 0; i < obj.variableMaxs.size(); ++i)
        {
            Double dv = (Double)obj.variableMaxs.elementAt(i);
            variableMaxs.addElement(new Double(dv.doubleValue()));
        }

        minTimeState = obj.minTimeState;
        maxTimeState = obj.maxTimeState;
        timeStateStride = obj.timeStateStride;
        combineRule = obj.combineRule;
        histogramType = obj.histogramType;
        histogramNumBins = obj.histogramNumBins;
        histogramAutoScaleNumBins = obj.histogramAutoScaleNumBins;
        histogramStartBin = obj.histogramStartBin;
        histogramEndBin = obj.histogramEndBin;
        histogramVariable = new String(obj.histogramVariable);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return SelectionProperties_numAdditionalAtts;
    }

    public boolean equals(SelectionProperties obj)
    {
        int i;

        // Compare the elements in the variables vector.
        boolean variables_equal = (obj.variables.size() == variables.size());
        for(i = 0; (i < variables.size()) && variables_equal; ++i)
        {
            // Make references to String from Object.
            String variables1 = (String)variables.elementAt(i);
            String variables2 = (String)obj.variables.elementAt(i);
            variables_equal = variables1.equals(variables2);
        }
        // Compare the elements in the variableMins vector.
        boolean variableMins_equal = (obj.variableMins.size() == variableMins.size());
        for(i = 0; (i < variableMins.size()) && variableMins_equal; ++i)
        {
            // Make references to Double from Object.
            Double variableMins1 = (Double)variableMins.elementAt(i);
            Double variableMins2 = (Double)obj.variableMins.elementAt(i);
            variableMins_equal = variableMins1.equals(variableMins2);
        }
        // Compare the elements in the variableMaxs vector.
        boolean variableMaxs_equal = (obj.variableMaxs.size() == variableMaxs.size());
        for(i = 0; (i < variableMaxs.size()) && variableMaxs_equal; ++i)
        {
            // Make references to Double from Object.
            Double variableMaxs1 = (Double)variableMaxs.elementAt(i);
            Double variableMaxs2 = (Double)obj.variableMaxs.elementAt(i);
            variableMaxs_equal = variableMaxs1.equals(variableMaxs2);
        }
        // Create the return value
        return ((name.equals(obj.name)) &&
                (source.equals(obj.source)) &&
                (host.equals(obj.host)) &&
                (selectionType == obj.selectionType) &&
                (idVariableType == obj.idVariableType) &&
                (idVariable.equals(obj.idVariable)) &&
                variables_equal &&
                variableMins_equal &&
                variableMaxs_equal &&
                (minTimeState == obj.minTimeState) &&
                (maxTimeState == obj.maxTimeState) &&
                (timeStateStride == obj.timeStateStride) &&
                (combineRule == obj.combineRule) &&
                (histogramType == obj.histogramType) &&
                (histogramNumBins == obj.histogramNumBins) &&
                (histogramAutoScaleNumBins == obj.histogramAutoScaleNumBins) &&
                (histogramStartBin == obj.histogramStartBin) &&
                (histogramEndBin == obj.histogramEndBin) &&
                (histogramVariable.equals(obj.histogramVariable)));
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetSource(String source_)
    {
        source = source_;
        Select(1);
    }

    public void SetHost(String host_)
    {
        host = host_;
        Select(2);
    }

    public void SetSelectionType(int selectionType_)
    {
        selectionType = selectionType_;
        Select(3);
    }

    public void SetIdVariableType(int idVariableType_)
    {
        idVariableType = idVariableType_;
        Select(4);
    }

    public void SetIdVariable(String idVariable_)
    {
        idVariable = idVariable_;
        Select(5);
    }

    public void SetVariables(Vector variables_)
    {
        variables = variables_;
        Select(6);
    }

    public void SetVariableMins(Vector variableMins_)
    {
        variableMins = variableMins_;
        Select(7);
    }

    public void SetVariableMaxs(Vector variableMaxs_)
    {
        variableMaxs = variableMaxs_;
        Select(8);
    }

    public void SetMinTimeState(int minTimeState_)
    {
        minTimeState = minTimeState_;
        Select(9);
    }

    public void SetMaxTimeState(int maxTimeState_)
    {
        maxTimeState = maxTimeState_;
        Select(10);
    }

    public void SetTimeStateStride(int timeStateStride_)
    {
        timeStateStride = timeStateStride_;
        Select(11);
    }

    public void SetCombineRule(int combineRule_)
    {
        combineRule = combineRule_;
        Select(12);
    }

    public void SetHistogramType(int histogramType_)
    {
        histogramType = histogramType_;
        Select(13);
    }

    public void SetHistogramNumBins(int histogramNumBins_)
    {
        histogramNumBins = histogramNumBins_;
        Select(14);
    }

    public void SetHistogramAutoScaleNumBins(boolean histogramAutoScaleNumBins_)
    {
        histogramAutoScaleNumBins = histogramAutoScaleNumBins_;
        Select(15);
    }

    public void SetHistogramStartBin(int histogramStartBin_)
    {
        histogramStartBin = histogramStartBin_;
        Select(16);
    }

    public void SetHistogramEndBin(int histogramEndBin_)
    {
        histogramEndBin = histogramEndBin_;
        Select(17);
    }

    public void SetHistogramVariable(String histogramVariable_)
    {
        histogramVariable = histogramVariable_;
        Select(18);
    }

    // Property getting methods
    public String  GetName() { return name; }
    public String  GetSource() { return source; }
    public String  GetHost() { return host; }
    public int     GetSelectionType() { return selectionType; }
    public int     GetIdVariableType() { return idVariableType; }
    public String  GetIdVariable() { return idVariable; }
    public Vector  GetVariables() { return variables; }
    public Vector  GetVariableMins() { return variableMins; }
    public Vector  GetVariableMaxs() { return variableMaxs; }
    public int     GetMinTimeState() { return minTimeState; }
    public int     GetMaxTimeState() { return maxTimeState; }
    public int     GetTimeStateStride() { return timeStateStride; }
    public int     GetCombineRule() { return combineRule; }
    public int     GetHistogramType() { return histogramType; }
    public int     GetHistogramNumBins() { return histogramNumBins; }
    public boolean GetHistogramAutoScaleNumBins() { return histogramAutoScaleNumBins; }
    public int     GetHistogramStartBin() { return histogramStartBin; }
    public int     GetHistogramEndBin() { return histogramEndBin; }
    public String  GetHistogramVariable() { return histogramVariable; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
            buf.WriteString(source);
        if(WriteSelect(2, buf))
            buf.WriteString(host);
        if(WriteSelect(3, buf))
            buf.WriteInt(selectionType);
        if(WriteSelect(4, buf))
            buf.WriteInt(idVariableType);
        if(WriteSelect(5, buf))
            buf.WriteString(idVariable);
        if(WriteSelect(6, buf))
            buf.WriteStringVector(variables);
        if(WriteSelect(7, buf))
            buf.WriteDoubleVector(variableMins);
        if(WriteSelect(8, buf))
            buf.WriteDoubleVector(variableMaxs);
        if(WriteSelect(9, buf))
            buf.WriteInt(minTimeState);
        if(WriteSelect(10, buf))
            buf.WriteInt(maxTimeState);
        if(WriteSelect(11, buf))
            buf.WriteInt(timeStateStride);
        if(WriteSelect(12, buf))
            buf.WriteInt(combineRule);
        if(WriteSelect(13, buf))
            buf.WriteInt(histogramType);
        if(WriteSelect(14, buf))
            buf.WriteInt(histogramNumBins);
        if(WriteSelect(15, buf))
            buf.WriteBool(histogramAutoScaleNumBins);
        if(WriteSelect(16, buf))
            buf.WriteInt(histogramStartBin);
        if(WriteSelect(17, buf))
            buf.WriteInt(histogramEndBin);
        if(WriteSelect(18, buf))
            buf.WriteString(histogramVariable);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetName(buf.ReadString());
            break;
        case 1:
            SetSource(buf.ReadString());
            break;
        case 2:
            SetHost(buf.ReadString());
            break;
        case 3:
            SetSelectionType(buf.ReadInt());
            break;
        case 4:
            SetIdVariableType(buf.ReadInt());
            break;
        case 5:
            SetIdVariable(buf.ReadString());
            break;
        case 6:
            SetVariables(buf.ReadStringVector());
            break;
        case 7:
            SetVariableMins(buf.ReadDoubleVector());
            break;
        case 8:
            SetVariableMaxs(buf.ReadDoubleVector());
            break;
        case 9:
            SetMinTimeState(buf.ReadInt());
            break;
        case 10:
            SetMaxTimeState(buf.ReadInt());
            break;
        case 11:
            SetTimeStateStride(buf.ReadInt());
            break;
        case 12:
            SetCombineRule(buf.ReadInt());
            break;
        case 13:
            SetHistogramType(buf.ReadInt());
            break;
        case 14:
            SetHistogramNumBins(buf.ReadInt());
            break;
        case 15:
            SetHistogramAutoScaleNumBins(buf.ReadBool());
            break;
        case 16:
            SetHistogramStartBin(buf.ReadInt());
            break;
        case 17:
            SetHistogramEndBin(buf.ReadInt());
            break;
        case 18:
            SetHistogramVariable(buf.ReadString());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("name", name, indent) + "\n";
        str = str + stringToString("source", source, indent) + "\n";
        str = str + stringToString("host", host, indent) + "\n";
        str = str + indent + "selectionType = ";
        if(selectionType == SELECTIONTYPE_BASICSELECTION)
            str = str + "SELECTIONTYPE_BASICSELECTION";
        if(selectionType == SELECTIONTYPE_CUMULATIVEQUERYSELECTION)
            str = str + "SELECTIONTYPE_CUMULATIVEQUERYSELECTION";
        str = str + "\n";
        str = str + indent + "idVariableType = ";
        if(idVariableType == IDVARIABLETYPE_USEZONEIDFORID)
            str = str + "IDVARIABLETYPE_USEZONEIDFORID";
        if(idVariableType == IDVARIABLETYPE_USEGLOBALZONEIDFORID)
            str = str + "IDVARIABLETYPE_USEGLOBALZONEIDFORID";
        if(idVariableType == IDVARIABLETYPE_USELOCATIONSFORID)
            str = str + "IDVARIABLETYPE_USELOCATIONSFORID";
        if(idVariableType == IDVARIABLETYPE_USEVARIABLEFORID)
            str = str + "IDVARIABLETYPE_USEVARIABLEFORID";
        str = str + "\n";
        str = str + stringToString("idVariable", idVariable, indent) + "\n";
        str = str + stringVectorToString("variables", variables, indent) + "\n";
        str = str + doubleVectorToString("variableMins", variableMins, indent) + "\n";
        str = str + doubleVectorToString("variableMaxs", variableMaxs, indent) + "\n";
        str = str + intToString("minTimeState", minTimeState, indent) + "\n";
        str = str + intToString("maxTimeState", maxTimeState, indent) + "\n";
        str = str + intToString("timeStateStride", timeStateStride, indent) + "\n";
        str = str + indent + "combineRule = ";
        if(combineRule == COMBINATIONTYPE_COMBINEAND)
            str = str + "COMBINATIONTYPE_COMBINEAND";
        if(combineRule == COMBINATIONTYPE_COMBINEOR)
            str = str + "COMBINATIONTYPE_COMBINEOR";
        str = str + "\n";
        str = str + indent + "histogramType = ";
        if(histogramType == HISTOGRAMTYPE_HISTOGRAMTIME)
            str = str + "HISTOGRAMTYPE_HISTOGRAMTIME";
        if(histogramType == HISTOGRAMTYPE_HISTOGRAMMATCHES)
            str = str + "HISTOGRAMTYPE_HISTOGRAMMATCHES";
        if(histogramType == HISTOGRAMTYPE_HISTOGRAMID)
            str = str + "HISTOGRAMTYPE_HISTOGRAMID";
        if(histogramType == HISTOGRAMTYPE_HISTOGRAMVARIABLE)
            str = str + "HISTOGRAMTYPE_HISTOGRAMVARIABLE";
        str = str + "\n";
        str = str + intToString("histogramNumBins", histogramNumBins, indent) + "\n";
        str = str + boolToString("histogramAutoScaleNumBins", histogramAutoScaleNumBins, indent) + "\n";
        str = str + intToString("histogramStartBin", histogramStartBin, indent) + "\n";
        str = str + intToString("histogramEndBin", histogramEndBin, indent) + "\n";
        str = str + stringToString("histogramVariable", histogramVariable, indent) + "\n";
        return str;
    }


    // Attributes
    private String  name;
    private String  source;
    private String  host;
    private int     selectionType;
    private int     idVariableType;
    private String  idVariable;
    private Vector  variables; // vector of String objects
    private Vector  variableMins; // vector of Double objects
    private Vector  variableMaxs; // vector of Double objects
    private int     minTimeState;
    private int     maxTimeState;
    private int     timeStateStride;
    private int     combineRule;
    private int     histogramType;
    private int     histogramNumBins;
    private boolean histogramAutoScaleNumBins;
    private int     histogramStartBin;
    private int     histogramEndBin;
    private String  histogramVariable;
}

