// ***************************************************************************
//
// Copyright (c) 2000 - 2011, Lawrence Livermore National Security, LLC
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
import java.lang.Integer;

// ****************************************************************************
// Class: QueryAttributes
//
// Purpose:
//    This class contains attributes used for query.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class QueryAttributes extends AttributeSubject
{
    private static int QueryAttributes_numAdditionalAtts = 19;

    // Enum values
    public final static int ELEMENTTYPE_ZONE = 0;
    public final static int ELEMENTTYPE_NODE = 1;

    public final static int VARTYPE_MESH = 0;
    public final static int VARTYPE_SCALAR = 1;
    public final static int VARTYPE_VECTOR = 2;
    public final static int VARTYPE_TENSOR = 3;
    public final static int VARTYPE_SYMMETRIC_TENSOR = 4;
    public final static int VARTYPE_ARRAY = 5;
    public final static int VARTYPE_LABEL = 6;
    public final static int VARTYPE_MATERIAL = 7;
    public final static int VARTYPE_SPECIES = 8;
    public final static int VARTYPE_CURVE = 9;
    public final static int VARTYPE_UNKNOWN = 10;

    public final static int DATATYPE_ACTUALDATA = 0;
    public final static int DATATYPE_ORIGINALDATA = 1;


    public QueryAttributes()
    {
        super(QueryAttributes_numAdditionalAtts);

        name = new String("");
        variables = new Vector();
        variables.addElement(new String("default"));
        resultsMessage = new String("");
        worldPoint = new double[3];
        worldPoint[0] = 0;
        worldPoint[1] = 0;
        worldPoint[2] = 0;
        domain = -1;
        element = -1;
        resultsValue = new Vector();
        resultsValue.addElement(new Double(0));
        elementType = ELEMENTTYPE_ZONE;
        timeStep = 0;
        varTypes = new Vector();
        dataType = DATATYPE_ACTUALDATA;
        pipeIndex = -1;
        useGlobalId = false;
        xUnits = new String("");
        yUnits = new String("");
        darg1 = new Vector();
        darg1.addElement(new Double(0));
        darg2 = new Vector();
        darg2.addElement(new Double(0));
        floatFormat = new String("%g");
        xmlResult = new String("");
    }

    public QueryAttributes(int nMoreFields)
    {
        super(QueryAttributes_numAdditionalAtts + nMoreFields);

        name = new String("");
        variables = new Vector();
        variables.addElement(new String("default"));
        resultsMessage = new String("");
        worldPoint = new double[3];
        worldPoint[0] = 0;
        worldPoint[1] = 0;
        worldPoint[2] = 0;
        domain = -1;
        element = -1;
        resultsValue = new Vector();
        resultsValue.addElement(new Double(0));
        elementType = ELEMENTTYPE_ZONE;
        timeStep = 0;
        varTypes = new Vector();
        dataType = DATATYPE_ACTUALDATA;
        pipeIndex = -1;
        useGlobalId = false;
        xUnits = new String("");
        yUnits = new String("");
        darg1 = new Vector();
        darg1.addElement(new Double(0));
        darg2 = new Vector();
        darg2.addElement(new Double(0));
        floatFormat = new String("%g");
        xmlResult = new String("");
    }

    public QueryAttributes(QueryAttributes obj)
    {
        super(QueryAttributes_numAdditionalAtts);

        int i;

        name = new String(obj.name);
        variables = new Vector(obj.variables.size());
        for(i = 0; i < obj.variables.size(); ++i)
            variables.addElement(new String((String)obj.variables.elementAt(i)));

        resultsMessage = new String(obj.resultsMessage);
        worldPoint = new double[3];
        worldPoint[0] = obj.worldPoint[0];
        worldPoint[1] = obj.worldPoint[1];
        worldPoint[2] = obj.worldPoint[2];

        domain = obj.domain;
        element = obj.element;
        resultsValue = new Vector(obj.resultsValue.size());
        for(i = 0; i < obj.resultsValue.size(); ++i)
        {
            Double dv = (Double)obj.resultsValue.elementAt(i);
            resultsValue.addElement(new Double(dv.doubleValue()));
        }

        elementType = obj.elementType;
        timeStep = obj.timeStep;
        varTypes = new Vector();
        for(i = 0; i < obj.varTypes.size(); ++i)
        {
            Integer iv = (Integer)obj.varTypes.elementAt(i);
            varTypes.addElement(new Integer(iv.intValue()));
        }
        dataType = obj.dataType;
        pipeIndex = obj.pipeIndex;
        useGlobalId = obj.useGlobalId;
        xUnits = new String(obj.xUnits);
        yUnits = new String(obj.yUnits);
        darg1 = new Vector(obj.darg1.size());
        for(i = 0; i < obj.darg1.size(); ++i)
        {
            Double dv = (Double)obj.darg1.elementAt(i);
            darg1.addElement(new Double(dv.doubleValue()));
        }

        darg2 = new Vector(obj.darg2.size());
        for(i = 0; i < obj.darg2.size(); ++i)
        {
            Double dv = (Double)obj.darg2.elementAt(i);
            darg2.addElement(new Double(dv.doubleValue()));
        }

        floatFormat = new String(obj.floatFormat);
        xmlResult = new String(obj.xmlResult);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return QueryAttributes_numAdditionalAtts;
    }

    public boolean equals(QueryAttributes obj)
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
        // Compare the worldPoint arrays.
        boolean worldPoint_equal = true;
        for(i = 0; i < 3 && worldPoint_equal; ++i)
            worldPoint_equal = (worldPoint[i] == obj.worldPoint[i]);

        // Compare the elements in the resultsValue vector.
        boolean resultsValue_equal = (obj.resultsValue.size() == resultsValue.size());
        for(i = 0; (i < resultsValue.size()) && resultsValue_equal; ++i)
        {
            // Make references to Double from Object.
            Double resultsValue1 = (Double)resultsValue.elementAt(i);
            Double resultsValue2 = (Double)obj.resultsValue.elementAt(i);
            resultsValue_equal = resultsValue1.equals(resultsValue2);
        }
        // Compare the elements in the varTypes vector.
        boolean varTypes_equal = (obj.varTypes.size() == varTypes.size());
        for(i = 0; (i < varTypes.size()) && varTypes_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer varTypes1 = (Integer)varTypes.elementAt(i);
            Integer varTypes2 = (Integer)obj.varTypes.elementAt(i);
            varTypes_equal = varTypes1.equals(varTypes2);
        }
        // Compare the elements in the darg1 vector.
        boolean darg1_equal = (obj.darg1.size() == darg1.size());
        for(i = 0; (i < darg1.size()) && darg1_equal; ++i)
        {
            // Make references to Double from Object.
            Double darg11 = (Double)darg1.elementAt(i);
            Double darg12 = (Double)obj.darg1.elementAt(i);
            darg1_equal = darg11.equals(darg12);
        }
        // Compare the elements in the darg2 vector.
        boolean darg2_equal = (obj.darg2.size() == darg2.size());
        for(i = 0; (i < darg2.size()) && darg2_equal; ++i)
        {
            // Make references to Double from Object.
            Double darg21 = (Double)darg2.elementAt(i);
            Double darg22 = (Double)obj.darg2.elementAt(i);
            darg2_equal = darg21.equals(darg22);
        }
        // Create the return value
        return ((name.equals(obj.name)) &&
                variables_equal &&
                (resultsMessage.equals(obj.resultsMessage)) &&
                worldPoint_equal &&
                (domain == obj.domain) &&
                (element == obj.element) &&
                resultsValue_equal &&
                (elementType == obj.elementType) &&
                (timeStep == obj.timeStep) &&
                varTypes_equal &&
                (dataType == obj.dataType) &&
                (pipeIndex == obj.pipeIndex) &&
                (useGlobalId == obj.useGlobalId) &&
                (xUnits.equals(obj.xUnits)) &&
                (yUnits.equals(obj.yUnits)) &&
                darg1_equal &&
                darg2_equal &&
                (floatFormat.equals(obj.floatFormat)) &&
                (xmlResult.equals(obj.xmlResult)));
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetVariables(Vector variables_)
    {
        variables = variables_;
        Select(1);
    }

    public void SetResultsMessage(String resultsMessage_)
    {
        resultsMessage = resultsMessage_;
        Select(2);
    }

    public void SetWorldPoint(double[] worldPoint_)
    {
        worldPoint[0] = worldPoint_[0];
        worldPoint[1] = worldPoint_[1];
        worldPoint[2] = worldPoint_[2];
        Select(3);
    }

    public void SetWorldPoint(double e0, double e1, double e2)
    {
        worldPoint[0] = e0;
        worldPoint[1] = e1;
        worldPoint[2] = e2;
        Select(3);
    }

    public void SetDomain(int domain_)
    {
        domain = domain_;
        Select(4);
    }

    public void SetElement(int element_)
    {
        element = element_;
        Select(5);
    }

    public void SetResultsValue(Vector resultsValue_)
    {
        resultsValue = resultsValue_;
        Select(6);
    }

    public void SetElementType(int elementType_)
    {
        elementType = elementType_;
        Select(7);
    }

    public void SetTimeStep(int timeStep_)
    {
        timeStep = timeStep_;
        Select(8);
    }

    public void SetVarTypes(Vector varTypes_)
    {
        varTypes = varTypes_;
        Select(9);
    }

    public void SetDataType(int dataType_)
    {
        dataType = dataType_;
        Select(10);
    }

    public void SetPipeIndex(int pipeIndex_)
    {
        pipeIndex = pipeIndex_;
        Select(11);
    }

    public void SetUseGlobalId(boolean useGlobalId_)
    {
        useGlobalId = useGlobalId_;
        Select(12);
    }

    public void SetXUnits(String xUnits_)
    {
        xUnits = xUnits_;
        Select(13);
    }

    public void SetYUnits(String yUnits_)
    {
        yUnits = yUnits_;
        Select(14);
    }

    public void SetDarg1(Vector darg1_)
    {
        darg1 = darg1_;
        Select(15);
    }

    public void SetDarg2(Vector darg2_)
    {
        darg2 = darg2_;
        Select(16);
    }

    public void SetFloatFormat(String floatFormat_)
    {
        floatFormat = floatFormat_;
        Select(17);
    }

    public void SetXmlResult(String xmlResult_)
    {
        xmlResult = xmlResult_;
        Select(18);
    }

    // Property getting methods
    public String   GetName() { return name; }
    public Vector   GetVariables() { return variables; }
    public String   GetResultsMessage() { return resultsMessage; }
    public double[] GetWorldPoint() { return worldPoint; }
    public int      GetDomain() { return domain; }
    public int      GetElement() { return element; }
    public Vector   GetResultsValue() { return resultsValue; }
    public int      GetElementType() { return elementType; }
    public int      GetTimeStep() { return timeStep; }
    public Vector   GetVarTypes() { return varTypes; }
    public int      GetDataType() { return dataType; }
    public int      GetPipeIndex() { return pipeIndex; }
    public boolean  GetUseGlobalId() { return useGlobalId; }
    public String   GetXUnits() { return xUnits; }
    public String   GetYUnits() { return yUnits; }
    public Vector   GetDarg1() { return darg1; }
    public Vector   GetDarg2() { return darg2; }
    public String   GetFloatFormat() { return floatFormat; }
    public String   GetXmlResult() { return xmlResult; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
            buf.WriteStringVector(variables);
        if(WriteSelect(2, buf))
            buf.WriteString(resultsMessage);
        if(WriteSelect(3, buf))
            buf.WriteDoubleArray(worldPoint);
        if(WriteSelect(4, buf))
            buf.WriteInt(domain);
        if(WriteSelect(5, buf))
            buf.WriteInt(element);
        if(WriteSelect(6, buf))
            buf.WriteDoubleVector(resultsValue);
        if(WriteSelect(7, buf))
            buf.WriteInt(elementType);
        if(WriteSelect(8, buf))
            buf.WriteInt(timeStep);
        if(WriteSelect(9, buf))
            buf.WriteIntVector(varTypes);
        if(WriteSelect(10, buf))
            buf.WriteInt(dataType);
        if(WriteSelect(11, buf))
            buf.WriteInt(pipeIndex);
        if(WriteSelect(12, buf))
            buf.WriteBool(useGlobalId);
        if(WriteSelect(13, buf))
            buf.WriteString(xUnits);
        if(WriteSelect(14, buf))
            buf.WriteString(yUnits);
        if(WriteSelect(15, buf))
            buf.WriteDoubleVector(darg1);
        if(WriteSelect(16, buf))
            buf.WriteDoubleVector(darg2);
        if(WriteSelect(17, buf))
            buf.WriteString(floatFormat);
        if(WriteSelect(18, buf))
            buf.WriteString(xmlResult);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetName(buf.ReadString());
            break;
        case 1:
            SetVariables(buf.ReadStringVector());
            break;
        case 2:
            SetResultsMessage(buf.ReadString());
            break;
        case 3:
            SetWorldPoint(buf.ReadDoubleArray());
            break;
        case 4:
            SetDomain(buf.ReadInt());
            break;
        case 5:
            SetElement(buf.ReadInt());
            break;
        case 6:
            SetResultsValue(buf.ReadDoubleVector());
            break;
        case 7:
            SetElementType(buf.ReadInt());
            break;
        case 8:
            SetTimeStep(buf.ReadInt());
            break;
        case 9:
            SetVarTypes(buf.ReadIntVector());
            break;
        case 10:
            SetDataType(buf.ReadInt());
            break;
        case 11:
            SetPipeIndex(buf.ReadInt());
            break;
        case 12:
            SetUseGlobalId(buf.ReadBool());
            break;
        case 13:
            SetXUnits(buf.ReadString());
            break;
        case 14:
            SetYUnits(buf.ReadString());
            break;
        case 15:
            SetDarg1(buf.ReadDoubleVector());
            break;
        case 16:
            SetDarg2(buf.ReadDoubleVector());
            break;
        case 17:
            SetFloatFormat(buf.ReadString());
            break;
        case 18:
            SetXmlResult(buf.ReadString());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("name", name, indent) + "\n";
        str = str + stringVectorToString("variables", variables, indent) + "\n";
        str = str + stringToString("resultsMessage", resultsMessage, indent) + "\n";
        str = str + doubleArrayToString("worldPoint", worldPoint, indent) + "\n";
        str = str + intToString("domain", domain, indent) + "\n";
        str = str + intToString("element", element, indent) + "\n";
        str = str + doubleVectorToString("resultsValue", resultsValue, indent) + "\n";
        str = str + indent + "elementType = ";
        if(elementType == ELEMENTTYPE_ZONE)
            str = str + "ELEMENTTYPE_ZONE";
        if(elementType == ELEMENTTYPE_NODE)
            str = str + "ELEMENTTYPE_NODE";
        str = str + "\n";
        str = str + intToString("timeStep", timeStep, indent) + "\n";
        str = str + intVectorToString("varTypes", varTypes, indent) + "\n";
        str = str + indent + "dataType = ";
        if(dataType == DATATYPE_ACTUALDATA)
            str = str + "DATATYPE_ACTUALDATA";
        if(dataType == DATATYPE_ORIGINALDATA)
            str = str + "DATATYPE_ORIGINALDATA";
        str = str + "\n";
        str = str + intToString("pipeIndex", pipeIndex, indent) + "\n";
        str = str + boolToString("useGlobalId", useGlobalId, indent) + "\n";
        str = str + stringToString("xUnits", xUnits, indent) + "\n";
        str = str + stringToString("yUnits", yUnits, indent) + "\n";
        str = str + doubleVectorToString("darg1", darg1, indent) + "\n";
        str = str + doubleVectorToString("darg2", darg2, indent) + "\n";
        str = str + stringToString("floatFormat", floatFormat, indent) + "\n";
        str = str + stringToString("xmlResult", xmlResult, indent) + "\n";
        return str;
    }


    // Attributes
    private String   name;
    private Vector   variables; // vector of String objects
    private String   resultsMessage;
    private double[] worldPoint;
    private int      domain;
    private int      element;
    private Vector   resultsValue; // vector of Double objects
    private int      elementType;
    private int      timeStep;
    private Vector   varTypes; // vector of Integer objects
    private int      dataType;
    private int      pipeIndex;
    private boolean  useGlobalId;
    private String   xUnits;
    private String   yUnits;
    private Vector   darg1; // vector of Double objects
    private Vector   darg2; // vector of Double objects
    private String   floatFormat;
    private String   xmlResult;
}

