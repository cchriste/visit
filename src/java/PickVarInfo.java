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
// Class: PickVarInfo
//
// Purpose:
//    This class contains PickVarInfo.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class PickVarInfo extends AttributeSubject
{
    private static int PickVarInfo_numAdditionalAtts = 14;

    // Enum values
    public final static int CENTERING_NODAL = 0;
    public final static int CENTERING_ZONAL = 1;
    public final static int CENTERING_NONE = 2;


    public PickVarInfo()
    {
        super(PickVarInfo_numAdditionalAtts);

        variableName = new String("");
        variableType = new String("");
        names = new Vector();
        values = new Vector();
        mixNames = new Vector();
        mixValues = new Vector();
        mixVar = false;
        centering = CENTERING_NONE;
        miscMessage = new String("");
        numMatsPerZone = new Vector();
        matNames = new Vector();
        numSpecsPerMat = new Vector();
        treatAsASCII = false;
        floatFormat = new String("%g");
    }

    public PickVarInfo(int nMoreFields)
    {
        super(PickVarInfo_numAdditionalAtts + nMoreFields);

        variableName = new String("");
        variableType = new String("");
        names = new Vector();
        values = new Vector();
        mixNames = new Vector();
        mixValues = new Vector();
        mixVar = false;
        centering = CENTERING_NONE;
        miscMessage = new String("");
        numMatsPerZone = new Vector();
        matNames = new Vector();
        numSpecsPerMat = new Vector();
        treatAsASCII = false;
        floatFormat = new String("%g");
    }

    public PickVarInfo(PickVarInfo obj)
    {
        super(PickVarInfo_numAdditionalAtts);

        int i;

        variableName = new String(obj.variableName);
        variableType = new String(obj.variableType);
        names = new Vector(obj.names.size());
        for(i = 0; i < obj.names.size(); ++i)
            names.addElement(new String((String)obj.names.elementAt(i)));

        values = new Vector(obj.values.size());
        for(i = 0; i < obj.values.size(); ++i)
        {
            Double dv = (Double)obj.values.elementAt(i);
            values.addElement(new Double(dv.doubleValue()));
        }

        mixNames = new Vector(obj.mixNames.size());
        for(i = 0; i < obj.mixNames.size(); ++i)
            mixNames.addElement(new String((String)obj.mixNames.elementAt(i)));

        mixValues = new Vector(obj.mixValues.size());
        for(i = 0; i < obj.mixValues.size(); ++i)
        {
            Double dv = (Double)obj.mixValues.elementAt(i);
            mixValues.addElement(new Double(dv.doubleValue()));
        }

        mixVar = obj.mixVar;
        centering = obj.centering;
        miscMessage = new String(obj.miscMessage);
        numMatsPerZone = new Vector();
        for(i = 0; i < obj.numMatsPerZone.size(); ++i)
        {
            Integer iv = (Integer)obj.numMatsPerZone.elementAt(i);
            numMatsPerZone.addElement(new Integer(iv.intValue()));
        }
        matNames = new Vector(obj.matNames.size());
        for(i = 0; i < obj.matNames.size(); ++i)
            matNames.addElement(new String((String)obj.matNames.elementAt(i)));

        numSpecsPerMat = new Vector();
        for(i = 0; i < obj.numSpecsPerMat.size(); ++i)
        {
            Integer iv = (Integer)obj.numSpecsPerMat.elementAt(i);
            numSpecsPerMat.addElement(new Integer(iv.intValue()));
        }
        treatAsASCII = obj.treatAsASCII;
        floatFormat = new String(obj.floatFormat);

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return PickVarInfo_numAdditionalAtts;
    }

    public boolean equals(PickVarInfo obj)
    {
        int i;

        // Compare the elements in the names vector.
        boolean names_equal = (obj.names.size() == names.size());
        for(i = 0; (i < names.size()) && names_equal; ++i)
        {
            // Make references to String from Object.
            String names1 = (String)names.elementAt(i);
            String names2 = (String)obj.names.elementAt(i);
            names_equal = names1.equals(names2);
        }
        // Compare the elements in the values vector.
        boolean values_equal = (obj.values.size() == values.size());
        for(i = 0; (i < values.size()) && values_equal; ++i)
        {
            // Make references to Double from Object.
            Double values1 = (Double)values.elementAt(i);
            Double values2 = (Double)obj.values.elementAt(i);
            values_equal = values1.equals(values2);
        }
        // Compare the elements in the mixNames vector.
        boolean mixNames_equal = (obj.mixNames.size() == mixNames.size());
        for(i = 0; (i < mixNames.size()) && mixNames_equal; ++i)
        {
            // Make references to String from Object.
            String mixNames1 = (String)mixNames.elementAt(i);
            String mixNames2 = (String)obj.mixNames.elementAt(i);
            mixNames_equal = mixNames1.equals(mixNames2);
        }
        // Compare the elements in the mixValues vector.
        boolean mixValues_equal = (obj.mixValues.size() == mixValues.size());
        for(i = 0; (i < mixValues.size()) && mixValues_equal; ++i)
        {
            // Make references to Double from Object.
            Double mixValues1 = (Double)mixValues.elementAt(i);
            Double mixValues2 = (Double)obj.mixValues.elementAt(i);
            mixValues_equal = mixValues1.equals(mixValues2);
        }
        // Compare the elements in the numMatsPerZone vector.
        boolean numMatsPerZone_equal = (obj.numMatsPerZone.size() == numMatsPerZone.size());
        for(i = 0; (i < numMatsPerZone.size()) && numMatsPerZone_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer numMatsPerZone1 = (Integer)numMatsPerZone.elementAt(i);
            Integer numMatsPerZone2 = (Integer)obj.numMatsPerZone.elementAt(i);
            numMatsPerZone_equal = numMatsPerZone1.equals(numMatsPerZone2);
        }
        // Compare the elements in the matNames vector.
        boolean matNames_equal = (obj.matNames.size() == matNames.size());
        for(i = 0; (i < matNames.size()) && matNames_equal; ++i)
        {
            // Make references to String from Object.
            String matNames1 = (String)matNames.elementAt(i);
            String matNames2 = (String)obj.matNames.elementAt(i);
            matNames_equal = matNames1.equals(matNames2);
        }
        // Compare the elements in the numSpecsPerMat vector.
        boolean numSpecsPerMat_equal = (obj.numSpecsPerMat.size() == numSpecsPerMat.size());
        for(i = 0; (i < numSpecsPerMat.size()) && numSpecsPerMat_equal; ++i)
        {
            // Make references to Integer from Object.
            Integer numSpecsPerMat1 = (Integer)numSpecsPerMat.elementAt(i);
            Integer numSpecsPerMat2 = (Integer)obj.numSpecsPerMat.elementAt(i);
            numSpecsPerMat_equal = numSpecsPerMat1.equals(numSpecsPerMat2);
        }
        // Create the return value
        return ((variableName.equals(obj.variableName)) &&
                (variableType.equals(obj.variableType)) &&
                names_equal &&
                values_equal &&
                mixNames_equal &&
                mixValues_equal &&
                (mixVar == obj.mixVar) &&
                (centering == obj.centering) &&
                (miscMessage.equals(obj.miscMessage)) &&
                numMatsPerZone_equal &&
                matNames_equal &&
                numSpecsPerMat_equal &&
                (treatAsASCII == obj.treatAsASCII) &&
                (floatFormat.equals(obj.floatFormat)));
    }

    // Property setting methods
    public void SetVariableName(String variableName_)
    {
        variableName = variableName_;
        Select(0);
    }

    public void SetVariableType(String variableType_)
    {
        variableType = variableType_;
        Select(1);
    }

    public void SetNames(Vector names_)
    {
        names = names_;
        Select(2);
    }

    public void SetValues(Vector values_)
    {
        values = values_;
        Select(3);
    }

    public void SetMixNames(Vector mixNames_)
    {
        mixNames = mixNames_;
        Select(4);
    }

    public void SetMixValues(Vector mixValues_)
    {
        mixValues = mixValues_;
        Select(5);
    }

    public void SetMixVar(boolean mixVar_)
    {
        mixVar = mixVar_;
        Select(6);
    }

    public void SetCentering(int centering_)
    {
        centering = centering_;
        Select(7);
    }

    public void SetMiscMessage(String miscMessage_)
    {
        miscMessage = miscMessage_;
        Select(8);
    }

    public void SetNumMatsPerZone(Vector numMatsPerZone_)
    {
        numMatsPerZone = numMatsPerZone_;
        Select(9);
    }

    public void SetMatNames(Vector matNames_)
    {
        matNames = matNames_;
        Select(10);
    }

    public void SetNumSpecsPerMat(Vector numSpecsPerMat_)
    {
        numSpecsPerMat = numSpecsPerMat_;
        Select(11);
    }

    public void SetTreatAsASCII(boolean treatAsASCII_)
    {
        treatAsASCII = treatAsASCII_;
        Select(12);
    }

    public void SetFloatFormat(String floatFormat_)
    {
        floatFormat = floatFormat_;
        Select(13);
    }

    // Property getting methods
    public String  GetVariableName() { return variableName; }
    public String  GetVariableType() { return variableType; }
    public Vector  GetNames() { return names; }
    public Vector  GetValues() { return values; }
    public Vector  GetMixNames() { return mixNames; }
    public Vector  GetMixValues() { return mixValues; }
    public boolean GetMixVar() { return mixVar; }
    public int     GetCentering() { return centering; }
    public String  GetMiscMessage() { return miscMessage; }
    public Vector  GetNumMatsPerZone() { return numMatsPerZone; }
    public Vector  GetMatNames() { return matNames; }
    public Vector  GetNumSpecsPerMat() { return numSpecsPerMat; }
    public boolean GetTreatAsASCII() { return treatAsASCII; }
    public String  GetFloatFormat() { return floatFormat; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(variableName);
        if(WriteSelect(1, buf))
            buf.WriteString(variableType);
        if(WriteSelect(2, buf))
            buf.WriteStringVector(names);
        if(WriteSelect(3, buf))
            buf.WriteDoubleVector(values);
        if(WriteSelect(4, buf))
            buf.WriteStringVector(mixNames);
        if(WriteSelect(5, buf))
            buf.WriteDoubleVector(mixValues);
        if(WriteSelect(6, buf))
            buf.WriteBool(mixVar);
        if(WriteSelect(7, buf))
            buf.WriteInt(centering);
        if(WriteSelect(8, buf))
            buf.WriteString(miscMessage);
        if(WriteSelect(9, buf))
            buf.WriteIntVector(numMatsPerZone);
        if(WriteSelect(10, buf))
            buf.WriteStringVector(matNames);
        if(WriteSelect(11, buf))
            buf.WriteIntVector(numSpecsPerMat);
        if(WriteSelect(12, buf))
            buf.WriteBool(treatAsASCII);
        if(WriteSelect(13, buf))
            buf.WriteString(floatFormat);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetVariableName(buf.ReadString());
            break;
        case 1:
            SetVariableType(buf.ReadString());
            break;
        case 2:
            SetNames(buf.ReadStringVector());
            break;
        case 3:
            SetValues(buf.ReadDoubleVector());
            break;
        case 4:
            SetMixNames(buf.ReadStringVector());
            break;
        case 5:
            SetMixValues(buf.ReadDoubleVector());
            break;
        case 6:
            SetMixVar(buf.ReadBool());
            break;
        case 7:
            SetCentering(buf.ReadInt());
            break;
        case 8:
            SetMiscMessage(buf.ReadString());
            break;
        case 9:
            SetNumMatsPerZone(buf.ReadIntVector());
            break;
        case 10:
            SetMatNames(buf.ReadStringVector());
            break;
        case 11:
            SetNumSpecsPerMat(buf.ReadIntVector());
            break;
        case 12:
            SetTreatAsASCII(buf.ReadBool());
            break;
        case 13:
            SetFloatFormat(buf.ReadString());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("variableName", variableName, indent) + "\n";
        str = str + stringToString("variableType", variableType, indent) + "\n";
        str = str + stringVectorToString("names", names, indent) + "\n";
        str = str + doubleVectorToString("values", values, indent) + "\n";
        str = str + stringVectorToString("mixNames", mixNames, indent) + "\n";
        str = str + doubleVectorToString("mixValues", mixValues, indent) + "\n";
        str = str + boolToString("mixVar", mixVar, indent) + "\n";
        str = str + indent + "centering = ";
        if(centering == CENTERING_NODAL)
            str = str + "CENTERING_NODAL";
        if(centering == CENTERING_ZONAL)
            str = str + "CENTERING_ZONAL";
        if(centering == CENTERING_NONE)
            str = str + "CENTERING_NONE";
        str = str + "\n";
        str = str + stringToString("miscMessage", miscMessage, indent) + "\n";
        str = str + intVectorToString("numMatsPerZone", numMatsPerZone, indent) + "\n";
        str = str + stringVectorToString("matNames", matNames, indent) + "\n";
        str = str + intVectorToString("numSpecsPerMat", numSpecsPerMat, indent) + "\n";
        str = str + boolToString("treatAsASCII", treatAsASCII, indent) + "\n";
        str = str + stringToString("floatFormat", floatFormat, indent) + "\n";
        return str;
    }


    // Attributes
    private String  variableName;
    private String  variableType;
    private Vector  names; // vector of String objects
    private Vector  values; // vector of Double objects
    private Vector  mixNames; // vector of String objects
    private Vector  mixValues; // vector of Double objects
    private boolean mixVar;
    private int     centering;
    private String  miscMessage;
    private Vector  numMatsPerZone; // vector of Integer objects
    private Vector  matNames; // vector of String objects
    private Vector  numSpecsPerMat; // vector of Integer objects
    private boolean treatAsASCII;
    private String  floatFormat;
}

