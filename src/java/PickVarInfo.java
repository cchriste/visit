// ***************************************************************************
//
// Copyright (c) 2000 - 2006, The Regents of the University of California
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
// Creation:   Mon Jan 8 20:53:15 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class PickVarInfo extends AttributeSubject
{
    // Enum values
    public final static int CENTERING_NODAL = 0;
    public final static int CENTERING_ZONAL = 1;
    public final static int CENTERING_NONE = 2;


    public PickVarInfo()
    {
        super(13);

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
    }

    public PickVarInfo(PickVarInfo obj)
    {
        super(13);

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

        SelectAll();
    }

    public boolean equals(PickVarInfo obj)
    {
        int i;

        // Create the return value
        return ((variableName == obj.variableName) &&
                (variableType == obj.variableType) &&
                (names == obj.names) &&
                (values == obj.values) &&
                (mixNames == obj.mixNames) &&
                (mixValues == obj.mixValues) &&
                (mixVar == obj.mixVar) &&
                (centering == obj.centering) &&
                (miscMessage == obj.miscMessage) &&
                (numMatsPerZone == obj.numMatsPerZone) &&
                (matNames == obj.matNames) &&
                (numSpecsPerMat == obj.numSpecsPerMat) &&
                (treatAsASCII == obj.treatAsASCII));
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
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
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
            }
        }
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
}

