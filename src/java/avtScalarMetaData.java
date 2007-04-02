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
import java.lang.Integer;

// ****************************************************************************
// Class: avtScalarMetaData
//
// Purpose:
//    Contains scalar metadata attributes
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Mar 14 17:56:08 PST 2007
//
// Modifications:
//   
// ****************************************************************************

public class avtScalarMetaData extends AttributeSubject
{
    public avtScalarMetaData()
    {
        super(14);

        name = new String("scalar");
        originalName = new String("");
        validVariable = true;
        meshName = new String("mesh");
        centering = 0;
        hasUnits = false;
        units = new String("");
        hasDataExtents = false;
        minDataExtents = 0;
        maxDataExtents = 0;
        treatAsASCII = false;
        isEnumeration = false;
        enumNames = new Vector();
        enumValues = new Vector();
    }

    public avtScalarMetaData(avtScalarMetaData obj)
    {
        super(14);

        int i;

        name = new String(obj.name);
        originalName = new String(obj.originalName);
        validVariable = obj.validVariable;
        meshName = new String(obj.meshName);
        centering = obj.centering;
        hasUnits = obj.hasUnits;
        units = new String(obj.units);
        hasDataExtents = obj.hasDataExtents;
        minDataExtents = obj.minDataExtents;
        maxDataExtents = obj.maxDataExtents;
        treatAsASCII = obj.treatAsASCII;
        isEnumeration = obj.isEnumeration;
        enumNames = new Vector(obj.enumNames.size());
        for(i = 0; i < obj.enumNames.size(); ++i)
            enumNames.addElement(new String((String)obj.enumNames.elementAt(i)));

        enumValues = new Vector();
        for(i = 0; i < obj.enumValues.size(); ++i)
        {
            Integer iv = (Integer)obj.enumValues.elementAt(i);
            enumValues.addElement(new Integer(iv.intValue()));
        }

        SelectAll();
    }

    public boolean equals(avtScalarMetaData obj)
    {
        int i;

        // Create the return value
        return ((name == obj.name) &&
                (originalName == obj.originalName) &&
                (validVariable == obj.validVariable) &&
                (meshName == obj.meshName) &&
                (centering == obj.centering) &&
                (hasUnits == obj.hasUnits) &&
                (units == obj.units) &&
                (hasDataExtents == obj.hasDataExtents) &&
                (minDataExtents == obj.minDataExtents) &&
                (maxDataExtents == obj.maxDataExtents) &&
                (treatAsASCII == obj.treatAsASCII) &&
                (isEnumeration == obj.isEnumeration) &&
                (enumNames == obj.enumNames) &&
                (enumValues == obj.enumValues));
    }

    // Property setting methods
    public void SetName(String name_)
    {
        name = name_;
        Select(0);
    }

    public void SetOriginalName(String originalName_)
    {
        originalName = originalName_;
        Select(1);
    }

    public void SetValidVariable(boolean validVariable_)
    {
        validVariable = validVariable_;
        Select(2);
    }

    public void SetMeshName(String meshName_)
    {
        meshName = meshName_;
        Select(3);
    }

    public void SetCentering(int centering_)
    {
        centering = centering_;
        Select(4);
    }

    public void SetHasUnits(boolean hasUnits_)
    {
        hasUnits = hasUnits_;
        Select(5);
    }

    public void SetUnits(String units_)
    {
        units = units_;
        Select(6);
    }

    public void SetHasDataExtents(boolean hasDataExtents_)
    {
        hasDataExtents = hasDataExtents_;
        Select(7);
    }

    public void SetMinDataExtents(double minDataExtents_)
    {
        minDataExtents = minDataExtents_;
        Select(8);
    }

    public void SetMaxDataExtents(double maxDataExtents_)
    {
        maxDataExtents = maxDataExtents_;
        Select(9);
    }

    public void SetTreatAsASCII(boolean treatAsASCII_)
    {
        treatAsASCII = treatAsASCII_;
        Select(10);
    }

    public void SetIsEnumeration(boolean isEnumeration_)
    {
        isEnumeration = isEnumeration_;
        Select(11);
    }

    public void SetEnumNames(Vector enumNames_)
    {
        enumNames = enumNames_;
        Select(12);
    }

    public void SetEnumValues(Vector enumValues_)
    {
        enumValues = enumValues_;
        Select(13);
    }

    // Property getting methods
    public String  GetName() { return name; }
    public String  GetOriginalName() { return originalName; }
    public boolean GetValidVariable() { return validVariable; }
    public String  GetMeshName() { return meshName; }
    public int     GetCentering() { return centering; }
    public boolean GetHasUnits() { return hasUnits; }
    public String  GetUnits() { return units; }
    public boolean GetHasDataExtents() { return hasDataExtents; }
    public double  GetMinDataExtents() { return minDataExtents; }
    public double  GetMaxDataExtents() { return maxDataExtents; }
    public boolean GetTreatAsASCII() { return treatAsASCII; }
    public boolean GetIsEnumeration() { return isEnumeration; }
    public Vector  GetEnumNames() { return enumNames; }
    public Vector  GetEnumValues() { return enumValues; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(name);
        if(WriteSelect(1, buf))
            buf.WriteString(originalName);
        if(WriteSelect(2, buf))
            buf.WriteBool(validVariable);
        if(WriteSelect(3, buf))
            buf.WriteString(meshName);
        if(WriteSelect(4, buf))
            buf.WriteInt(centering);
        if(WriteSelect(5, buf))
            buf.WriteBool(hasUnits);
        if(WriteSelect(6, buf))
            buf.WriteString(units);
        if(WriteSelect(7, buf))
            buf.WriteBool(hasDataExtents);
        if(WriteSelect(8, buf))
            buf.WriteDouble(minDataExtents);
        if(WriteSelect(9, buf))
            buf.WriteDouble(maxDataExtents);
        if(WriteSelect(10, buf))
            buf.WriteBool(treatAsASCII);
        if(WriteSelect(11, buf))
            buf.WriteBool(isEnumeration);
        if(WriteSelect(12, buf))
            buf.WriteStringVector(enumNames);
        if(WriteSelect(13, buf))
            buf.WriteIntVector(enumValues);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetName(buf.ReadString());
                break;
            case 1:
                SetOriginalName(buf.ReadString());
                break;
            case 2:
                SetValidVariable(buf.ReadBool());
                break;
            case 3:
                SetMeshName(buf.ReadString());
                break;
            case 4:
                SetCentering(buf.ReadInt());
                break;
            case 5:
                SetHasUnits(buf.ReadBool());
                break;
            case 6:
                SetUnits(buf.ReadString());
                break;
            case 7:
                SetHasDataExtents(buf.ReadBool());
                break;
            case 8:
                SetMinDataExtents(buf.ReadDouble());
                break;
            case 9:
                SetMaxDataExtents(buf.ReadDouble());
                break;
            case 10:
                SetTreatAsASCII(buf.ReadBool());
                break;
            case 11:
                SetIsEnumeration(buf.ReadBool());
                break;
            case 12:
                SetEnumNames(buf.ReadStringVector());
                break;
            case 13:
                SetEnumValues(buf.ReadIntVector());
                break;
            }
        }
    }


    // Attributes
    private String  name;
    private String  originalName;
    private boolean validVariable;
    private String  meshName;
    private int     centering;
    private boolean hasUnits;
    private String  units;
    private boolean hasDataExtents;
    private double  minDataExtents;
    private double  maxDataExtents;
    private boolean treatAsASCII;
    private boolean isEnumeration;
    private Vector  enumNames; // vector of String objects
    private Vector  enumValues; // vector of Integer objects
}

