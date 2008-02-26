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


// ****************************************************************************
// Class: avtSymmetricTensorMetaData
//
// Purpose:
//    Contains symmetricTensor metadata attributes
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Mon Feb 25 15:10:59 PST 2008
//
// Modifications:
//   
// ****************************************************************************

public class avtSymmetricTensorMetaData extends AttributeSubject
{
    public avtSymmetricTensorMetaData()
    {
        super(8);

        name = new String("symmetricTensor");
        originalName = new String("");
        validVariable = true;
        meshName = new String("mesh");
        centering = 0;
        hasUnits = false;
        units = new String("");
        dim = 0;
    }

    public avtSymmetricTensorMetaData(avtSymmetricTensorMetaData obj)
    {
        super(8);

        name = new String(obj.name);
        originalName = new String(obj.originalName);
        validVariable = obj.validVariable;
        meshName = new String(obj.meshName);
        centering = obj.centering;
        hasUnits = obj.hasUnits;
        units = new String(obj.units);
        dim = obj.dim;

        SelectAll();
    }

    public boolean equals(avtSymmetricTensorMetaData obj)
    {
        // Create the return value
        return ((name == obj.name) &&
                (originalName == obj.originalName) &&
                (validVariable == obj.validVariable) &&
                (meshName == obj.meshName) &&
                (centering == obj.centering) &&
                (hasUnits == obj.hasUnits) &&
                (units == obj.units) &&
                (dim == obj.dim));
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

    public void SetDim(int dim_)
    {
        dim = dim_;
        Select(7);
    }

    // Property getting methods
    public String  GetName() { return name; }
    public String  GetOriginalName() { return originalName; }
    public boolean GetValidVariable() { return validVariable; }
    public String  GetMeshName() { return meshName; }
    public int     GetCentering() { return centering; }
    public boolean GetHasUnits() { return hasUnits; }
    public String  GetUnits() { return units; }
    public int     GetDim() { return dim; }

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
            buf.WriteInt(dim);
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
                SetDim(buf.ReadInt());
                break;
            }
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("name", name, indent) + "\n";
        str = str + stringToString("originalName", originalName, indent) + "\n";
        str = str + boolToString("validVariable", validVariable, indent) + "\n";
        str = str + stringToString("meshName", meshName, indent) + "\n";
        str = str + intToString("centering", centering, indent) + "\n";
        str = str + boolToString("hasUnits", hasUnits, indent) + "\n";
        str = str + stringToString("units", units, indent) + "\n";
        str = str + intToString("dim", dim, indent) + "\n";
        return str;
    }


    // Attributes
    private String  name;
    private String  originalName;
    private boolean validVariable;
    private String  meshName;
    private int     centering;
    private boolean hasUnits;
    private String  units;
    private int     dim;
}

