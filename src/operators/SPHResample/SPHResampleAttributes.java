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

package llnl.visit.operators;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;

// ****************************************************************************
// Class: SPHResampleAttributes
//
// Purpose:
//    
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class SPHResampleAttributes extends AttributeSubject implements Plugin
{
    private static int SPHResampleAttributes_numAdditionalAtts = 12;

    public SPHResampleAttributes()
    {
        super(SPHResampleAttributes_numAdditionalAtts);

        minX = 0f;
        maxX = 1f;
        xnum = 10;
        minY = 0f;
        maxY = 1f;
        ynum = 10;
        minZ = 0f;
        maxZ = 1f;
        znum = 10;
        tensorSupportVariable = new String("H");
        weightVariable = new String("mass");
        RK = true;
    }

    public SPHResampleAttributes(int nMoreFields)
    {
        super(SPHResampleAttributes_numAdditionalAtts + nMoreFields);

        minX = 0f;
        maxX = 1f;
        xnum = 10;
        minY = 0f;
        maxY = 1f;
        ynum = 10;
        minZ = 0f;
        maxZ = 1f;
        znum = 10;
        tensorSupportVariable = new String("H");
        weightVariable = new String("mass");
        RK = true;
    }

    public SPHResampleAttributes(SPHResampleAttributes obj)
    {
        super(SPHResampleAttributes_numAdditionalAtts);

        minX = obj.minX;
        maxX = obj.maxX;
        xnum = obj.xnum;
        minY = obj.minY;
        maxY = obj.maxY;
        ynum = obj.ynum;
        minZ = obj.minZ;
        maxZ = obj.maxZ;
        znum = obj.znum;
        tensorSupportVariable = new String(obj.tensorSupportVariable);
        weightVariable = new String(obj.weightVariable);
        RK = obj.RK;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return SPHResampleAttributes_numAdditionalAtts;
    }

    public boolean equals(SPHResampleAttributes obj)
    {
        // Create the return value
        return ((minX == obj.minX) &&
                (maxX == obj.maxX) &&
                (xnum == obj.xnum) &&
                (minY == obj.minY) &&
                (maxY == obj.maxY) &&
                (ynum == obj.ynum) &&
                (minZ == obj.minZ) &&
                (maxZ == obj.maxZ) &&
                (znum == obj.znum) &&
                (tensorSupportVariable.equals(obj.tensorSupportVariable)) &&
                (weightVariable.equals(obj.weightVariable)) &&
                (RK == obj.RK));
    }

    public String GetName() { return "SPHResample"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetMinX(float minX_)
    {
        minX = minX_;
        Select(0);
    }

    public void SetMaxX(float maxX_)
    {
        maxX = maxX_;
        Select(1);
    }

    public void SetXnum(int xnum_)
    {
        xnum = xnum_;
        Select(2);
    }

    public void SetMinY(float minY_)
    {
        minY = minY_;
        Select(3);
    }

    public void SetMaxY(float maxY_)
    {
        maxY = maxY_;
        Select(4);
    }

    public void SetYnum(int ynum_)
    {
        ynum = ynum_;
        Select(5);
    }

    public void SetMinZ(float minZ_)
    {
        minZ = minZ_;
        Select(6);
    }

    public void SetMaxZ(float maxZ_)
    {
        maxZ = maxZ_;
        Select(7);
    }

    public void SetZnum(int znum_)
    {
        znum = znum_;
        Select(8);
    }

    public void SetTensorSupportVariable(String tensorSupportVariable_)
    {
        tensorSupportVariable = tensorSupportVariable_;
        Select(9);
    }

    public void SetWeightVariable(String weightVariable_)
    {
        weightVariable = weightVariable_;
        Select(10);
    }

    public void SetRK(boolean RK_)
    {
        RK = RK_;
        Select(11);
    }

    // Property getting methods
    public float   GetMinX() { return minX; }
    public float   GetMaxX() { return maxX; }
    public int     GetXnum() { return xnum; }
    public float   GetMinY() { return minY; }
    public float   GetMaxY() { return maxY; }
    public int     GetYnum() { return ynum; }
    public float   GetMinZ() { return minZ; }
    public float   GetMaxZ() { return maxZ; }
    public int     GetZnum() { return znum; }
    public String  GetTensorSupportVariable() { return tensorSupportVariable; }
    public String  GetWeightVariable() { return weightVariable; }
    public boolean GetRK() { return RK; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteFloat(minX);
        if(WriteSelect(1, buf))
            buf.WriteFloat(maxX);
        if(WriteSelect(2, buf))
            buf.WriteInt(xnum);
        if(WriteSelect(3, buf))
            buf.WriteFloat(minY);
        if(WriteSelect(4, buf))
            buf.WriteFloat(maxY);
        if(WriteSelect(5, buf))
            buf.WriteInt(ynum);
        if(WriteSelect(6, buf))
            buf.WriteFloat(minZ);
        if(WriteSelect(7, buf))
            buf.WriteFloat(maxZ);
        if(WriteSelect(8, buf))
            buf.WriteInt(znum);
        if(WriteSelect(9, buf))
            buf.WriteString(tensorSupportVariable);
        if(WriteSelect(10, buf))
            buf.WriteString(weightVariable);
        if(WriteSelect(11, buf))
            buf.WriteBool(RK);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetMinX(buf.ReadFloat());
            break;
        case 1:
            SetMaxX(buf.ReadFloat());
            break;
        case 2:
            SetXnum(buf.ReadInt());
            break;
        case 3:
            SetMinY(buf.ReadFloat());
            break;
        case 4:
            SetMaxY(buf.ReadFloat());
            break;
        case 5:
            SetYnum(buf.ReadInt());
            break;
        case 6:
            SetMinZ(buf.ReadFloat());
            break;
        case 7:
            SetMaxZ(buf.ReadFloat());
            break;
        case 8:
            SetZnum(buf.ReadInt());
            break;
        case 9:
            SetTensorSupportVariable(buf.ReadString());
            break;
        case 10:
            SetWeightVariable(buf.ReadString());
            break;
        case 11:
            SetRK(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + floatToString("minX", minX, indent) + "\n";
        str = str + floatToString("maxX", maxX, indent) + "\n";
        str = str + intToString("xnum", xnum, indent) + "\n";
        str = str + floatToString("minY", minY, indent) + "\n";
        str = str + floatToString("maxY", maxY, indent) + "\n";
        str = str + intToString("ynum", ynum, indent) + "\n";
        str = str + floatToString("minZ", minZ, indent) + "\n";
        str = str + floatToString("maxZ", maxZ, indent) + "\n";
        str = str + intToString("znum", znum, indent) + "\n";
        str = str + stringToString("tensorSupportVariable", tensorSupportVariable, indent) + "\n";
        str = str + stringToString("weightVariable", weightVariable, indent) + "\n";
        str = str + boolToString("RK", RK, indent) + "\n";
        return str;
    }


    // Attributes
    private float   minX;
    private float   maxX;
    private int     xnum;
    private float   minY;
    private float   maxY;
    private int     ynum;
    private float   minZ;
    private float   maxZ;
    private int     znum;
    private String  tensorSupportVariable;
    private String  weightVariable;
    private boolean RK;
}

