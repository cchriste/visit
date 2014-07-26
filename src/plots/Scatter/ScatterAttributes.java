// ***************************************************************************
//
// Copyright (c) 2000 - 2014, Lawrence Livermore National Security, LLC
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

package llnl.visit.plots;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;
import llnl.visit.ColorAttribute;

// ****************************************************************************
// Class: ScatterAttributes
//
// Purpose:
//    Attributes for the scatter plot
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class ScatterAttributes extends AttributeSubject implements Plugin
{
    private static int ScatterAttributes_numAdditionalAtts = 41;

    // Enum values
    public final static int SCALING_LINEAR = 0;
    public final static int SCALING_LOG = 1;
    public final static int SCALING_SKEW = 2;

    public final static int COLORINGMETHOD_COLORBYFOREGROUNDCOLOR = 0;
    public final static int COLORINGMETHOD_COLORBYSINGLECOLOR = 1;
    public final static int COLORINGMETHOD_COLORBYCOLORTABLE = 2;

    public final static int POINTTYPE_BOX = 0;
    public final static int POINTTYPE_AXIS = 1;
    public final static int POINTTYPE_ICOSAHEDRON = 2;
    public final static int POINTTYPE_OCTAHEDRON = 3;
    public final static int POINTTYPE_TETRAHEDRON = 4;
    public final static int POINTTYPE_SPHEREGEOMETRY = 5;
    public final static int POINTTYPE_POINT = 6;
    public final static int POINTTYPE_SPHERE = 7;

    public final static int VARIABLEROLE_COORDINATE0 = 0;
    public final static int VARIABLEROLE_COORDINATE1 = 1;
    public final static int VARIABLEROLE_COORDINATE2 = 2;
    public final static int VARIABLEROLE_COLOR = 3;
    public final static int VARIABLEROLE_NONE = 4;


    public ScatterAttributes()
    {
        super(ScatterAttributes_numAdditionalAtts);

        var1 = new String("default");
        var1Role = VARIABLEROLE_COORDINATE0;
        var1MinFlag = false;
        var1MaxFlag = false;
        var1Min = 0;
        var1Max = 1;
        var1Scaling = SCALING_LINEAR;
        var1SkewFactor = 1;
        var2Role = VARIABLEROLE_COORDINATE1;
        var2 = new String("default");
        var2MinFlag = false;
        var2MaxFlag = false;
        var2Min = 0;
        var2Max = 1;
        var2Scaling = SCALING_LINEAR;
        var2SkewFactor = 1;
        var3Role = VARIABLEROLE_NONE;
        var3 = new String("default");
        var3MinFlag = false;
        var3MaxFlag = false;
        var3Min = 0;
        var3Max = 1;
        var3Scaling = SCALING_LINEAR;
        var3SkewFactor = 1;
        var4Role = VARIABLEROLE_NONE;
        var4 = new String("default");
        var4MinFlag = false;
        var4MaxFlag = false;
        var4Min = 0;
        var4Max = 1;
        var4Scaling = SCALING_LINEAR;
        var4SkewFactor = 1;
        pointSize = 0.05;
        pointSizePixels = 1;
        pointType = POINTTYPE_POINT;
        scaleCube = true;
        colorType = COLORINGMETHOD_COLORBYFOREGROUNDCOLOR;
        singleColor = new ColorAttribute(255, 0, 0);
        colorTableName = new String("Default");
        invertColorTable = false;
        legendFlag = true;
    }

    public ScatterAttributes(int nMoreFields)
    {
        super(ScatterAttributes_numAdditionalAtts + nMoreFields);

        var1 = new String("default");
        var1Role = VARIABLEROLE_COORDINATE0;
        var1MinFlag = false;
        var1MaxFlag = false;
        var1Min = 0;
        var1Max = 1;
        var1Scaling = SCALING_LINEAR;
        var1SkewFactor = 1;
        var2Role = VARIABLEROLE_COORDINATE1;
        var2 = new String("default");
        var2MinFlag = false;
        var2MaxFlag = false;
        var2Min = 0;
        var2Max = 1;
        var2Scaling = SCALING_LINEAR;
        var2SkewFactor = 1;
        var3Role = VARIABLEROLE_NONE;
        var3 = new String("default");
        var3MinFlag = false;
        var3MaxFlag = false;
        var3Min = 0;
        var3Max = 1;
        var3Scaling = SCALING_LINEAR;
        var3SkewFactor = 1;
        var4Role = VARIABLEROLE_NONE;
        var4 = new String("default");
        var4MinFlag = false;
        var4MaxFlag = false;
        var4Min = 0;
        var4Max = 1;
        var4Scaling = SCALING_LINEAR;
        var4SkewFactor = 1;
        pointSize = 0.05;
        pointSizePixels = 1;
        pointType = POINTTYPE_POINT;
        scaleCube = true;
        colorType = COLORINGMETHOD_COLORBYFOREGROUNDCOLOR;
        singleColor = new ColorAttribute(255, 0, 0);
        colorTableName = new String("Default");
        invertColorTable = false;
        legendFlag = true;
    }

    public ScatterAttributes(ScatterAttributes obj)
    {
        super(ScatterAttributes_numAdditionalAtts);

        var1 = new String(obj.var1);
        var1Role = obj.var1Role;
        var1MinFlag = obj.var1MinFlag;
        var1MaxFlag = obj.var1MaxFlag;
        var1Min = obj.var1Min;
        var1Max = obj.var1Max;
        var1Scaling = obj.var1Scaling;
        var1SkewFactor = obj.var1SkewFactor;
        var2Role = obj.var2Role;
        var2 = new String(obj.var2);
        var2MinFlag = obj.var2MinFlag;
        var2MaxFlag = obj.var2MaxFlag;
        var2Min = obj.var2Min;
        var2Max = obj.var2Max;
        var2Scaling = obj.var2Scaling;
        var2SkewFactor = obj.var2SkewFactor;
        var3Role = obj.var3Role;
        var3 = new String(obj.var3);
        var3MinFlag = obj.var3MinFlag;
        var3MaxFlag = obj.var3MaxFlag;
        var3Min = obj.var3Min;
        var3Max = obj.var3Max;
        var3Scaling = obj.var3Scaling;
        var3SkewFactor = obj.var3SkewFactor;
        var4Role = obj.var4Role;
        var4 = new String(obj.var4);
        var4MinFlag = obj.var4MinFlag;
        var4MaxFlag = obj.var4MaxFlag;
        var4Min = obj.var4Min;
        var4Max = obj.var4Max;
        var4Scaling = obj.var4Scaling;
        var4SkewFactor = obj.var4SkewFactor;
        pointSize = obj.pointSize;
        pointSizePixels = obj.pointSizePixels;
        pointType = obj.pointType;
        scaleCube = obj.scaleCube;
        colorType = obj.colorType;
        singleColor = new ColorAttribute(obj.singleColor);
        colorTableName = new String(obj.colorTableName);
        invertColorTable = obj.invertColorTable;
        legendFlag = obj.legendFlag;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return ScatterAttributes_numAdditionalAtts;
    }

    public boolean equals(ScatterAttributes obj)
    {
        // Create the return value
        return ((var1.equals(obj.var1)) &&
                (var1Role == obj.var1Role) &&
                (var1MinFlag == obj.var1MinFlag) &&
                (var1MaxFlag == obj.var1MaxFlag) &&
                (var1Min == obj.var1Min) &&
                (var1Max == obj.var1Max) &&
                (var1Scaling == obj.var1Scaling) &&
                (var1SkewFactor == obj.var1SkewFactor) &&
                (var2Role == obj.var2Role) &&
                (var2.equals(obj.var2)) &&
                (var2MinFlag == obj.var2MinFlag) &&
                (var2MaxFlag == obj.var2MaxFlag) &&
                (var2Min == obj.var2Min) &&
                (var2Max == obj.var2Max) &&
                (var2Scaling == obj.var2Scaling) &&
                (var2SkewFactor == obj.var2SkewFactor) &&
                (var3Role == obj.var3Role) &&
                (var3.equals(obj.var3)) &&
                (var3MinFlag == obj.var3MinFlag) &&
                (var3MaxFlag == obj.var3MaxFlag) &&
                (var3Min == obj.var3Min) &&
                (var3Max == obj.var3Max) &&
                (var3Scaling == obj.var3Scaling) &&
                (var3SkewFactor == obj.var3SkewFactor) &&
                (var4Role == obj.var4Role) &&
                (var4.equals(obj.var4)) &&
                (var4MinFlag == obj.var4MinFlag) &&
                (var4MaxFlag == obj.var4MaxFlag) &&
                (var4Min == obj.var4Min) &&
                (var4Max == obj.var4Max) &&
                (var4Scaling == obj.var4Scaling) &&
                (var4SkewFactor == obj.var4SkewFactor) &&
                (pointSize == obj.pointSize) &&
                (pointSizePixels == obj.pointSizePixels) &&
                (pointType == obj.pointType) &&
                (scaleCube == obj.scaleCube) &&
                (colorType == obj.colorType) &&
                (singleColor == obj.singleColor) &&
                (colorTableName.equals(obj.colorTableName)) &&
                (invertColorTable == obj.invertColorTable) &&
                (legendFlag == obj.legendFlag));
    }

    public String GetName() { return "Scatter"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetVar1(String var1_)
    {
        var1 = var1_;
        Select(0);
    }

    public void SetVar1Role(int var1Role_)
    {
        var1Role = var1Role_;
        Select(1);
    }

    public void SetVar1MinFlag(boolean var1MinFlag_)
    {
        var1MinFlag = var1MinFlag_;
        Select(2);
    }

    public void SetVar1MaxFlag(boolean var1MaxFlag_)
    {
        var1MaxFlag = var1MaxFlag_;
        Select(3);
    }

    public void SetVar1Min(double var1Min_)
    {
        var1Min = var1Min_;
        Select(4);
    }

    public void SetVar1Max(double var1Max_)
    {
        var1Max = var1Max_;
        Select(5);
    }

    public void SetVar1Scaling(int var1Scaling_)
    {
        var1Scaling = var1Scaling_;
        Select(6);
    }

    public void SetVar1SkewFactor(double var1SkewFactor_)
    {
        var1SkewFactor = var1SkewFactor_;
        Select(7);
    }

    public void SetVar2Role(int var2Role_)
    {
        var2Role = var2Role_;
        Select(8);
    }

    public void SetVar2(String var2_)
    {
        var2 = var2_;
        Select(9);
    }

    public void SetVar2MinFlag(boolean var2MinFlag_)
    {
        var2MinFlag = var2MinFlag_;
        Select(10);
    }

    public void SetVar2MaxFlag(boolean var2MaxFlag_)
    {
        var2MaxFlag = var2MaxFlag_;
        Select(11);
    }

    public void SetVar2Min(double var2Min_)
    {
        var2Min = var2Min_;
        Select(12);
    }

    public void SetVar2Max(double var2Max_)
    {
        var2Max = var2Max_;
        Select(13);
    }

    public void SetVar2Scaling(int var2Scaling_)
    {
        var2Scaling = var2Scaling_;
        Select(14);
    }

    public void SetVar2SkewFactor(double var2SkewFactor_)
    {
        var2SkewFactor = var2SkewFactor_;
        Select(15);
    }

    public void SetVar3Role(int var3Role_)
    {
        var3Role = var3Role_;
        Select(16);
    }

    public void SetVar3(String var3_)
    {
        var3 = var3_;
        Select(17);
    }

    public void SetVar3MinFlag(boolean var3MinFlag_)
    {
        var3MinFlag = var3MinFlag_;
        Select(18);
    }

    public void SetVar3MaxFlag(boolean var3MaxFlag_)
    {
        var3MaxFlag = var3MaxFlag_;
        Select(19);
    }

    public void SetVar3Min(double var3Min_)
    {
        var3Min = var3Min_;
        Select(20);
    }

    public void SetVar3Max(double var3Max_)
    {
        var3Max = var3Max_;
        Select(21);
    }

    public void SetVar3Scaling(int var3Scaling_)
    {
        var3Scaling = var3Scaling_;
        Select(22);
    }

    public void SetVar3SkewFactor(double var3SkewFactor_)
    {
        var3SkewFactor = var3SkewFactor_;
        Select(23);
    }

    public void SetVar4Role(int var4Role_)
    {
        var4Role = var4Role_;
        Select(24);
    }

    public void SetVar4(String var4_)
    {
        var4 = var4_;
        Select(25);
    }

    public void SetVar4MinFlag(boolean var4MinFlag_)
    {
        var4MinFlag = var4MinFlag_;
        Select(26);
    }

    public void SetVar4MaxFlag(boolean var4MaxFlag_)
    {
        var4MaxFlag = var4MaxFlag_;
        Select(27);
    }

    public void SetVar4Min(double var4Min_)
    {
        var4Min = var4Min_;
        Select(28);
    }

    public void SetVar4Max(double var4Max_)
    {
        var4Max = var4Max_;
        Select(29);
    }

    public void SetVar4Scaling(int var4Scaling_)
    {
        var4Scaling = var4Scaling_;
        Select(30);
    }

    public void SetVar4SkewFactor(double var4SkewFactor_)
    {
        var4SkewFactor = var4SkewFactor_;
        Select(31);
    }

    public void SetPointSize(double pointSize_)
    {
        pointSize = pointSize_;
        Select(32);
    }

    public void SetPointSizePixels(int pointSizePixels_)
    {
        pointSizePixels = pointSizePixels_;
        Select(33);
    }

    public void SetPointType(int pointType_)
    {
        pointType = pointType_;
        Select(34);
    }

    public void SetScaleCube(boolean scaleCube_)
    {
        scaleCube = scaleCube_;
        Select(35);
    }

    public void SetColorType(int colorType_)
    {
        colorType = colorType_;
        Select(36);
    }

    public void SetSingleColor(ColorAttribute singleColor_)
    {
        singleColor = singleColor_;
        Select(37);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(38);
    }

    public void SetInvertColorTable(boolean invertColorTable_)
    {
        invertColorTable = invertColorTable_;
        Select(39);
    }

    public void SetLegendFlag(boolean legendFlag_)
    {
        legendFlag = legendFlag_;
        Select(40);
    }

    // Property getting methods
    public String         GetVar1() { return var1; }
    public int            GetVar1Role() { return var1Role; }
    public boolean        GetVar1MinFlag() { return var1MinFlag; }
    public boolean        GetVar1MaxFlag() { return var1MaxFlag; }
    public double         GetVar1Min() { return var1Min; }
    public double         GetVar1Max() { return var1Max; }
    public int            GetVar1Scaling() { return var1Scaling; }
    public double         GetVar1SkewFactor() { return var1SkewFactor; }
    public int            GetVar2Role() { return var2Role; }
    public String         GetVar2() { return var2; }
    public boolean        GetVar2MinFlag() { return var2MinFlag; }
    public boolean        GetVar2MaxFlag() { return var2MaxFlag; }
    public double         GetVar2Min() { return var2Min; }
    public double         GetVar2Max() { return var2Max; }
    public int            GetVar2Scaling() { return var2Scaling; }
    public double         GetVar2SkewFactor() { return var2SkewFactor; }
    public int            GetVar3Role() { return var3Role; }
    public String         GetVar3() { return var3; }
    public boolean        GetVar3MinFlag() { return var3MinFlag; }
    public boolean        GetVar3MaxFlag() { return var3MaxFlag; }
    public double         GetVar3Min() { return var3Min; }
    public double         GetVar3Max() { return var3Max; }
    public int            GetVar3Scaling() { return var3Scaling; }
    public double         GetVar3SkewFactor() { return var3SkewFactor; }
    public int            GetVar4Role() { return var4Role; }
    public String         GetVar4() { return var4; }
    public boolean        GetVar4MinFlag() { return var4MinFlag; }
    public boolean        GetVar4MaxFlag() { return var4MaxFlag; }
    public double         GetVar4Min() { return var4Min; }
    public double         GetVar4Max() { return var4Max; }
    public int            GetVar4Scaling() { return var4Scaling; }
    public double         GetVar4SkewFactor() { return var4SkewFactor; }
    public double         GetPointSize() { return pointSize; }
    public int            GetPointSizePixels() { return pointSizePixels; }
    public int            GetPointType() { return pointType; }
    public boolean        GetScaleCube() { return scaleCube; }
    public int            GetColorType() { return colorType; }
    public ColorAttribute GetSingleColor() { return singleColor; }
    public String         GetColorTableName() { return colorTableName; }
    public boolean        GetInvertColorTable() { return invertColorTable; }
    public boolean        GetLegendFlag() { return legendFlag; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteString(var1);
        if(WriteSelect(1, buf))
            buf.WriteInt(var1Role);
        if(WriteSelect(2, buf))
            buf.WriteBool(var1MinFlag);
        if(WriteSelect(3, buf))
            buf.WriteBool(var1MaxFlag);
        if(WriteSelect(4, buf))
            buf.WriteDouble(var1Min);
        if(WriteSelect(5, buf))
            buf.WriteDouble(var1Max);
        if(WriteSelect(6, buf))
            buf.WriteInt(var1Scaling);
        if(WriteSelect(7, buf))
            buf.WriteDouble(var1SkewFactor);
        if(WriteSelect(8, buf))
            buf.WriteInt(var2Role);
        if(WriteSelect(9, buf))
            buf.WriteString(var2);
        if(WriteSelect(10, buf))
            buf.WriteBool(var2MinFlag);
        if(WriteSelect(11, buf))
            buf.WriteBool(var2MaxFlag);
        if(WriteSelect(12, buf))
            buf.WriteDouble(var2Min);
        if(WriteSelect(13, buf))
            buf.WriteDouble(var2Max);
        if(WriteSelect(14, buf))
            buf.WriteInt(var2Scaling);
        if(WriteSelect(15, buf))
            buf.WriteDouble(var2SkewFactor);
        if(WriteSelect(16, buf))
            buf.WriteInt(var3Role);
        if(WriteSelect(17, buf))
            buf.WriteString(var3);
        if(WriteSelect(18, buf))
            buf.WriteBool(var3MinFlag);
        if(WriteSelect(19, buf))
            buf.WriteBool(var3MaxFlag);
        if(WriteSelect(20, buf))
            buf.WriteDouble(var3Min);
        if(WriteSelect(21, buf))
            buf.WriteDouble(var3Max);
        if(WriteSelect(22, buf))
            buf.WriteInt(var3Scaling);
        if(WriteSelect(23, buf))
            buf.WriteDouble(var3SkewFactor);
        if(WriteSelect(24, buf))
            buf.WriteInt(var4Role);
        if(WriteSelect(25, buf))
            buf.WriteString(var4);
        if(WriteSelect(26, buf))
            buf.WriteBool(var4MinFlag);
        if(WriteSelect(27, buf))
            buf.WriteBool(var4MaxFlag);
        if(WriteSelect(28, buf))
            buf.WriteDouble(var4Min);
        if(WriteSelect(29, buf))
            buf.WriteDouble(var4Max);
        if(WriteSelect(30, buf))
            buf.WriteInt(var4Scaling);
        if(WriteSelect(31, buf))
            buf.WriteDouble(var4SkewFactor);
        if(WriteSelect(32, buf))
            buf.WriteDouble(pointSize);
        if(WriteSelect(33, buf))
            buf.WriteInt(pointSizePixels);
        if(WriteSelect(34, buf))
            buf.WriteInt(pointType);
        if(WriteSelect(35, buf))
            buf.WriteBool(scaleCube);
        if(WriteSelect(36, buf))
            buf.WriteInt(colorType);
        if(WriteSelect(37, buf))
            singleColor.Write(buf);
        if(WriteSelect(38, buf))
            buf.WriteString(colorTableName);
        if(WriteSelect(39, buf))
            buf.WriteBool(invertColorTable);
        if(WriteSelect(40, buf))
            buf.WriteBool(legendFlag);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetVar1(buf.ReadString());
            break;
        case 1:
            SetVar1Role(buf.ReadInt());
            break;
        case 2:
            SetVar1MinFlag(buf.ReadBool());
            break;
        case 3:
            SetVar1MaxFlag(buf.ReadBool());
            break;
        case 4:
            SetVar1Min(buf.ReadDouble());
            break;
        case 5:
            SetVar1Max(buf.ReadDouble());
            break;
        case 6:
            SetVar1Scaling(buf.ReadInt());
            break;
        case 7:
            SetVar1SkewFactor(buf.ReadDouble());
            break;
        case 8:
            SetVar2Role(buf.ReadInt());
            break;
        case 9:
            SetVar2(buf.ReadString());
            break;
        case 10:
            SetVar2MinFlag(buf.ReadBool());
            break;
        case 11:
            SetVar2MaxFlag(buf.ReadBool());
            break;
        case 12:
            SetVar2Min(buf.ReadDouble());
            break;
        case 13:
            SetVar2Max(buf.ReadDouble());
            break;
        case 14:
            SetVar2Scaling(buf.ReadInt());
            break;
        case 15:
            SetVar2SkewFactor(buf.ReadDouble());
            break;
        case 16:
            SetVar3Role(buf.ReadInt());
            break;
        case 17:
            SetVar3(buf.ReadString());
            break;
        case 18:
            SetVar3MinFlag(buf.ReadBool());
            break;
        case 19:
            SetVar3MaxFlag(buf.ReadBool());
            break;
        case 20:
            SetVar3Min(buf.ReadDouble());
            break;
        case 21:
            SetVar3Max(buf.ReadDouble());
            break;
        case 22:
            SetVar3Scaling(buf.ReadInt());
            break;
        case 23:
            SetVar3SkewFactor(buf.ReadDouble());
            break;
        case 24:
            SetVar4Role(buf.ReadInt());
            break;
        case 25:
            SetVar4(buf.ReadString());
            break;
        case 26:
            SetVar4MinFlag(buf.ReadBool());
            break;
        case 27:
            SetVar4MaxFlag(buf.ReadBool());
            break;
        case 28:
            SetVar4Min(buf.ReadDouble());
            break;
        case 29:
            SetVar4Max(buf.ReadDouble());
            break;
        case 30:
            SetVar4Scaling(buf.ReadInt());
            break;
        case 31:
            SetVar4SkewFactor(buf.ReadDouble());
            break;
        case 32:
            SetPointSize(buf.ReadDouble());
            break;
        case 33:
            SetPointSizePixels(buf.ReadInt());
            break;
        case 34:
            SetPointType(buf.ReadInt());
            break;
        case 35:
            SetScaleCube(buf.ReadBool());
            break;
        case 36:
            SetColorType(buf.ReadInt());
            break;
        case 37:
            singleColor.Read(buf);
            Select(37);
            break;
        case 38:
            SetColorTableName(buf.ReadString());
            break;
        case 39:
            SetInvertColorTable(buf.ReadBool());
            break;
        case 40:
            SetLegendFlag(buf.ReadBool());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringToString("var1", var1, indent) + "\n";
        str = str + indent + "var1Role = ";
        if(var1Role == VARIABLEROLE_COORDINATE0)
            str = str + "VARIABLEROLE_COORDINATE0";
        if(var1Role == VARIABLEROLE_COORDINATE1)
            str = str + "VARIABLEROLE_COORDINATE1";
        if(var1Role == VARIABLEROLE_COORDINATE2)
            str = str + "VARIABLEROLE_COORDINATE2";
        if(var1Role == VARIABLEROLE_COLOR)
            str = str + "VARIABLEROLE_COLOR";
        if(var1Role == VARIABLEROLE_NONE)
            str = str + "VARIABLEROLE_NONE";
        str = str + "\n";
        str = str + boolToString("var1MinFlag", var1MinFlag, indent) + "\n";
        str = str + boolToString("var1MaxFlag", var1MaxFlag, indent) + "\n";
        str = str + doubleToString("var1Min", var1Min, indent) + "\n";
        str = str + doubleToString("var1Max", var1Max, indent) + "\n";
        str = str + indent + "var1Scaling = ";
        if(var1Scaling == SCALING_LINEAR)
            str = str + "SCALING_LINEAR";
        if(var1Scaling == SCALING_LOG)
            str = str + "SCALING_LOG";
        if(var1Scaling == SCALING_SKEW)
            str = str + "SCALING_SKEW";
        str = str + "\n";
        str = str + doubleToString("var1SkewFactor", var1SkewFactor, indent) + "\n";
        str = str + indent + "var2Role = ";
        if(var2Role == VARIABLEROLE_COORDINATE0)
            str = str + "VARIABLEROLE_COORDINATE0";
        if(var2Role == VARIABLEROLE_COORDINATE1)
            str = str + "VARIABLEROLE_COORDINATE1";
        if(var2Role == VARIABLEROLE_COORDINATE2)
            str = str + "VARIABLEROLE_COORDINATE2";
        if(var2Role == VARIABLEROLE_COLOR)
            str = str + "VARIABLEROLE_COLOR";
        if(var2Role == VARIABLEROLE_NONE)
            str = str + "VARIABLEROLE_NONE";
        str = str + "\n";
        str = str + stringToString("var2", var2, indent) + "\n";
        str = str + boolToString("var2MinFlag", var2MinFlag, indent) + "\n";
        str = str + boolToString("var2MaxFlag", var2MaxFlag, indent) + "\n";
        str = str + doubleToString("var2Min", var2Min, indent) + "\n";
        str = str + doubleToString("var2Max", var2Max, indent) + "\n";
        str = str + indent + "var2Scaling = ";
        if(var2Scaling == SCALING_LINEAR)
            str = str + "SCALING_LINEAR";
        if(var2Scaling == SCALING_LOG)
            str = str + "SCALING_LOG";
        if(var2Scaling == SCALING_SKEW)
            str = str + "SCALING_SKEW";
        str = str + "\n";
        str = str + doubleToString("var2SkewFactor", var2SkewFactor, indent) + "\n";
        str = str + indent + "var3Role = ";
        if(var3Role == VARIABLEROLE_COORDINATE0)
            str = str + "VARIABLEROLE_COORDINATE0";
        if(var3Role == VARIABLEROLE_COORDINATE1)
            str = str + "VARIABLEROLE_COORDINATE1";
        if(var3Role == VARIABLEROLE_COORDINATE2)
            str = str + "VARIABLEROLE_COORDINATE2";
        if(var3Role == VARIABLEROLE_COLOR)
            str = str + "VARIABLEROLE_COLOR";
        if(var3Role == VARIABLEROLE_NONE)
            str = str + "VARIABLEROLE_NONE";
        str = str + "\n";
        str = str + stringToString("var3", var3, indent) + "\n";
        str = str + boolToString("var3MinFlag", var3MinFlag, indent) + "\n";
        str = str + boolToString("var3MaxFlag", var3MaxFlag, indent) + "\n";
        str = str + doubleToString("var3Min", var3Min, indent) + "\n";
        str = str + doubleToString("var3Max", var3Max, indent) + "\n";
        str = str + indent + "var3Scaling = ";
        if(var3Scaling == SCALING_LINEAR)
            str = str + "SCALING_LINEAR";
        if(var3Scaling == SCALING_LOG)
            str = str + "SCALING_LOG";
        if(var3Scaling == SCALING_SKEW)
            str = str + "SCALING_SKEW";
        str = str + "\n";
        str = str + doubleToString("var3SkewFactor", var3SkewFactor, indent) + "\n";
        str = str + indent + "var4Role = ";
        if(var4Role == VARIABLEROLE_COORDINATE0)
            str = str + "VARIABLEROLE_COORDINATE0";
        if(var4Role == VARIABLEROLE_COORDINATE1)
            str = str + "VARIABLEROLE_COORDINATE1";
        if(var4Role == VARIABLEROLE_COORDINATE2)
            str = str + "VARIABLEROLE_COORDINATE2";
        if(var4Role == VARIABLEROLE_COLOR)
            str = str + "VARIABLEROLE_COLOR";
        if(var4Role == VARIABLEROLE_NONE)
            str = str + "VARIABLEROLE_NONE";
        str = str + "\n";
        str = str + stringToString("var4", var4, indent) + "\n";
        str = str + boolToString("var4MinFlag", var4MinFlag, indent) + "\n";
        str = str + boolToString("var4MaxFlag", var4MaxFlag, indent) + "\n";
        str = str + doubleToString("var4Min", var4Min, indent) + "\n";
        str = str + doubleToString("var4Max", var4Max, indent) + "\n";
        str = str + indent + "var4Scaling = ";
        if(var4Scaling == SCALING_LINEAR)
            str = str + "SCALING_LINEAR";
        if(var4Scaling == SCALING_LOG)
            str = str + "SCALING_LOG";
        if(var4Scaling == SCALING_SKEW)
            str = str + "SCALING_SKEW";
        str = str + "\n";
        str = str + doubleToString("var4SkewFactor", var4SkewFactor, indent) + "\n";
        str = str + doubleToString("pointSize", pointSize, indent) + "\n";
        str = str + intToString("pointSizePixels", pointSizePixels, indent) + "\n";
        str = str + indent + "pointType = ";
        if(pointType == POINTTYPE_BOX)
            str = str + "POINTTYPE_BOX";
        if(pointType == POINTTYPE_AXIS)
            str = str + "POINTTYPE_AXIS";
        if(pointType == POINTTYPE_ICOSAHEDRON)
            str = str + "POINTTYPE_ICOSAHEDRON";
        if(pointType == POINTTYPE_OCTAHEDRON)
            str = str + "POINTTYPE_OCTAHEDRON";
        if(pointType == POINTTYPE_TETRAHEDRON)
            str = str + "POINTTYPE_TETRAHEDRON";
        if(pointType == POINTTYPE_SPHEREGEOMETRY)
            str = str + "POINTTYPE_SPHEREGEOMETRY";
        if(pointType == POINTTYPE_POINT)
            str = str + "POINTTYPE_POINT";
        if(pointType == POINTTYPE_SPHERE)
            str = str + "POINTTYPE_SPHERE";
        str = str + "\n";
        str = str + boolToString("scaleCube", scaleCube, indent) + "\n";
        str = str + indent + "colorType = ";
        if(colorType == COLORINGMETHOD_COLORBYFOREGROUNDCOLOR)
            str = str + "COLORINGMETHOD_COLORBYFOREGROUNDCOLOR";
        if(colorType == COLORINGMETHOD_COLORBYSINGLECOLOR)
            str = str + "COLORINGMETHOD_COLORBYSINGLECOLOR";
        if(colorType == COLORINGMETHOD_COLORBYCOLORTABLE)
            str = str + "COLORINGMETHOD_COLORBYCOLORTABLE";
        str = str + "\n";
        str = str + indent + "singleColor = {" + singleColor.Red() + ", " + singleColor.Green() + ", " + singleColor.Blue() + ", " + singleColor.Alpha() + "}\n";
        str = str + stringToString("colorTableName", colorTableName, indent) + "\n";
        str = str + boolToString("invertColorTable", invertColorTable, indent) + "\n";
        str = str + boolToString("legendFlag", legendFlag, indent) + "\n";
        return str;
    }


    // Attributes
    private String         var1;
    private int            var1Role;
    private boolean        var1MinFlag;
    private boolean        var1MaxFlag;
    private double         var1Min;
    private double         var1Max;
    private int            var1Scaling;
    private double         var1SkewFactor;
    private int            var2Role;
    private String         var2;
    private boolean        var2MinFlag;
    private boolean        var2MaxFlag;
    private double         var2Min;
    private double         var2Max;
    private int            var2Scaling;
    private double         var2SkewFactor;
    private int            var3Role;
    private String         var3;
    private boolean        var3MinFlag;
    private boolean        var3MaxFlag;
    private double         var3Min;
    private double         var3Max;
    private int            var3Scaling;
    private double         var3SkewFactor;
    private int            var4Role;
    private String         var4;
    private boolean        var4MinFlag;
    private boolean        var4MaxFlag;
    private double         var4Min;
    private double         var4Max;
    private int            var4Scaling;
    private double         var4SkewFactor;
    private double         pointSize;
    private int            pointSizePixels;
    private int            pointType;
    private boolean        scaleCube;
    private int            colorType;
    private ColorAttribute singleColor;
    private String         colorTableName;
    private boolean        invertColorTable;
    private boolean        legendFlag;
}

