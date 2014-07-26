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
import llnl.visit.ColorAttributeList;
import java.util.Vector;

// ****************************************************************************
// Class: SubsetAttributes
//
// Purpose:
//    This class contains the plot attributes for the subset boundary plot.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

public class SubsetAttributes extends AttributeSubject implements Plugin
{
    private static int SubsetAttributes_numAdditionalAtts = 20;

    // Enum values
    public final static int SUBSET_TYPE_DOMAIN = 0;
    public final static int SUBSET_TYPE_GROUP = 1;
    public final static int SUBSET_TYPE_MATERIAL = 2;
    public final static int SUBSET_TYPE_ENUMSCALAR = 3;
    public final static int SUBSET_TYPE_MESH = 4;
    public final static int SUBSET_TYPE_UNKNOWN = 5;

    public final static int COLORINGMETHOD_COLORBYSINGLECOLOR = 0;
    public final static int COLORINGMETHOD_COLORBYMULTIPLECOLORS = 1;
    public final static int COLORINGMETHOD_COLORBYCOLORTABLE = 2;

    public final static int POINTTYPE_BOX = 0;
    public final static int POINTTYPE_AXIS = 1;
    public final static int POINTTYPE_ICOSAHEDRON = 2;
    public final static int POINTTYPE_OCTAHEDRON = 3;
    public final static int POINTTYPE_TETRAHEDRON = 4;
    public final static int POINTTYPE_SPHEREGEOMETRY = 5;
    public final static int POINTTYPE_POINT = 6;
    public final static int POINTTYPE_SPHERE = 7;


    public SubsetAttributes()
    {
        super(SubsetAttributes_numAdditionalAtts);

        colorType = COLORINGMETHOD_COLORBYMULTIPLECOLORS;
        colorTableName = new String("Default");
        invertColorTable = false;
        filledFlag = true;
        legendFlag = true;
        lineStyle = 0;
        lineWidth = 0;
        singleColor = new ColorAttribute();
        multiColor = new ColorAttributeList();
        subsetNames = new Vector();
        subsetType = SUBSET_TYPE_UNKNOWN;
        opacity = 1;
        wireframe = false;
        drawInternal = false;
        smoothingLevel = 0;
        pointSize = 0.05;
        pointType = POINTTYPE_POINT;
        pointSizeVarEnabled = false;
        pointSizeVar = new String("default");
        pointSizePixels = 2;
    }

    public SubsetAttributes(int nMoreFields)
    {
        super(SubsetAttributes_numAdditionalAtts + nMoreFields);

        colorType = COLORINGMETHOD_COLORBYMULTIPLECOLORS;
        colorTableName = new String("Default");
        invertColorTable = false;
        filledFlag = true;
        legendFlag = true;
        lineStyle = 0;
        lineWidth = 0;
        singleColor = new ColorAttribute();
        multiColor = new ColorAttributeList();
        subsetNames = new Vector();
        subsetType = SUBSET_TYPE_UNKNOWN;
        opacity = 1;
        wireframe = false;
        drawInternal = false;
        smoothingLevel = 0;
        pointSize = 0.05;
        pointType = POINTTYPE_POINT;
        pointSizeVarEnabled = false;
        pointSizeVar = new String("default");
        pointSizePixels = 2;
    }

    public SubsetAttributes(SubsetAttributes obj)
    {
        super(SubsetAttributes_numAdditionalAtts);

        int i;

        colorType = obj.colorType;
        colorTableName = new String(obj.colorTableName);
        invertColorTable = obj.invertColorTable;
        filledFlag = obj.filledFlag;
        legendFlag = obj.legendFlag;
        lineStyle = obj.lineStyle;
        lineWidth = obj.lineWidth;
        singleColor = new ColorAttribute(obj.singleColor);
        multiColor = new ColorAttributeList(obj.multiColor);
        subsetNames = new Vector(obj.subsetNames.size());
        for(i = 0; i < obj.subsetNames.size(); ++i)
            subsetNames.addElement(new String((String)obj.subsetNames.elementAt(i)));

        subsetType = obj.subsetType;
        opacity = obj.opacity;
        wireframe = obj.wireframe;
        drawInternal = obj.drawInternal;
        smoothingLevel = obj.smoothingLevel;
        pointSize = obj.pointSize;
        pointType = obj.pointType;
        pointSizeVarEnabled = obj.pointSizeVarEnabled;
        pointSizeVar = new String(obj.pointSizeVar);
        pointSizePixels = obj.pointSizePixels;

        SelectAll();
    }

    public int Offset()
    {
        return super.Offset() + super.GetNumAdditionalAttributes();
    }

    public int GetNumAdditionalAttributes()
    {
        return SubsetAttributes_numAdditionalAtts;
    }

    public boolean equals(SubsetAttributes obj)
    {
        int i;

        // Compare the elements in the subsetNames vector.
        boolean subsetNames_equal = (obj.subsetNames.size() == subsetNames.size());
        for(i = 0; (i < subsetNames.size()) && subsetNames_equal; ++i)
        {
            // Make references to String from Object.
            String subsetNames1 = (String)subsetNames.elementAt(i);
            String subsetNames2 = (String)obj.subsetNames.elementAt(i);
            subsetNames_equal = subsetNames1.equals(subsetNames2);
        }
        // Create the return value
        return ((colorType == obj.colorType) &&
                (colorTableName.equals(obj.colorTableName)) &&
                (invertColorTable == obj.invertColorTable) &&
                (filledFlag == obj.filledFlag) &&
                (legendFlag == obj.legendFlag) &&
                (lineStyle == obj.lineStyle) &&
                (lineWidth == obj.lineWidth) &&
                (singleColor == obj.singleColor) &&
                (multiColor.equals(obj.multiColor)) &&
                subsetNames_equal &&
                (subsetType == obj.subsetType) &&
                (opacity == obj.opacity) &&
                (wireframe == obj.wireframe) &&
                (drawInternal == obj.drawInternal) &&
                (smoothingLevel == obj.smoothingLevel) &&
                (pointSize == obj.pointSize) &&
                (pointType == obj.pointType) &&
                (pointSizeVarEnabled == obj.pointSizeVarEnabled) &&
                (pointSizeVar.equals(obj.pointSizeVar)) &&
                (pointSizePixels == obj.pointSizePixels));
    }

    public String GetName() { return "Subset"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetColorType(int colorType_)
    {
        colorType = colorType_;
        Select(0);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(1);
    }

    public void SetInvertColorTable(boolean invertColorTable_)
    {
        invertColorTable = invertColorTable_;
        Select(2);
    }

    public void SetFilledFlag(boolean filledFlag_)
    {
        filledFlag = filledFlag_;
        Select(3);
    }

    public void SetLegendFlag(boolean legendFlag_)
    {
        legendFlag = legendFlag_;
        Select(4);
    }

    public void SetLineStyle(int lineStyle_)
    {
        lineStyle = lineStyle_;
        Select(5);
    }

    public void SetLineWidth(int lineWidth_)
    {
        lineWidth = lineWidth_;
        Select(6);
    }

    public void SetSingleColor(ColorAttribute singleColor_)
    {
        singleColor = singleColor_;
        Select(7);
    }

    public void SetMultiColor(ColorAttributeList multiColor_)
    {
        multiColor = multiColor_;
        Select(8);
    }

    public void SetSubsetNames(Vector subsetNames_)
    {
        subsetNames = subsetNames_;
        Select(9);
    }

    public void SetSubsetType(int subsetType_)
    {
        subsetType = subsetType_;
        Select(10);
    }

    public void SetOpacity(double opacity_)
    {
        opacity = opacity_;
        Select(11);
    }

    public void SetWireframe(boolean wireframe_)
    {
        wireframe = wireframe_;
        Select(12);
    }

    public void SetDrawInternal(boolean drawInternal_)
    {
        drawInternal = drawInternal_;
        Select(13);
    }

    public void SetSmoothingLevel(int smoothingLevel_)
    {
        smoothingLevel = smoothingLevel_;
        Select(14);
    }

    public void SetPointSize(double pointSize_)
    {
        pointSize = pointSize_;
        Select(15);
    }

    public void SetPointType(int pointType_)
    {
        pointType = pointType_;
        Select(16);
    }

    public void SetPointSizeVarEnabled(boolean pointSizeVarEnabled_)
    {
        pointSizeVarEnabled = pointSizeVarEnabled_;
        Select(17);
    }

    public void SetPointSizeVar(String pointSizeVar_)
    {
        pointSizeVar = pointSizeVar_;
        Select(18);
    }

    public void SetPointSizePixels(int pointSizePixels_)
    {
        pointSizePixels = pointSizePixels_;
        Select(19);
    }

    // Property getting methods
    public int                GetColorType() { return colorType; }
    public String             GetColorTableName() { return colorTableName; }
    public boolean            GetInvertColorTable() { return invertColorTable; }
    public boolean            GetFilledFlag() { return filledFlag; }
    public boolean            GetLegendFlag() { return legendFlag; }
    public int                GetLineStyle() { return lineStyle; }
    public int                GetLineWidth() { return lineWidth; }
    public ColorAttribute     GetSingleColor() { return singleColor; }
    public ColorAttributeList GetMultiColor() { return multiColor; }
    public Vector             GetSubsetNames() { return subsetNames; }
    public int                GetSubsetType() { return subsetType; }
    public double             GetOpacity() { return opacity; }
    public boolean            GetWireframe() { return wireframe; }
    public boolean            GetDrawInternal() { return drawInternal; }
    public int                GetSmoothingLevel() { return smoothingLevel; }
    public double             GetPointSize() { return pointSize; }
    public int                GetPointType() { return pointType; }
    public boolean            GetPointSizeVarEnabled() { return pointSizeVarEnabled; }
    public String             GetPointSizeVar() { return pointSizeVar; }
    public int                GetPointSizePixels() { return pointSizePixels; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(colorType);
        if(WriteSelect(1, buf))
            buf.WriteString(colorTableName);
        if(WriteSelect(2, buf))
            buf.WriteBool(invertColorTable);
        if(WriteSelect(3, buf))
            buf.WriteBool(filledFlag);
        if(WriteSelect(4, buf))
            buf.WriteBool(legendFlag);
        if(WriteSelect(5, buf))
            buf.WriteInt(lineStyle);
        if(WriteSelect(6, buf))
            buf.WriteInt(lineWidth);
        if(WriteSelect(7, buf))
            singleColor.Write(buf);
        if(WriteSelect(8, buf))
            multiColor.Write(buf);
        if(WriteSelect(9, buf))
            buf.WriteStringVector(subsetNames);
        if(WriteSelect(10, buf))
            buf.WriteInt(subsetType);
        if(WriteSelect(11, buf))
            buf.WriteDouble(opacity);
        if(WriteSelect(12, buf))
            buf.WriteBool(wireframe);
        if(WriteSelect(13, buf))
            buf.WriteBool(drawInternal);
        if(WriteSelect(14, buf))
            buf.WriteInt(smoothingLevel);
        if(WriteSelect(15, buf))
            buf.WriteDouble(pointSize);
        if(WriteSelect(16, buf))
            buf.WriteInt(pointType);
        if(WriteSelect(17, buf))
            buf.WriteBool(pointSizeVarEnabled);
        if(WriteSelect(18, buf))
            buf.WriteString(pointSizeVar);
        if(WriteSelect(19, buf))
            buf.WriteInt(pointSizePixels);
    }

    public void ReadAtts(int index, CommunicationBuffer buf)
    {
        switch(index)
        {
        case 0:
            SetColorType(buf.ReadInt());
            break;
        case 1:
            SetColorTableName(buf.ReadString());
            break;
        case 2:
            SetInvertColorTable(buf.ReadBool());
            break;
        case 3:
            SetFilledFlag(buf.ReadBool());
            break;
        case 4:
            SetLegendFlag(buf.ReadBool());
            break;
        case 5:
            SetLineStyle(buf.ReadInt());
            break;
        case 6:
            SetLineWidth(buf.ReadInt());
            break;
        case 7:
            singleColor.Read(buf);
            Select(7);
            break;
        case 8:
            multiColor.Read(buf);
            Select(8);
            break;
        case 9:
            SetSubsetNames(buf.ReadStringVector());
            break;
        case 10:
            SetSubsetType(buf.ReadInt());
            break;
        case 11:
            SetOpacity(buf.ReadDouble());
            break;
        case 12:
            SetWireframe(buf.ReadBool());
            break;
        case 13:
            SetDrawInternal(buf.ReadBool());
            break;
        case 14:
            SetSmoothingLevel(buf.ReadInt());
            break;
        case 15:
            SetPointSize(buf.ReadDouble());
            break;
        case 16:
            SetPointType(buf.ReadInt());
            break;
        case 17:
            SetPointSizeVarEnabled(buf.ReadBool());
            break;
        case 18:
            SetPointSizeVar(buf.ReadString());
            break;
        case 19:
            SetPointSizePixels(buf.ReadInt());
            break;
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "colorType = ";
        if(colorType == COLORINGMETHOD_COLORBYSINGLECOLOR)
            str = str + "COLORINGMETHOD_COLORBYSINGLECOLOR";
        if(colorType == COLORINGMETHOD_COLORBYMULTIPLECOLORS)
            str = str + "COLORINGMETHOD_COLORBYMULTIPLECOLORS";
        if(colorType == COLORINGMETHOD_COLORBYCOLORTABLE)
            str = str + "COLORINGMETHOD_COLORBYCOLORTABLE";
        str = str + "\n";
        str = str + stringToString("colorTableName", colorTableName, indent) + "\n";
        str = str + boolToString("invertColorTable", invertColorTable, indent) + "\n";
        str = str + boolToString("filledFlag", filledFlag, indent) + "\n";
        str = str + boolToString("legendFlag", legendFlag, indent) + "\n";
        str = str + intToString("lineStyle", lineStyle, indent) + "\n";
        str = str + intToString("lineWidth", lineWidth, indent) + "\n";
        str = str + indent + "singleColor = {" + singleColor.Red() + ", " + singleColor.Green() + ", " + singleColor.Blue() + ", " + singleColor.Alpha() + "}\n";
        str = str + indent + "multiColor = {\n" + multiColor.toString(indent + "    ") + indent + "}\n";
        str = str + stringVectorToString("subsetNames", subsetNames, indent) + "\n";
        str = str + indent + "subsetType = ";
        if(subsetType == SUBSET_TYPE_DOMAIN)
            str = str + "SUBSET_TYPE_DOMAIN";
        if(subsetType == SUBSET_TYPE_GROUP)
            str = str + "SUBSET_TYPE_GROUP";
        if(subsetType == SUBSET_TYPE_MATERIAL)
            str = str + "SUBSET_TYPE_MATERIAL";
        if(subsetType == SUBSET_TYPE_ENUMSCALAR)
            str = str + "SUBSET_TYPE_ENUMSCALAR";
        if(subsetType == SUBSET_TYPE_MESH)
            str = str + "SUBSET_TYPE_MESH";
        if(subsetType == SUBSET_TYPE_UNKNOWN)
            str = str + "SUBSET_TYPE_UNKNOWN";
        str = str + "\n";
        str = str + doubleToString("opacity", opacity, indent) + "\n";
        str = str + boolToString("wireframe", wireframe, indent) + "\n";
        str = str + boolToString("drawInternal", drawInternal, indent) + "\n";
        str = str + intToString("smoothingLevel", smoothingLevel, indent) + "\n";
        str = str + doubleToString("pointSize", pointSize, indent) + "\n";
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
        str = str + boolToString("pointSizeVarEnabled", pointSizeVarEnabled, indent) + "\n";
        str = str + stringToString("pointSizeVar", pointSizeVar, indent) + "\n";
        str = str + intToString("pointSizePixels", pointSizePixels, indent) + "\n";
        return str;
    }


    // Attributes
    private int                colorType;
    private String             colorTableName;
    private boolean            invertColorTable;
    private boolean            filledFlag;
    private boolean            legendFlag;
    private int                lineStyle;
    private int                lineWidth;
    private ColorAttribute     singleColor;
    private ColorAttributeList multiColor;
    private Vector             subsetNames; // vector of String objects
    private int                subsetType;
    private double             opacity;
    private boolean            wireframe;
    private boolean            drawInternal;
    private int                smoothingLevel;
    private double             pointSize;
    private int                pointType;
    private boolean            pointSizeVarEnabled;
    private String             pointSizeVar;
    private int                pointSizePixels;
}

